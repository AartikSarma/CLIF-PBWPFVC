#' Clean + waterfall-fill the CLIF resp_support table
#'
#' R port of the Python reference pipeline. Builds an hourly scaffold,
#' applies device/mode heuristics, creates hierarchical episode IDs,
#' and performs numeric waterfall fill inside each mode_name_id block.
#'
#' @param resp_support A data.frame / tibble of the raw CLIF respiratory-support
#'   table (timestamps assumed UTC).
#' @param id_col Character. Encounter-level identifier column.
#'   Default \code{"hospitalization_id"}.
#' @param bfill Logical. If TRUE numeric setters are back-filled after
#'   forward-fill; if FALSE (default) only forward-fill is used.
#' @param verbose Logical. Print progress banners when TRUE.
#'
#' @return A tibble with hourly scaffold rows inserted, device/mode heuristics
#'   applied, hierarchical episode IDs, numeric waterfall fill, tracheostomy
#'   flag forward-filled, and one unique row per (id_col, recorded_dttm).

library(tidyverse)
library(data.table)
library(lubridate)

process_resp_support_waterfall <- function(resp_support,
                                           id_col = "hospitalization_id",
                                           bfill = FALSE,
                                           verbose = TRUE) {

  p <- if (verbose) message else function(...) invisible(NULL)

  # Helper: forward-fill (+ optional back-fill) a vector
  fb <- function(x) {
    out <- vctrs::vec_fill_missing(x, direction = "down")
    if (bfill) out <- vctrs::vec_fill_missing(out, direction = "up")
    out
  }

  # Helper: change-ID — cumulative sum of transitions within groups
  change_id <- function(col, id) {
    filled <- replace_na(col, "missing")
    tibble(.id = id, .val = filled) |>
      mutate(.changed = .val != lag(.val, default = ""),
             .changed = replace_na(.changed, TRUE),
             .block = cumsum(.changed),
             .data.table.aware = NULL,
             .id_change = .val != lag(.val, default = ""),
             .id_change = replace_na(.id_change, TRUE)) |>
      pull(.block) -> blocks
    # Redo properly: per-group cumsum
    dt <- data.table(.id = id, .val = replace_na(col, "missing"))
    dt[, .result := cumsum(c(TRUE, .val[-1] != .val[-.N])), by = .id]
    as.integer(dt$.result)
  }

  # ------------------------------------------------------------------ #
  # Phase 0 - set-up & hourly scaffold                                 #
  # ------------------------------------------------------------------ #
  p("Phase 0: initialise & create hourly scaffold")
  rs <- as_tibble(resp_support)

  # Lower-case categorical strings
  cat_cols <- intersect(
    c("device_category", "device_name", "mode_category", "mode_name"),
    names(rs)
  )
  rs <- rs |> mutate(across(all_of(cat_cols), tolower))

  # Numeric coercion
  num_cols <- intersect(
    c("tracheostomy", "fio2_set", "lpm_set", "peep_set",
      "tidal_volume_set", "resp_rate_set", "resp_rate_obs",
      "pressure_support_set", "peak_inspiratory_pressure_set"),
    names(rs)
  )
  rs <- rs |> mutate(across(all_of(num_cols), \(x) suppressWarnings(as.numeric(x))))

  # FiO2 scaling: if documented as 40 -> 0.40

  if ("fio2_set" %in% names(rs)) {
    fio2_mean <- mean(rs$fio2_set, na.rm = TRUE)
    if (!is.na(fio2_mean) && fio2_mean > 1.0) {
      rs <- rs |> mutate(fio2_set = if_else(fio2_set > 1, fio2_set / 100, fio2_set))
      p("  Scaled FiO2 values > 1 down by /100")
    }
  }

  # Build hourly scaffold
  p("  Building hourly scaffold")
  bounds <- rs |>
    filter(!is.na(recorded_dttm)) |>
    group_by(.data[[id_col]]) |>
    summarise(
      tmin_h = floor_date(min(recorded_dttm), "hour"),
      tmax_h = floor_date(max(recorded_dttm), "hour"),
      .groups = "drop"
    )

  scaffold <- bounds |>
    rowwise() |>
    reframe(
      !!id_col := .data[[id_col]],
      recorded_dttm = seq(tmin_h, tmax_h, by = "1 hour") + seconds(59 * 60 + 59)
    ) |>
    mutate(
      recorded_date = as.Date(recorded_dttm),
      recorded_hour = hour(recorded_dttm),
      is_scaffold = TRUE
    )

  if (verbose) p(sprintf("  Scaffold rows created: %s", format(nrow(scaffold), big.mark = ",")))

  rs <- rs |>
    mutate(
      recorded_date = as.Date(recorded_dttm),
      recorded_hour = hour(recorded_dttm)
    )

  # ------------------------------------------------------------------ #
  # Phase 1 - heuristic device / mode inference                        #
  # ------------------------------------------------------------------ #
  p("Phase 1: heuristic inference of device & mode")

  # Most-frequent fall-back labels
  device_counts <- rs |>
    filter(!is.na(device_name), !is.na(device_category)) |>
    count(device_name, device_category, sort = TRUE)

  most_common_imv_name <- device_counts |>
    filter(device_category == "imv") |>
    slice(1) |>
    pull(device_name)
  most_common_imv_name <- if (length(most_common_imv_name) == 0) "ventilator" else most_common_imv_name

  most_common_nippv_name <- device_counts |>
    filter(device_category == "nippv") |>
    slice(1) |>
    pull(device_name)
  most_common_nippv_name <- if (length(most_common_nippv_name) == 0) "bipap" else most_common_nippv_name

  mode_counts <- rs |>
    filter(!is.na(mode_name), !is.na(mode_category)) |>
    count(mode_name, mode_category, sort = TRUE)

  most_common_cmv_name <- mode_counts |>
    filter(mode_category == "assist control-volume control") |>
    slice(1) |>
    pull(mode_name)
  most_common_cmv_name <- if (length(most_common_cmv_name) == 0) "AC/VC" else most_common_cmv_name

  # 1-a: IMV from mode_category
  imv_mode_pattern <- "assist control-volume control|simv|pressure control"
  mask_1a <- is.na(rs$device_category) &
    is.na(rs$device_name) &
    str_detect(replace_na(rs$mode_category, ""), imv_mode_pattern)
  rs$device_category[mask_1a] <- "imv"
  rs$device_name[mask_1a] <- most_common_imv_name

  # 1-b: IMV look-behind/ahead
  rs <- rs |> arrange(.data[[id_col]], recorded_dttm)
  rs <- rs |>
    group_by(.data[[id_col]]) |>
    mutate(
      prev_cat = lag(device_category),
      next_cat = lead(device_category)
    ) |>
    ungroup()

  imv_like <- is.na(rs$device_category) &
    (replace_na(rs$prev_cat, "") == "imv" | replace_na(rs$next_cat, "") == "imv") &
    replace_na(rs$peep_set > 1, FALSE) &
    replace_na(rs$resp_rate_set > 1, FALSE) &
    replace_na(rs$tidal_volume_set > 1, FALSE)
  rs$device_category[imv_like] <- "imv"
  rs$device_name[imv_like] <- most_common_imv_name

  # 1-c: NIPPV heuristics (recompute prev/next after 1-b changes)
  rs <- rs |>
    group_by(.data[[id_col]]) |>
    mutate(
      prev_cat = lag(device_category),
      next_cat = lead(device_category)
    ) |>
    ungroup()

  nippv_like <- is.na(rs$device_category) &
    (replace_na(rs$prev_cat, "") == "nippv" | replace_na(rs$next_cat, "") == "nippv") &
    replace_na(rs$peak_inspiratory_pressure_set > 1, FALSE) &
    replace_na(rs$pressure_support_set > 1, FALSE)
  rs$device_category[nippv_like] <- "nippv"
  rs$device_name[nippv_like & is.na(rs$device_name)] <- most_common_nippv_name

  rs <- rs |> select(-prev_cat, -next_cat)

  # 1-d: Clean duplicates & empty rows
  rs <- rs |> arrange(.data[[id_col]], recorded_dttm)
  rs <- rs |>
    group_by(.data[[id_col]], recorded_dttm) |>
    mutate(dup_count = n()) |>
    ungroup() |>
    filter(!(dup_count > 1 & device_category == "nippv"))

  rs <- rs |>
    group_by(.data[[id_col]], recorded_dttm) |>
    mutate(dup_count = n()) |>
    ungroup() |>
    filter(!(dup_count > 1 & is.na(device_category))) |>
    select(-dup_count)

  # 1-e: Nasal-cannula rows must never carry PEEP
  if ("peep_set" %in% names(rs)) {
    mask_bad_nc <- replace_na(rs$device_category == "nasal cannula" & rs$peep_set > 0, FALSE)
    if (any(mask_bad_nc)) {
      rs$device_category[mask_bad_nc] <- NA_character_
      p(sprintf("  %s rows had PEEP>0 on nasal cannula — device_category reset",
                format(sum(mask_bad_nc), big.mark = ",")))
    }
  }

  # Drop rows with nothing useful
  all_na_cols <- intersect(
    c("device_category", "device_name", "mode_category", "mode_name",
      "tracheostomy", "fio2_set", "lpm_set", "peep_set", "tidal_volume_set",
      "resp_rate_set", "resp_rate_obs", "pressure_support_set",
      "peak_inspiratory_pressure_set"),
    names(rs)
  )
  rs <- rs |> filter(!if_all(all_of(all_na_cols), is.na))

  # Unique per timestamp
  rs <- rs |> distinct(.data[[id_col]], recorded_dttm, .keep_all = TRUE)

  # Merge scaffold
  rs <- rs |> mutate(is_scaffold = FALSE)
  rs <- bind_rows(rs, scaffold) |>
    arrange(.data[[id_col]], recorded_dttm, recorded_date, recorded_hour)

  # ------------------------------------------------------------------ #
  # Phase 2 - hierarchical IDs                                         #
  # ------------------------------------------------------------------ #
  p("Phase 2: build hierarchical IDs")

  # Forward-fill device_category per encounter
  rs <- rs |>
    arrange(.data[[id_col]], recorded_dttm) |>
    group_by(.data[[id_col]]) |>
    mutate(device_category = vctrs::vec_fill_missing(device_category, direction = "down")) |>
    ungroup()

  rs$device_cat_id <- change_id(rs$device_category, rs[[id_col]])

  # Forward(/back)-fill device_name within (id, device_cat_id)
  rs <- rs |>
    arrange(recorded_dttm) |>
    group_by(.data[[id_col]], device_cat_id) |>
    mutate(device_name = fb(device_name)) |>
    ungroup()

  rs$device_id <- change_id(rs$device_name, rs[[id_col]])

  # Forward(/back)-fill mode_category within (id, device_id)
  rs <- rs |>
    arrange(.data[[id_col]], recorded_dttm) |>
    group_by(.data[[id_col]], device_id) |>
    mutate(mode_category = fb(mode_category)) |>
    ungroup()

  rs$mode_cat_id <- change_id(replace_na(rs$mode_category, "missing"), rs[[id_col]])

  # Forward(/back)-fill mode_name within (id, mode_cat_id)
  rs <- rs |>
    group_by(.data[[id_col]], mode_cat_id) |>
    mutate(mode_name = fb(mode_name)) |>
    ungroup()

  rs$mode_name_id <- change_id(replace_na(rs$mode_name, "missing"), rs[[id_col]])

  # ------------------------------------------------------------------ #
  # Phase 3 - numeric waterfall                                        #
  # ------------------------------------------------------------------ #
  fill_type <- if (bfill) "bi-directional" else "forward-only"
  p(sprintf("Phase 3: %s numeric fill inside mode_name_id blocks", fill_type))

  # FiO2 default for room-air
  if ("fio2_set" %in% names(rs)) {
    rs <- rs |>
      mutate(fio2_set = if_else(
        replace_na(device_category == "room air", FALSE) & is.na(fio2_set),
        0.21, fio2_set
      ))
  }

  # Tidal-volume clean-up
  if ("tidal_volume_set" %in% names(rs)) {
    bad_tv <- (
      (replace_na(rs$mode_category == "pressure support/cpap", FALSE) &
         !is.na(rs$pressure_support_set)) |
      (is.na(rs$mode_category) &
         replace_na(str_detect(rs$device_name, "trach"), FALSE)) |
      (replace_na(rs$mode_category == "pressure support/cpap", FALSE) &
         replace_na(str_detect(rs$device_name, "trach"), FALSE))
    )
    rs$tidal_volume_set[bad_tv] <- NA_real_
  }

  num_cols_fill <- intersect(
    c("fio2_set", "lpm_set", "peep_set", "tidal_volume_set",
      "pressure_support_set", "resp_rate_set", "resp_rate_obs",
      "peak_inspiratory_pressure_set"),
    names(rs)
  )

  n_enc <- n_distinct(rs[[id_col]])
  p(sprintf("  applying waterfall fill to %s encounters", format(n_enc, big.mark = ",")))

  # Waterfall fill within mode_name_id blocks
  # Use data.table for performance
  dt <- as.data.table(rs)

  for (col_name in num_cols_fill) {
    # Within each (id, mode_name_id) block, handle trach collar sub-breaks
    dt[, (col_name) := {
      vals <- get(col_name)
      dc <- device_category
      if (any(dc == "trach collar", na.rm = TRUE)) {
        breaker <- cumsum(replace_na(dc == "trach collar", FALSE))
        unsplit(lapply(split(vals, breaker), fb), breaker)
      } else {
        fb(vals)
      }
    }, by = c(id_col, "mode_name_id")]
  }

  rs <- as_tibble(dt)

  # T-piece -> classify as blow by
  tpiece <- is.na(rs$mode_category) &
    replace_na(str_detect(rs$device_name, "t-piece"), FALSE)
  rs$mode_category[tpiece] <- "blow by"

  # Tracheostomy flag forward-fill per encounter
  if ("tracheostomy" %in% names(rs)) {
    rs <- rs |>
      group_by(.data[[id_col]]) |>
      mutate(tracheostomy = vctrs::vec_fill_missing(tracheostomy, direction = "down")) |>
      ungroup()
  }

  # ------------------------------------------------------------------ #
  # Phase 4 - final tidy-up                                            #
  # ------------------------------------------------------------------ #
  p("Phase 4: final dedup & ordering")
  rs <- rs |>
    distinct() |>
    arrange(.data[[id_col]], recorded_dttm) |>
    select(-any_of(c("recorded_date", "recorded_hour")))

  p("[OK] Respiratory-support waterfall complete.")
  rs
}
