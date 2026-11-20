# Data-zone helpers for local CLIF projects (see SOP.md). R twin of clif_paths.py.
#
#   1. Core CLIF tables stay where they are: read them from data_directory, never copy them.
#   2. Intermediate files go to intermediate_directory, by default a folder OUTSIDE the repo.
#   3. Intermediate files carry no dates, geography or direct identifiers unless a column is
#      explicitly allowed (and the reason is written down in the README).
#   4. Final files go to results_directory; only reviewed aggregates are copied to
#      output/final_no_phi/ with release_result().
#
# Usage (run from the project root):
#   source(file.path("utils", "clif_paths.R"))
#   cfg <- load_config()
#
# Needs jsonlite and arrow.

PHI_NAME_PATTERNS <- c(
  "date or time"      = "dttm|(^|_)date($|_)|birth|(^|_)dob($|_)",
  "geography"         = "zip|census|county|fips|state_code|latitude|longitude|address|geocode",
  "direct identifier" = "mrn|ssn|first_name|last_name|patient_name|phone|email"
)
AGE_PATTERN <- "(^|_)age($|_)"
AGE_CAP <- 90  # HIPAA Safe Harbor: ages over 89 are reported as 90+
RELEASE_DIR <- file.path("output", "final_no_phi")

.norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

.is_within <- function(path, parent) {
  path <- .norm(path); parent <- .norm(parent)
  path == parent || startsWith(path, paste0(sub("/$", "", parent), "/"))
}

.expand <- function(value) {
  # Expand ~ and $VAR / ${VAR}
  value <- path.expand(value)
  for (v in regmatches(value, gregexpr("\\$\\{?[A-Za-z_][A-Za-z0-9_]*\\}?", value))[[1]]) {
    name <- gsub("[${}]", "", v)
    value <- sub(v, Sys.getenv(name), value, fixed = TRUE)
  }
  value
}

load_config <- function(path = file.path("config", "config.json")) {
  if (!file.exists(path)) {
    stop(path, " not found: copy config/config_template.json to config/config.json")
  }
  cfg <- jsonlite::fromJSON(path)
  root <- .norm(getwd())
  workspace <- file.path(dirname(root), "clif_workspace", basename(root))
  data_dir <- .norm(.expand(cfg$data_directory))
  if (!dir.exists(data_dir)) stop("data_directory does not exist: ", data_dir)
  cfg$data_directory <- data_dir

  defaults <- c(output_directory = "logs", intermediate_directory = "intermediate",
                results_directory = "results")
  for (key in names(defaults)) {
    value <- cfg[[key]]
    folder <- if (is.null(value) || !nzchar(value)) file.path(workspace, defaults[[key]]) else .expand(value)
    folder <- .norm(folder)
    if (.is_within(folder, root)) {
      stop(key, " (", folder, ") is inside the project folder; point it outside the repo")
    }
    if (.is_within(folder, data_dir)) {
      stop(key, " (", folder, ") is inside data_directory; core tables are read-only")
    }
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    cfg[[key]] <- .norm(folder)
  }
  if (is.null(cfg$min_cell_size)) cfg$min_cell_size <- 10
  if (is.null(cfg$data_release)) cfg$data_release <- "unspecified"
  cfg
}

# Return a named character vector: column -> reason, for columns that look like PHI.
phi_columns <- function(df, allow = character()) {
  found <- character()
  for (col in setdiff(names(df), allow)) {
    x <- df[[col]]
    if (inherits(x, c("POSIXt", "Date"))) {
      found[col] <- "date or time type"
      next
    }
    hit <- names(PHI_NAME_PATTERNS)[vapply(PHI_NAME_PATTERNS, grepl, logical(1), x = col,
                                           ignore.case = TRUE)]
    if (length(hit)) {
      found[col] <- hit[1]
    } else if (grepl(AGE_PATTERN, col, ignore.case = TRUE) && is.numeric(x) &&
               isTRUE(max(x, na.rm = TRUE) > AGE_CAP)) {
      found[col] <- sprintf("age above %d; top-code it with top_code_age()", AGE_CAP)
    }
  }
  found
}

# Replace every *_dttm column with time since index_col (units: "hours", "days", "mins").
# New columns are named <col>_offset_h / _offset_d / _offset_min; originals are dropped.
to_offsets <- function(df, index_col, units = "hours", drop_index = TRUE) {
  suffix <- c(hours = "h", days = "d", mins = "min")[[units]]
  index <- df[[index_col]]
  for (col in setdiff(grep("_dttm$", names(df), value = TRUE), index_col)) {
    new <- sub("_dttm$", paste0("_offset_", suffix), col)
    df[[new]] <- as.numeric(difftime(df[[col]], index, units = units))
    df[[col]] <- NULL
  }
  if (drop_index) df[[index_col]] <- NULL
  df
}

top_code_age <- function(df, col = "age_at_admission", cap = AGE_CAP) {
  df[[col]] <- pmin(df[[col]], cap)
  df
}

write_intermediate <- function(df, name, cfg, allow = character()) {
  flagged <- phi_columns(df, allow)
  if (length(flagged)) {
    stop("Not writing ", name, ": possible PHI in ",
         paste0(names(flagged), " (", flagged, ")", collapse = ", "),
         ". Reduce it in memory (to_offsets, drop geography) or pass allow = c(...) ",
         "and record the reason in the README.")
  }
  path <- file.path(cfg$intermediate_directory, paste0(name, ".parquet"))
  arrow::write_parquet(df, path)
  invisible(path)
}

write_result <- function(df, name, cfg) {
  path <- file.path(cfg$results_directory, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

small_cells <- function(csv_path, min_n) {
  df <- utils::read.csv(csv_path, check.names = FALSE)
  bad <- character()
  for (col in grep("^(n|count|n_.*|.*_n|.*_count)$", names(df), value = TRUE, ignore.case = TRUE)) {
    vals <- suppressWarnings(as.numeric(df[[col]]))
    rows <- which(vals > 0 & vals < min_n)
    bad <- c(bad, sprintf("%s row %d: %g", col, rows, vals[rows]))
  }
  bad
}

# Copy a reviewed result to output/final_no_phi/<site>_<name>. A person must still review it.
release_result <- function(name, cfg) {
  src <- file.path(cfg$results_directory, name)
  if (!file.exists(src)) stop("not found: ", src)
  if (grepl("\\.parquet$", name)) stop("Parquet files are not released; release aggregate CSVs or figures")
  if (grepl("\\.csv$", name)) {
    bad <- small_cells(src, cfg$min_cell_size)
    if (length(bad)) stop(name, " has cells below n=", cfg$min_cell_size, ": ",
                          paste(utils::head(bad, 10), collapse = "; "))
  }
  dir.create(RELEASE_DIR, recursive = TRUE, showWarnings = FALSE)
  dst <- file.path(RELEASE_DIR, paste0(cfg$site_name, "_", basename(src)))
  file.copy(src, dst, overwrite = TRUE)
  invisible(dst)
}

write_manifest <- function(cfg, counts = list()) {
  commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                     error = function(e) "unknown", warning = function(w) "unknown")
  dirty <- tryCatch(length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0,
                    error = function(e) NA)
  manifest <- list(
    run_date = format(Sys.Date()),
    git_commit = commit[1],
    uncommitted_changes = dirty,
    site_name = cfg$site_name,
    clif_version = cfg$clif_version,
    data_release = cfg$data_release,
    r_version = R.version.string,
    counts = counts
  )
  path <- file.path(cfg$results_directory, "run_manifest.json")
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}
