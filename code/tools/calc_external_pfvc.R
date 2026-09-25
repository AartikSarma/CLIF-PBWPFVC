#!/usr/bin/env Rscript
# =============================================================================
# calc_external_pfvc.R -- GLI-2012 PFVC, PFVC at age 25 and Devine PBW for an external trial
#                         table (e.g. ARMA), by study arm, on the same scale as script 03.
# =============================================================================
# Given a trial's individual data (age, sex, race, height, study arm), computes each patient's
# predicted FVC (race-specific GLI-2012 at the actual age, pfvc), predicted FVC at age 25
# (pfvc_age25: GLI's height, sex and race scaling without its age decline) and Devine PBW, and
# reports the mean and 25th/50th/75th percentiles per arm. With a delivered tidal volume it also
# reports VT/PBW, VT/PFVC and VT/PFVC at age 25 per arm. The computation mirrors script 03
# (rspiro::pred_GLI; race -> GLI ethnicity 1/2/5; height in cm; valid ages 3-95), so the values
# are directly comparable to the cohort's pfvc, pfvc_age25, pbw, vtpbw and vtpfvc.
#
# Input CSV columns (case-insensitive; common synonyms auto-detected):
#   age    years
#   sex    Male/Female | M/F | 1=male,2=female
#   race   White/Black/Other   (GLI: White->1 Caucasian, Black->2 African-American, else->5 Other)
#   height cm   (inches / metres auto-detected & converted, with a printed note)
#   arm    any label (e.g. "6 mL/kg", "12 mL/kg")
#   tidal_volume  OPTIONAL delivered VT (mL or L auto-detected). If present, also reports the
#                 delivered strain per arm: VT/PBW (mL/kg), VT/PFVC (%), VT/PFVC at age 25 (%).
# Usage:  uvr run code/tools/calc_external_pfvc.R <input.csv> [output.csv]
#
# Every assumption (column match, race mapping, height and VT units, dropped rows) is printed for audit.
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(rspiro) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: uvr run code/tools/calc_external_pfvc.R <input.csv> [output.csv]")
in_path  <- args[1]
out_path <- if (length(args) >= 2) args[2] else sub("\\.csv$", "_pfvc_by_arm.csv", in_path)
raw <- readr::read_csv(in_path, show_col_types = FALSE)
cat("Read", nrow(raw), "rows,", ncol(raw), "cols from", in_path, "\n")

pick <- function(cands) {
  hit <- names(raw)[tolower(trimws(names(raw))) %in% cands]
  if (length(hit) == 0) stop("No column for {", paste(cands, collapse = "/"),
                             "}. Available: ", paste(names(raw), collapse = ", "))
  hit[1]
}
c_age <- pick(c("age", "age_at_admission", "age_years", "ageyrs"))
c_sex <- pick(c("sex", "sex_category", "gender"))
c_race<- pick(c("race", "race_category", "ethnicity", "race_ethnicity"))
c_ht  <- pick(c("height", "height_cm", "ht", "ht_cm", "heightcm", "height_in"))
c_arm <- pick(c("arm", "study_arm", "group", "treatment", "trt", "randomization", "tidal_volume_group"))
pick_opt <- function(cands) { hit <- names(raw)[tolower(trimws(names(raw))) %in% cands]
  if (length(hit)) hit[1] else NA_character_ }
c_vt <- pick_opt(c("tidal_volume", "tidal_volume_set", "tidal_volume_ml", "vt", "vt_ml",
                   "vt_set", "vt_delivered", "delivered_vt", "vtdelivered", "tidalvolume"))
cat(sprintf("Columns -> age:%s  sex:%s  race:%s  height:%s  arm:%s  tidal_volume:%s\n",
            c_age, c_sex, c_race, c_ht, c_arm, if (is.na(c_vt)) "(none)" else c_vt))

df <- tibble(age = as.numeric(raw[[c_age]]), sex_raw = as.character(raw[[c_sex]]),
             race_raw = as.character(raw[[c_race]]), height = as.numeric(raw[[c_ht]]),
             arm = as.character(raw[[c_arm]]),
             vt = if (!is.na(c_vt)) as.numeric(raw[[c_vt]]) else NA_real_)

# sex -> GLI gender (1 male, 2 female); error (do not guess) on unmapped values
s <- tolower(trimws(df$sex_raw))
df$gender <- case_when(s %in% c("male", "m", "1", "man") ~ 1L,
                       s %in% c("female", "f", "2", "woman") ~ 2L, TRUE ~ NA_integer_)
if (any(is.na(df$gender) & !is.na(df$sex_raw) & df$sex_raw != ""))
  stop("Unmapped sex values: ", paste(unique(df$sex_raw[is.na(df$gender)]), collapse = ", "))

# race -> GLI ethnicity (1 Caucasian, 2 African-American, 5 Other), the script-03 convention
r <- tolower(trimws(df$race_raw))
df$ethnicity <- case_when(
  r %in% c("white", "caucasian", "european", "1") ~ 1L,
  r %in% c("black", "african american", "african-american", "aa", "2") ~ 2L, TRUE ~ 5L)
cat("\nRace -> GLI ethnicity mapping (verify this is correct for your coding):\n")
print(as.data.frame(df %>% count(race_raw, ethnicity) %>% arrange(ethnicity)), row.names = FALSE)

# height -> cm (auto-detect unit; error rather than silently mis-scale)
med_h <- median(df$height, na.rm = TRUE)
unit <- case_when(med_h >= 120 & med_h <= 230 ~ "cm", med_h >= 1.2 & med_h <= 2.4 ~ "m",
                  med_h >= 48 & med_h <= 90 ~ "in", TRUE ~ "unknown")
if (unit == "unknown") stop("Cannot infer height unit from median ", round(med_h, 1),
                            " -- expected cm (~170), m (~1.7), or inches (~67). Pass height in cm.")
df$height_cm <- df$height * c(cm = 1, m = 100, `in` = 2.54)[[unit]]
cat(sprintf("\nHeight unit detected: %s (median %.1f -> %.1f cm)\n", unit, med_h, median(df$height_cm, na.rm = TRUE)))

# drop GLI/Devine-uncomputable rows, report count + reason
n0 <- nrow(df)
keep <- !(is.na(df$age) | df$age < 3 | df$age > 95 | is.na(df$gender) |
            is.na(df$height_cm) | df$height_cm < 120 | df$height_cm > 230 | is.na(df$arm) | df$arm == "")
if (sum(!keep) > 0)
  cat(sprintf("Dropping %d/%d rows (missing/out-of-range: age 3-95, sex, height 120-230cm, or arm)\n", sum(!keep), n0))
df <- df[keep, ]
stopifnot(nrow(df) > 0)

# GLI PFVC (actual age) and pfvc_age25 (age 25), Devine PBW -- identical to script 03
df <- df %>% mutate(
  pfvc       = pred_GLI(age = age,          height = height_cm / 100, gender = gender, ethnicity = ethnicity, param = "FVC"),
  pfvc_age25 = pred_GLI(age = rep(25, n()), height = height_cm / 100, gender = gender, ethnicity = ethnicity, param = "FVC"),
  pbw        = if_else(gender == 1L, 50 + 2.3 * (height_cm / 2.54 - 60), 45.5 + 2.3 * (height_cm / 2.54 - 60)))

# delivered strain on each normalizer scale, IF a tidal-volume column was supplied
metrics <- c("pfvc", "pfvc_age25", "pbw")
if (!is.na(c_vt) && any(is.finite(df$vt))) {
  med_v <- median(df$vt, na.rm = TRUE)
  vunit <- if (med_v >= 100 & med_v <= 2000) "mL" else if (med_v >= 0.1 & med_v <= 2) "L" else "unknown"
  if (vunit == "unknown") stop("Cannot infer tidal-volume unit from median ", round(med_v, 1),
                               " -- expected mL (~450) or L (~0.45).")
  df$vt_ml <- if (vunit == "L") df$vt * 1000 else df$vt
  cat(sprintf("Tidal volume unit detected: %s (median %.1f -> %.0f mL)\n", vunit, med_v, median(df$vt_ml, na.rm = TRUE)))
  df <- df %>% mutate(vtpbw        = vt_ml / pbw,                 # mL/kg PBW (sanity: ~6 / ~12 by arm)
                      vtpfvc       = vt_ml / pfvc * 0.1,          # % strain, PFVC scale (engine units)
                      vtpfvc_age25 = vt_ml / pfvc_age25 * 0.1)    # % strain, pfvc_age25 scale
  metrics <- c(metrics, "vtpbw", "vtpfvc", "vtpfvc_age25")
}

# per-arm (+ pooled) mean & quartiles for each metric (NA-robust; n = non-missing used)
summ_one <- function(d, arm_lab) bind_rows(lapply(metrics, function(m) {
  x <- d[[m]]; x <- x[is.finite(x)]
  tibble(arm = arm_lab, metric = m, n = length(x), mean = mean(x),
         p25 = quantile(x, .25, names = FALSE), p50 = median(x), p75 = quantile(x, .75, names = FALSE)) }))
parts <- c(unname(lapply(split(df, df$arm), function(d) summ_one(d, d$arm[1]))),
           list(summ_one(df, "All (pooled)")))   # unname: split() names would nest into columns
out <- bind_rows(parts) %>% mutate(across(c(mean, p25, p50, p75), ~round(., 3)))
write_csv(out, out_path)

cat("\n=== pfvc / pfvc_age25 / PBW by study arm (mean + quartiles, litres / kg) ===\n")
print(as.data.frame(out), row.names = FALSE)
cat("\nBy randomization pfvc / pfvc_age25 / PBW should be about equal across arms (balance check).\n")
if ("vtpfvc" %in% metrics) {
  cat("SANITY: vtpbw should be ~6 (low arm) and ~12 (high arm). vtpfvc / vtpfvc_age25 are the\n",
      "delivered strain on each scale, comparable to the cohort's vtpfvc (script 03).\n")
} else {
  cat("No tidal-volume column found -- send delivered VT to also get VT/PFVC and VT/PFVC at age 25 per arm.\n")
}
message("Wrote ", out_path)
