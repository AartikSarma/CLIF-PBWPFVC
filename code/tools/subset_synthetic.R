# =============================================================================
# Script subset_synthetic: a smaller synthetic CLIF dataset for fast test loops
# =============================================================================
#
# Developer tool, not part of the pipeline. Samples a fraction of the patients in
# a synthetic CLIF folder and writes every table restricted to those patients to
# a sibling folder, so the pipeline can be exercised in a
# fraction of the time while code is still changing. Point the pipeline at the
# subset without editing config.json:
#
#   Rscript code/tools/subset_synthetic.R --frac 0.3 --seed 1
#   PBWPFVC_TABLES_PATH=~/Research/synthetic_clif/synth_clif_10k_sub30 \
#     Rscript code/01_cohort_identification.R      # and so on for 02, 03
#
# Tables keyed by hospitalization_id are filtered on the sampled hospitalizations,
# tables keyed by patient_id on the sampled patients, and tables with neither key
# (microbiology_susceptibility) are copied whole. Never run this on real data; it
# reads the folder named in config.json (or PBWPFVC_TABLES_PATH) and refuses any
# site whose name is not synthetic_clif.
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(here)
})
source("utils/config.R")

if (!identical(config$site_name, "synthetic_clif"))
  stop("subset_synthetic.R only subsets the synthetic site; config site_name is '",
       config$site_name, "'.")
if (!identical(config$file_type, "parquet"))
  stop("subset_synthetic.R handles parquet folders only.")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- which(args == flag)
  if (length(hit) == 1L && hit < length(args)) args[hit + 1L] else default
}
frac <- as.numeric(arg_value("--frac", "0.3"))
seed <- as.integer(arg_value("--seed", "1"))
stopifnot(is.finite(frac), frac > 0, frac < 1, is.finite(seed))

src <- path.expand(config$tables_path)
out <- arg_value("--out", paste0(sub("/+$", "", src), "_sub", round(100 * frac)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
message("Subsetting ", src, " -> ", out, " (fraction ", frac, ", seed ", seed, ")")

patients <- read_parquet(file.path(src, "clif_patient.parquet"))
set.seed(seed)
keep_patients <- sample(unique(patients$patient_id), size = ceiling(frac * n_distinct(patients$patient_id)))
hosp <- read_parquet(file.path(src, "clif_hospitalization.parquet")) %>%
  filter(patient_id %in% keep_patients)
keep_hosp <- unique(hosp$hospitalization_id)
message("  patients: ", length(keep_patients), " of ", n_distinct(patients$patient_id),
        "; hospitalizations: ", length(keep_hosp))

for (f in list.files(src, pattern = "^clif_.*\\.parquet$")) {
  ds   <- open_dataset(file.path(src, f))
  cols <- names(ds$schema)
  sub <- if ("hospitalization_id" %in% cols) {
    ds %>% filter(hospitalization_id %in% keep_hosp) %>% collect()
  } else if ("patient_id" %in% cols) {
    ds %>% filter(patient_id %in% keep_patients) %>% collect()
  } else {
    ds %>% collect()
  }
  write_parquet(sub, file.path(out, f))
  message(sprintf("  %-42s %9d rows", f, nrow(sub)))
}
message("Subset written to ", out)
