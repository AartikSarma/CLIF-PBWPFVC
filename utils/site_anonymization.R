# =============================================================================
# site_anonymization.R: size-ordered, de-identified cohort labels
# =============================================================================
# Every cross-cohort (pooled) figure and table labels the contributing cohorts
# "Site A", "Site B", ... ordered from the LARGEST analytic cohort to the
# smallest, so no shared output names a contributing institution.
#
# The size ranking comes from the per-site attrition logs (written by script 03),
# whose last step is the analytic cohort N. Using that single source keeps a given
# cohort's alias identical across every pooled script (code/pooling/), even though a
# pooled table may key sites by the `site` column written into each results table or
# by the results FOLDER name. The two can differ (a folder named "EU" may hold
# results written by site "emory"), so both spellings are registered as keys for the
# same alias.
#
# The real-name -> alias crosswalk is only ever printed to the console
# (print_site_alias_key); it is deliberately NOT written to the pooled output
# folder, which is the folder that gets shared.

library(tidyverse)

# Alias suffixes: A..Z, then AA, AB, ... for consortia larger than 26 cohorts.
site_alias_suffixes <- function(n_sites) {
  suffixes <- LETTERS
  if (n_sites > length(suffixes)) {
    suffixes <- c(suffixes, as.vector(t(outer(LETTERS, LETTERS, paste0))))
  }
  if (n_sites > length(suffixes)) {
    stop("More cohorts (", n_sites, ") than available site aliases (",
         length(suffixes), ").")
  }
  suffixes[seq_len(n_sites)]
}

# Analytic cohort size for one site folder, read from its attrition log. Accepts
# either a flat site folder or one whose tables sit in a nested `final/`.
read_site_cohort_size <- function(site_dir) {
  log_files <- c(Sys.glob(file.path(site_dir, "attrition_log_*.csv")),
                 Sys.glob(file.path(site_dir, "final", "attrition_log_*.csv")),
                 Sys.glob(file.path(site_dir, "cross_sectional", "attrition_log_*.csv")),
                 Sys.glob(file.path(site_dir, "final", "cross_sectional", "attrition_log_*.csv")))
  if (length(log_files) == 0) {
    stop("No attrition_log_*.csv in ", site_dir, ". The pooled figures need it to ",
         "rank cohorts by size for anonymized site labels; re-run scripts 01-03 for ",
         "that cohort and copy its results in.")
  }
  attrition <- read_csv(log_files[[1]], show_col_types = FALSE)
  missing_cols <- setdiff(c("site", "step_order", "n_remaining"), names(attrition))
  if (length(missing_cols) > 0) {
    stop(log_files[[1]], " is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  site_names <- unique(as.character(attrition$site))
  if (length(site_names) != 1) {
    stop(log_files[[1]], " covers ", length(site_names), " site names (",
         paste(site_names, collapse = ", "), "); expected exactly one.")
  }
  final_step <- attrition %>% slice_max(step_order, n = 1, with_ties = FALSE)
  tibble(folder = basename(site_dir), site = site_names,
         cohort_n = as.numeric(final_step$n_remaining))
}

# Build the crosswalk for a set of per-site results folders. Returns a list with
#   $table   one row per cohort (folder, site, cohort_n, site_label), largest first
#   $aliases named character vector mapping BOTH the folder name and the in-file
#            site name to that cohort's "Site X" label
build_site_aliases <- function(site_dirs) {
  sizes <- map_dfr(site_dirs, read_site_cohort_size) %>%
    # Largest cohort first; ties broken by name so the labels are reproducible.
    arrange(desc(cohort_n), site) %>%
    mutate(site_label = paste("Site", site_alias_suffixes(n())))

  duplicated_keys <- sizes$site[duplicated(sizes$site)]
  if (length(duplicated_keys) > 0) {
    stop("The same site name appears in more than one results folder: ",
         paste(unique(duplicated_keys), collapse = ", "))
  }

  aliases <- c(setNames(sizes$site_label, sizes$site),
               setNames(sizes$site_label, sizes$folder))
  # A folder whose name equals its site name contributes the same pair twice.
  aliases <- aliases[!duplicated(names(aliases))]

  list(table = sizes, aliases = aliases)
}

# Map site names (or results-folder names) to their anonymized labels. Unknown
# values are an error rather than a silent pass-through, so a mislabelled table
# can never leak a real cohort name into a pooled figure.
anonymize_site <- function(x, aliases) {
  keys <- as.character(x)
  labels <- unname(aliases[keys])
  if (any(is.na(labels))) {
    stop("No anonymized label for site(s): ",
         paste(unique(keys[is.na(labels)]), collapse = ", "),
         ". Known sites: ", paste(unique(aliases), collapse = ", "))
  }
  labels
}

# Convenience wrapper for a data frame carrying a `site` column.
anonymize_site_column <- function(df, aliases) {
  if (!"site" %in% names(df)) {
    stop("Table has no `site` column to anonymize.")
  }
  mutate(df, site = anonymize_site(site, aliases))
}

# Console-only crosswalk so the coordinator can still read the pooled outputs.
print_site_alias_key <- function(site_aliases) {
  message("Anonymized cohort labels (largest cohort first; console only, ",
          "not written to the pooled outputs):")
  walk(seq_len(nrow(site_aliases$table)), function(i) {
    row <- site_aliases$table[i, ]
    message("  ", row$site_label, " = ", row$site,
            " (folder ", row$folder, "; N = ", format(row$cohort_n, big.mark = ","), ")")
  })
}
