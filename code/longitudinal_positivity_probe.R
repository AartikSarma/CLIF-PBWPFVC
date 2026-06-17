# =============================================================================
# PROBE: decision-level positivity for a longitudinal stress-limiting TTE
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# NOT part of the pipeline. A cheap go/no-go check BEFORE building the full
# longitudinal target trial emulation (MSM / clone-censor-weight). The TTE
# escapes the cross-sectional positivity wall ONLY if, at the day level, the
# clinical decision to de-escalate tidal volume is NOT deterministic in measured
# state. This builds the VT/PFVC trajectory from the cleaned waterfall, defines
# the de-escalation decision among above-threshold patient-days, fits the
# decision propensity, and reports overlap (the analog of the cross-sectional
# c=0.996 / ESS~0 diagnostic, now per decision).
#
# READS ONLY: resp_support_waterfall_clean.parquet (longitudinal VT settings) +
# analysis_cross_sectional.parquet (baseline PFVC + demographics + severity).
# Needs NO mortality, so the synthetic-mortality bug is irrelevant here.
# =============================================================================

library(tidyverse)
library(arrow)
library(here)

source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

VTPFVC_TARGET <- 11      # stress threshold: VT/PFVC < 11% of predicted FVC
DEESC_FRAC    <- 0.05    # de-escalation = next-day VT reduced by >= 5% (relative)
MAX_DAY       <- 13      # first 14 vent-days
message("Site: ", site_name)

# =============================================================================
# 1. Baseline PFVC + demographics + t0, and the longitudinal VT settings
# =============================================================================
base <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet")) %>%
  transmute(hospitalization_id, t0 = recorded_dttm, pfvc, age_at_admission,
            sex_category, race_category, sofa_total)

wf <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  select(hospitalization_id, recorded_dttm, tidal_volume_set, fio2_set, peep_set,
         resp_rate_set, plateau_pressure_obs) %>%
  filter(!is.na(tidal_volume_set), tidal_volume_set > 0) %>%
  inner_join(base, by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_DAY,
         !is.na(pfvc), pfvc > 0,
         !is.na(age_at_admission), !is.na(sex_category), !is.na(race_category)) %>%
  mutate(vtpfvc = tidal_volume_set / pfvc * 0.1)          # same units as script 03

message("Ventilated patient-timepoints (days 0-", MAX_DAY, "): ", nrow(wf),
        " across ", n_distinct(wf$hospitalization_id), " hospitalizations")

# =============================================================================
# 2. Collapse to one row per (hospitalization, vent-day); define the decision
# =============================================================================
# plateau_pressure_obs is NOT forward-filled (project rule): daily value = median
# of RECORDED values that day, NA if none.
daily <- wf %>%
  group_by(hospitalization_id, vent_day, age_at_admission, sex_category,
           race_category, sofa_total) %>%
  summarise(vtpfvc  = median(vtpfvc, na.rm = TRUE),
            fio2    = median(fio2_set, na.rm = TRUE),
            peep    = median(peep_set, na.rm = TRUE),
            rr      = median(resp_rate_set, na.rm = TRUE),
            plateau = median(plateau_pressure_obs, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(hospitalization_id, vent_day) %>%
  group_by(hospitalization_id) %>%
  mutate(vtpfvc_next = lead(vtpfvc),
         has_next    = !is.na(vtpfvc_next)) %>%
  ungroup()

# --- full time-varying confounders: S/F + MAP (vitals), vasopressor (meds) -----
# (lactate is NOT in the pulled labs, so it cannot be added here -- noted as a
# limitation; the full TTE would add it where available.)
vit <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category %in% c("spo2", "map")) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_DAY) %>%
  group_by(hospitalization_id, vent_day, vital_category) %>%
  summarise(v = median(vital_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = vital_category, values_from = v)        # -> spo2, map

med <- read_parquet(file.path(output_dir, "cohort_meds.parquet")) %>%
  filter(med_group == "vasoactives") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(admin_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_DAY) %>%
  distinct(hospitalization_id, vent_day) %>% mutate(on_pressor = 1L)

daily <- daily %>%
  left_join(vit, by = c("hospitalization_id", "vent_day")) %>%
  left_join(med, by = c("hospitalization_id", "vent_day")) %>%
  mutate(on_pressor = coalesce(on_pressor, 0L),
         fio2_frac  = if_else(fio2 > 1.5, fio2 / 100, fio2),       # robust to %/fraction
         sf         = spo2 / fio2_frac)                            # S/F ratio

# DECISION POINTS = above-threshold days with a next day observed. Outcome =
# did the clinician de-escalate VT (>= DEESC_FRAC relative reduction) next day?
# Complete-case on the FULL confounder set so the spec comparison is apples-to-apples.
decision <- daily %>%
  filter(has_next, vtpfvc >= VTPFVC_TARGET, is.finite(fio2), is.finite(peep),
         is.finite(rr), is.finite(vtpfvc), is.finite(sf), is.finite(map)) %>%
  mutate(de_escalate = as.integer(vtpfvc_next < vtpfvc * (1 - DEESC_FRAC)))
message("Decision points (above-threshold patient-days w/ next day): ", nrow(decision),
        " | de-escalation rate = ", round(100 * mean(decision$de_escalate), 1), "%")

# =============================================================================
# 3. Decision propensity + overlap (the positivity verdict)
# =============================================================================
# If demographics/clinical state make the decision near-deterministic (mass at
# 0/1, AUC ~1, ESS ~0) -> positivity dead even longitudinally. If overlap exists
# -> the TTE is alive and worth building.
auc_fn <- function(score, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
overlap_for <- function(df, cov_str, spec) {
  fit <- glm(as.formula(paste("de_escalate ~", cov_str)), data = df, family = binomial)
  ps  <- as.numeric(predict(fit, type = "response"))
  p_m <- mean(df$de_escalate)
  sw  <- ifelse(df$de_escalate == 1, p_m / ps, (1 - p_m) / (1 - ps))
  tibble(spec = spec, n = nrow(df), deesc_rate = p_m,
         ps_min = min(ps), ps_max = max(ps),
         frac_extreme = mean(ps < 0.05 | ps > 0.95),
         ess_frac = (sum(sw)^2 / sum(sw^2)) / nrow(df),
         auc = auc_fn(ps, df$de_escalate))
}
# minimal resp-state spec -> + full time-varying confounders -> + demographics.
# The key test: does the STRICTER conditioning (full set) push AUC toward
# determinism / collapse overlap, or does positivity survive?
clin_min  <- "vtpfvc + fio2 + peep + rr + vent_day + sofa_total"
clin_full <- paste(clin_min, "+ sf + map + on_pressor")
clin      <- paste(clin_full, "+ age_at_admission + sex_category + race_category")
ov <- bind_rows(
  overlap_for(decision, clin_min,  "resp_state_min"),
  overlap_for(decision, clin_full, "full_confounders"),
  overlap_for(decision, clin,      "full_plus_demographics"))
write_csv(ov, file.path(final_dir, paste0("probe_decision_overlap_", site_name, ".csv")))
cat("\n=== DECISION-LEVEL OVERLAP (positivity probe) ===\n")
print(as.data.frame(ov %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)

# verdict heuristic
v <- ov %>% filter(spec == "full_plus_demographics")
verdict <- if (v$frac_extreme < 0.15 && v$ess_frac > 0.4 && v$auc < 0.9)
  "ALIVE: decision-level overlap present -> longitudinal TTE worth building" else if
  (v$frac_extreme > 0.4 || v$ess_frac < 0.1 || v$auc > 0.97)
  "DEAD: decision near-deterministic in state -> positivity fails longitudinally too" else
  "MARGINAL: partial overlap -> trim extreme strata / refine decision def before committing"
message("\nVERDICT (+demographics spec): ", verdict)

# propensity overlap figure (de-escalate vs not), +demographics spec
fit_d <- glm(as.formula(paste("de_escalate ~", clin)), data = decision, family = binomial)
plotdat <- decision %>% mutate(ps = as.numeric(predict(fit_d, type = "response")),
                               arm = ifelse(de_escalate == 1, "de-escalated", "not de-escalated"))
p <- ggplot(plotdat, aes(ps, fill = arm, colour = arm)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(0.05, 0.95), linetype = 3, colour = "grey50") +
  scale_fill_manual(values = c("de-escalated" = okabe[3], "not de-escalated" = okabe[1]),
                    aesthetics = c("fill", "colour"), name = NULL) +
  labs(x = "Propensity to de-escalate VT (next day) | state + demographics", y = "Density",
       title = "Decision-level positivity probe: overlap on the de-escalation decision",
       subtitle = paste0(site_name, " - overlap here (unlike the cross-sectional c=0.996) ",
                         "means the longitudinal TTE escapes the positivity wall")) +
  theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("probe_decision_overlap_", site_name, ".pdf")),
       p, width = 9, height = 4.5)
message("Wrote probe overlap table + figure to ", final_dir)
