# =============================================================================
# Script 05c: Age-correction sensitivity bracket (PBW / FVC_age25 / PFVC)
# =============================================================================
# The PBW-vs-PFVC difference is ENTIRELY the age channel (the two mechanics plots: PFVC
# differs from PBW by the discordance pbw/pfvc, which is ~deterministic in age). And the age
# channel is under-determined on THREE independent fronts: statistically (PFVC is collinear
# with age, so VT/PFVC's independent effect is unidentifiable -- the 05 PART 2f SE-explosion),
# physiologically (no normative strain-denominator data in the old/critically ill -- the
# FRC/EELV literature is middle-aged, where FRC ~ FVC and the proxies degenerate), and
# mechanically (E_rs-based plots share E_rs, so they can only SHOW the age re-weighting, not
# adjudicate it). No single analysis crowns a denominator.
#
# So this script does not try. It BRACKETS the age correction with three normalizers spanning
# the spectrum -- PBW (no age correction) / FVC_age25 (GLI height/sex/race structure but age
# pinned to 25, age-flat) / PFVC (full age correction) -- and reports 05's discordance and
# prognostic analyses across the bracket. The headline becomes a RANGE, not a point:
#   * stable across the bracket  => the clinical conclusion does not depend on the unknowable
#                                   age slope, and the whole age-identifiability worry is moot;
#   * fragile across the bracket => we have located the exact boundary of what the data support.
# FVC_age25 is the structural END of the bracket, agnostic about whether age-flat is correct --
# its job is to BOUND, not to be true. Derived in script 03 (pfvc_age25 + vtpfvc_age25 etc.).
#
# Decomposition used throughout (exact, by construction):
#   log(PBW/PFVC) = log(PBW/FVC_age25)  +  log(FVC_age25/PFVC)
#                   \___ structural ___/     \____ age ____/
#   log(VT/PFVC)  = log(VT/FVC_age25)   +  log(FVC_age25/PFVC)   (the age term = "disc_age")
# so "does PFVC beat FVC_age25" == "does the age correction earn its keep over structural size".
# Reads the same cross-sectional frame as 05; in-hospital mortality (deceased), as in 05.
# =============================================================================
library(tidyverse); library(arrow); library(here); library(splines)
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

auc_fn <- function(y, p) {                       # Mann-Whitney AUC (= C-statistic), as in 05
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc_optimism <- function(d, var, B) {            # apparent C minus bootstrap optimism (05 PART 2d)
  f <- reformulate(var, "deceased")
  app <- auc_fn(d$deceased, fitted(glm(f, data = d, family = binomial)))
  opt <- numeric(0)
  for (b in seq_len(B)) {
    bd <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
    m  <- tryCatch(glm(f, data = bd, family = binomial), error = function(e) NULL)
    if (is.null(m)) next
    opt <- c(opt, auc_fn(bd$deceased, fitted(m)) -
                  auc_fn(d$deceased, predict(m, newdata = d, type = "response")))
  }
  o <- mean(opt, na.rm = TRUE); if (!is.finite(o)) o <- 0
  tibble(c_apparent = app, c_corrected = app - o)
}

# VT cohort: full included cohort (NOT plateau-gated -- VT needs no dp). MP/Ers families
# inherit plateau-gating via their metric, exactly as in 05.
demog <- function(df) df %>%
  mutate(sex_category  = factor(sex_category,  levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         age10 = age_at_admission / 10, sf10 = sf_ratio / 10)
vt_base <- cross_sectional %>%
  filter(!is.na(pbw), pbw > 0, !is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0,
         !is.na(vtpbw), !is.na(vtpfvc), !is.na(vtpfvc_age25),
         !is.na(bmi), !is.na(sofa_total), !is.na(sf_ratio), !is.na(deceased)) %>%
  demog() %>%
  mutate(disc_total  = log(pbw / pfvc),          # = structural + age, by construction
         disc_struct = log(pbw / pfvc_age25),    # PBW vs structural size (GLI, age-flat)
         disc_age    = log(pfvc_age25 / pfvc),   # the age correction itself (FVC_age25 -> PFVC)
         l_vt25      = log(vtpfvc_age25))
mp_data  <- vt_base %>% filter(!is.na(mechanical_power), mechanical_power > 0,
                               !is.na(mp_pbw), !is.na(mp_pfvc), !is.na(mp_pfvc_age25))
ers_data <- vt_base %>% filter(!is.na(ers), ers > 0,
                               !is.na(ers_pbw), !is.na(ers_pfvc), !is.na(ers_pfvc_age25))
cat(sprintf("05c: VT cohort %d; MP %d; Ers %d; in-hospital mortality %.1f%%\n",
            nrow(vt_base), nrow(mp_data), nrow(ers_data), 100 * mean(vt_base$deceased)))

# =============================================================================
# A -- discordance decomposition: how much of PBW/PFVC is structural vs age?
# =============================================================================
# Decompose the cross-patient VARIATION (not the level: PBW/PFVC mixes kg and L, so its
# absolute level is a meaningless unit offset; the demographic bias lives in the spread).
# disc_total = disc_struct + disc_age exactly, so cov(total, age)/var(total) + cov(total,
# struct)/var(total) = 1 -- a clean variance partition of the size bias into its two legs.
r2 <- function(f) summary(lm(f, data = vt_base))$r.squared
vtot <- var(vt_base$disc_total)
decomp <- tibble(
  disc_struct_med = exp(median(vt_base$disc_struct)),   # reported for transparency (unit-laden)
  disc_age_med    = exp(median(vt_base$disc_age)),       # the age leg (FVC_age25/PFVC) -- unit-free
  var_age_share_pct    = 100 * cov(vt_base$disc_total, vt_base$disc_age)    / vtot,
  var_struct_share_pct = 100 * cov(vt_base$disc_total, vt_base$disc_struct) / vtot,
  age_R2_on_age        = r2(disc_age ~ ns(age10, 3)),                 # is the age leg just age? (~1)
  struct_R2_on_demo    = r2(disc_struct ~ sex_category + race_category + ns(age10, 3)))
write_csv(decomp, file.path(final_dir, paste0("norm_bracket_decomp_", site_name, ".csv")))

# =============================================================================
# B -- VT-tertile reclassification across the bracket (who moves, and on which leg)
# =============================================================================
tert <- vt_base %>% transmute(
  hospitalization_id, sex_category, race_category,
  age_grp = cut(age_at_admission, c(-Inf, 50, 65, 80, Inf), labels = c("<50","50-64","65-79","80+")),
  t_pbw  = ntile(vtpbw, 3), t_a25 = ntile(vtpfvc_age25, 3), t_pfvc = ntile(vtpfvc, 3))
reclass_leg <- function(a, b, leg) tert %>%
  summarise(leg = leg, pct = 100 * mean(.data[[a]] != .data[[b]])) %>%
  bind_cols(tert %>% group_by(age_grp) %>%
              summarise(p = 100 * mean(.data[[a]] != .data[[b]]), .groups = "drop") %>%
              pivot_wider(names_from = age_grp, values_from = p, names_prefix = "age_"))
reclass <- bind_rows(
  reclass_leg("t_pbw",  "t_a25",  "PBW -> FVC_age25 (structural leg)"),
  reclass_leg("t_a25",  "t_pfvc", "FVC_age25 -> PFVC (age leg)"),
  reclass_leg("t_pbw",  "t_pfvc", "PBW -> PFVC (total)"))
write_csv(reclass, file.path(final_dir, paste0("norm_bracket_reclassification_", site_name, ".csv")))

# =============================================================================
# C -- prognostic discrimination across the bracket (optimism-corrected univariate C)
# =============================================================================
B_DISC <- if (identical(site_name, "synthetic_clif")) 100L else 300L
set.seed(20260620)
sweep <- tribble(
  ~family, ~norm,        ~var,             ~data,
  "VT",    "PBW",        "vtpbw",          "vt",
  "VT",    "FVC_age25",  "vtpfvc_age25",   "vt",
  "VT",    "PFVC",       "vtpfvc",         "vt",
  "MP",    "PBW",        "mp_pbw",         "mp",
  "MP",    "FVC_age25",  "mp_pfvc_age25",  "mp",
  "MP",    "PFVC",       "mp_pfvc",        "mp",
  "Ers",   "PBW",        "ers_pbw",        "ers",
  "Ers",   "FVC_age25",  "ers_pfvc_age25", "ers",
  "Ers",   "PFVC",       "ers_pfvc",       "ers")
dsets <- list(vt = vt_base, mp = mp_data, ers = ers_data)
disc <- pmap_dfr(sweep, function(family, norm, var, data) {
  d <- dsets[[data]] %>% filter(is.finite(.data[[var]]))
  bind_cols(tibble(family = family, norm = factor(norm, c("PBW","FVC_age25","PFVC")), metric = var, n = nrow(d)),
            auc_optimism(d, var, B_DISC))
}) %>% arrange(family, norm)
write_csv(disc, file.path(final_dir, paste0("norm_bracket_discrimination_", site_name, ".csv")))

# =============================================================================
# D -- does the AGE correction add prognostic info OVER structural size? (the adjudicator)
# =============================================================================
# log(VT/PFVC) = log(VT/FVC_age25) + disc_age, so adding disc_age to a model with the
# structural-normalized VT recovers the full PFVC normalization. LRT of that addition =
# "does PFVC's age slope carry mortality signal beyond structural size + the covariates?"
# WITHOUT age adjustment it can; WITH age10 in the model it cannot if it is just age
# (the unidentifiability) -- the adjusted-vs-unadjusted gap localizes the channel.
lrt_age <- function(adjust) {
  cov <- if (adjust) " + vtpbw + sofa_total + sf10 + bmi + age10 + sex_category + race_category" else ""
  m0 <- glm(as.formula(paste0("deceased ~ l_vt25", cov)), vt_base, family = binomial)
  m1 <- glm(as.formula(paste0("deceased ~ l_vt25 + disc_age", cov)), vt_base, family = binomial)
  an <- anova(m0, m1, test = "LRT")
  tibble(adjustment = if (adjust) "adjusted (+age/sex/race)" else "unadjusted",
         lrt_p = an$`Pr(>Chi)`[2],
         dC = auc_fn(vt_base$deceased, fitted(m1)) - auc_fn(vt_base$deceased, fitted(m0)),
         or_disc_age = exp(coef(m1)[["disc_age"]]))
}
age_add <- bind_rows(lrt_age(FALSE), lrt_age(TRUE))
write_csv(age_add, file.path(final_dir, paste0("norm_bracket_age_increment_", site_name, ".csv")))

# ---- console summary --------------------------------------------------------
cat("\n=== 05c age-correction bracket (PBW / FVC_age25 / PFVC) ===\n")
cat(sprintf("A. Variation in the PBW/PFVC size bias: %.0f%% structural (height/sex/race), %.0f%% age slope\n",
            decomp$var_struct_share_pct, decomp$var_age_share_pct))
cat(sprintf("   age leg median FVC_age25/PFVC = %.2f; age leg R2-on-age = %.2f (>~0.9 => the age leg IS age)\n",
            decomp$disc_age_med, decomp$age_R2_on_age))
cat("\nB. VT-tertile reclassification (% moving), by leg:\n")
print(as.data.frame(reclass %>% mutate(across(where(is.numeric), ~round(.x, 1)))), row.names = FALSE)
cat("\nC. Optimism-corrected discrimination across the bracket (the RANGE is the headline):\n")
disc %>% group_by(family) %>% group_walk(function(g, k) {
  br <- sprintf("%s=%.3f", g$norm, g$c_corrected)
  cat(sprintf("   [%-3s] %s   | structural gain (FVC_age25-PBW)=%+.3f, age gain (PFVC-FVC_age25)=%+.3f\n",
              k$family, paste(br, collapse = "  "),
              g$c_corrected[g$norm == "FVC_age25"] - g$c_corrected[g$norm == "PBW"],
              g$c_corrected[g$norm == "PFVC"] - g$c_corrected[g$norm == "FVC_age25"]))
})
cat("\nD. Does the age correction add over structural size? (disc_age LRT)\n")
print(as.data.frame(age_add %>% mutate(lrt_p = signif(lrt_p, 2),
        dC = round(dC, 4), or_disc_age = round(or_disc_age, 2))), row.names = FALSE)
cat("   adjusted p NS / dC~0 while unadjusted significant => the age slope is age, not size (unidentifiable).\n")
message("Wrote norm_bracket_{decomp,reclassification,discrimination,age_increment}_", site_name, ".csv to ", final_dir)
