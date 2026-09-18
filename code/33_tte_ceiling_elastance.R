# =============================================================================
# 33_tte_ceiling_elastance: effect modification by MEASURED mis-sizing (specific elastance)
# =============================================================================
# Every modifier used so far -- PBW/PFVC discordance, PFVC itself -- is a deterministic
# function of age, sex, race and height (tte_within_demographic_hte: R^2 ~ 0.99), so "who benefits" has been a
# demographic statement and the emulation cannot say whether the target is lung size or age.
# This leaf runs the SAME clone-censor-weight design and the SAME estimator on a modifier that
# is MEASURED rather than predicted:
#
#   specific elastance  E_spec = Ers x PFVC   (cmH2O)
#
# Ers = driving pressure / tidal volume at the index timepoint. If PFVC is right about this
# patient's aerated volume, Chiumello's result puts E_spec near a physiologic constant; ABOVE it
# the aerated lung is SMALLER than PFVC predicts (the baby-lung direction), BELOW it larger. So
# E_spec is a per-patient, measured index of mis-sizing, pointing the same way as high PBW/PFVC
# but carrying information PBW/PFVC cannot: Ers is not a function of demographics, so E_spec
# VARIES WITHIN demographic strata. Read A quantifies exactly that (the R^2 of each modifier on
# age/sex/race/height, the analogue of the 4j identifying-variation table) -- it is the reason
# the leaf exists, and it should be reported whatever the CATE shows.
#
# Reads, all on the primary (bite-matched head-to-head) design unless PBWPFVC_CEIL_DIAG_DESIGN
# says otherwise:
#   A. identifying variation: R^2 of log E_spec vs log(PBW/PFVC) on demographics, and the
#      cross-tab of their tertiles (are these the same patients?)
#   B. tertile RDs + the Stiff-minus-Compliant gradient (SOFA- and BMI-adjusted MSM)
#   C. continuous CATE over log E_spec, absolute (risk difference) AND relative (risk ratio):
#      the modifier marks sicker patients, so a sloped RD with a flat RR is baseline risk
#   D. positivity by E_spec tertile x arm
#
# CAVEATS, all of which belong in the manuscript. E_spec needs a recorded plateau, so this runs
# on a subset (reported). Ers is RESPIRATORY-SYSTEM elastance: a stiff chest wall (obesity,
# abdominal pressure) raises it without the lung being small, so BMI is in the adjustment set and
# a BMI-stratified read is reported. Only the BASELINE value is used -- later elastance is on the
# causal path from the policy. And Chiumello's constancy is a GROUP MEAN with a wide individual
# spread, so E_spec is noisy per patient: read the smooth curve and the tertiles, never a
# threshold. Env: PBWPFVC_CEIL_DIAG_DESIGN (ceiling|cap), PBWPFVC_TTE_EXPO_FAMILY, PBWPFVC_NBOOT.
# =============================================================================
library(here)
source(here::here("code", "33_tte_ceiling_common.R"))

DESIGN <- Sys.getenv("PBWPFVC_CEIL_DIAG_DESIGN", "ceiling")
stopifnot(DESIGN %in% c("ceiling", "cap"))
prefix <- design_prefix(DESIGN)
tag    <- function(x) paste0(prefix, "ers_", x, "_", site_name)
ERS_LEVELS <- c("Compliant", "Mid", "Stiff")   # low / mid / high specific elastance

# --- the modifier ------------------------------------------------------------------------
# E_spec = Ers x PFVC, at the index timepoint (baseline; never a post-policy value).
emod <- base %>%
  transmute(hospitalization_id, ers, bmi, pfvc, pbw, age10, sex_category, race_category,
            sofa_total, height_cm, spec_ers = ers * pfvc, disc = pbw / pfvc) %>%
  # BMI is in the adjustment set (chest wall), so require it: a missing covariate would
  # propagate NA through the standardization grid and silently null out every estimate.
  filter(is.finite(spec_ers), spec_ers > 0, is.finite(disc), is.finite(bmi), is.finite(sofa_total))
n_all <- nrow(base)
message(sprintf("33_tte_ceiling_elastance [%s/%s]: %d of %d patients have a baseline plateau (E_spec); median E_spec %.1f cmH2O",
                FAM, DESIGN, nrow(emod), n_all, median(emod$spec_ers)))
if (nrow(emod) < 300) stop("33_tte_ceiling_elastance: fewer than 300 patients with a baseline E_spec -- not estimable here.")
emod <- emod %>%
  mutate(ers_grp = cut(spec_ers, c(-Inf, quantile(spec_ers, c(1/3, 2/3)), Inf), labels = ERS_LEVELS),
         disc_grp3 = cut(disc, c(-Inf, quantile(disc, c(1/3, 2/3)), Inf), labels = DISC_LEVELS),
         bmi_grp = cut(bmi, c(-Inf, 30, Inf), labels = c("BMI < 30", "BMI >= 30")),
         l_ers = log(spec_ers), l_disc = log(disc))

# =============================================================================
# A. Identifying variation: is the measured modifier demographic, like the predicted one?
# =============================================================================
r2_on_demo <- function(v, extra = NULL) {
  rhs <- paste(c("ns(age10, 4)", "sex_category", "race_category", "height_cm", extra), collapse = " + ")
  summary(lm(as.formula(paste(v, "~", rhs)), data = emod))$r.squared
}
idvar <- tibble(
  modifier = c("log specific elastance (Ers x PFVC, measured)", "log PBW/PFVC discordance (predicted)"),
  r2_on_demographics = c(r2_on_demo("l_ers"), r2_on_demo("l_disc")),
  r2_on_demographics_plus_bmi = c(r2_on_demo("l_ers", "bmi"), r2_on_demo("l_disc", "bmi")),
  sd_log = c(sd(emod$l_ers), sd(emod$l_disc)), n = nrow(emod)) %>%
  mutate(residual_variance_share = 1 - r2_on_demographics, family = FAM, design = DESIGN, site = site_name)
crosstab <- emod %>% count(ers_grp, disc_grp3) %>%
  group_by(ers_grp) %>% mutate(frac_of_ers_tertile = n / sum(n)) %>% ungroup() %>%
  mutate(family = FAM, design = DESIGN, site = site_name)
write_csv(idvar,    file.path(final_dir, paste0(tag("identifying_variation"), ".csv")))
write_csv(crosstab, file.path(final_dir, paste0(tag("crosstab"), ".csv")))

# =============================================================================
# B-D. The same design and estimator as the discordance reads, keyed on E_spec
# =============================================================================
des  <- build_ceiling_design(DESIGN, keep_pday = TRUE)
long <- des$long %>% inner_join(emod %>% select(hospitalization_id, ers_grp, l_ers, sofa_total, bmi, bmi_grp),
                                by = "hospitalization_id")
ids  <- unique(long$hospitalization_id)
prof <- long %>% distinct(hospitalization_id, ers_grp, l_ers, sofa_total, bmi)
message(sprintf("  design '%s': %d patients carry both the design and a baseline E_spec", DESIGN, length(ids)))

FORM_E <- died ~ arm * ers_grp + arm * ns(day, 4) + ers_grp * ns(day, 4) + sofa_total + bmi
std_rd_e <- function(fit, pf, g) {
  pp <- if (g == "All") pf else pf %>% filter(as.character(ers_grp) == g)
  if (nrow(pp) < 50) return(NA_real_)
  cells <- pp %>% count(ers_grp, sofa_total, bmi, name = "wt")
  grid  <- tidyr::crossing(cells, day = 1:HORIZON, arm = arm_f())
  grid$haz <- predict(fit, grid, type = "response")
  ci <- grid %>% group_by(arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
    group_by(arm) %>% arrange(day) %>% summarise(cif = 1 - prod(1 - h), .groups = "drop")
  ci$cif[ci$arm == "strain_limiting"] - ci$cif[ci$arm == "permissive"]
}
# continuous: fixed spline basis on log E_spec, standardized over the marginal (SOFA, BMI) cells
KD   <- as.integer(Sys.getenv("PBWPFVC_DISC_CATE_DF", "3"))
Bspl <- ns(prof$l_ers, df = KD); zc <- paste0("z", seq_len(KD))
add_z <- function(d) { m <- predict(Bspl, d$l_ers); for (j in seq_len(KD)) d[[zc[j]]] <- m[, j]; d }
long_z <- add_z(long)
FORM_C <- as.formula(paste0("died ~ arm*(", paste(zc, collapse = "+"), ") + arm*ns(day,4) + (",
                            paste(zc, collapse = "+"), ")*ns(day,4) + sofa_total + bmi"))
GRID  <- seq(quantile(prof$l_ers, .025), quantile(prof$l_ers, .975), length.out = 40)
qe    <- quantile(prof$l_ers, c(.10, .90))
EVAL  <- sort(unique(c(GRID, qe)))
cells <- prof %>% count(sofa_total, bmi, name = "wt")
grid_c <- tidyr::crossing(l_ers = EVAL, cells, day = 1:HORIZON, arm = arm_f()) %>% add_z()
cate_from <- function(fit) {
  g <- grid_c; g$haz <- predict(fit, g, type = "response")
  g %>% group_by(l_ers, arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
    group_by(l_ers, arm) %>% arrange(day) %>% summarise(cif = 1 - prod(1 - h), .groups = "drop") %>%
    group_by(l_ers) %>%
    summarise(rd = cif[arm == "strain_limiting"] - cif[arm == "permissive"],
              rr = cif[arm == "strain_limiting"] / cif[arm == "permissive"], .groups = "drop") %>%
    arrange(l_ers)
}
fit_t <- suppressWarnings(glm(FORM_E, data = long,   family = binomial, weights = ipcw))
fit_c <- suppressWarnings(glm(FORM_C, data = long_z, family = binomial, weights = ipcw))
pt_t  <- vapply(ERS_LEVELS, function(g) std_rd_e(fit_t, prof, g), numeric(1))
grad_pt <- unname(pt_t["Stiff"] - pt_t["Compliant"])
c0 <- cate_from(fit_c)
i10 <- which.min(abs(c0$l_ers - qe[1])); i90 <- which.min(abs(c0$l_ers - qe[2]))
slope_pt <- c0$rd[i90] - c0$rd[i10]; rr_ratio_pt <- c0$rr[i90] / c0$rr[i10]

# --- one cluster bootstrap over everything ------------------------------------------------
NC <- length(EVAL)
tmpl <- c(setNames(rep(NA_real_, 3), ERS_LEVELS), gradient = NA_real_, slope = NA_real_,
          rr_ratio = NA_real_, setNames(rep(NA_real_, NC), paste0("c", seq_len(NC))),
          setNames(rep(NA_real_, NC), paste0("r", seq_len(NC))))
boot_one <- function() {
  samp <- tibble(hospitalization_id = sample(ids, replace = TRUE))
  out <- tmpl
  lb <- long_z %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  pb <- samp %>% left_join(prof, by = "hospitalization_id")
  ft <- suppressWarnings(tryCatch(glm(FORM_E, data = lb, family = binomial, weights = ipcw), error = function(e) NULL))
  if (!is.null(ft)) {
    r <- vapply(ERS_LEVELS, function(g) std_rd_e(ft, pb, g), numeric(1))
    out[ERS_LEVELS] <- r; out["gradient"] <- unname(r["Stiff"] - r["Compliant"])
  }
  fc <- suppressWarnings(tryCatch(glm(FORM_C, data = lb, family = binomial, weights = ipcw), error = function(e) NULL))
  if (!is.null(fc)) {
    cv <- tryCatch(cate_from(fc), error = function(e) NULL)
    if (!is.null(cv) && nrow(cv) == NC) {
      out[paste0("c", seq_len(NC))] <- cv$rd; out[paste0("r", seq_len(NC))] <- cv$rr
      out["slope"] <- cv$rd[i90] - cv$rd[i10]; out["rr_ratio"] <- cv$rr[i90] / cv$rr[i10]
    }
  }
  out
}
nc_used <- min(N_CORES, N_BOOT)
message("33_tte_ceiling_elastance cluster bootstrap (", N_BOOT, " reps across ", nc_used, " core(s)) ...")
chunks <- split(seq_len(N_BOOT), cut(seq_len(N_BOOT), min(20L, N_BOOT), labels = FALSE))
bl <- vector("list", N_BOOT); done <- 0L
if (nc_used > 1) {
  cl <- makeCluster(nc_used, type = "PSOCK")
  clusterEvalQ(cl, { library(tidyverse); library(splines) })
  clusterExport(cl, envir = .GlobalEnv, varlist = c("HORIZON", "arm_f"))
  clusterExport(cl, envir = environment(), varlist = c(
    "ids", "long_z", "prof", "tmpl", "boot_one", "FORM_E", "FORM_C", "cate_from", "std_rd_e",
    "grid_c", "NC", "ERS_LEVELS", "i10", "i90"))
  clusterSetRNGStream(cl, 20260911L)
  tryCatch(for (ch in chunks) {
    bl[ch] <- parLapply(cl, ch, function(b) boot_one()); done <- done + length(ch)
    message(sprintf("  bootstrap %d/%d", done, N_BOOT))
  }, finally = stopCluster(cl))
} else for (ch in chunks) { for (b in ch) bl[[b]] <- boot_one(); done <- done + length(ch)
  message(sprintf("  bootstrap %d/%d", done, N_BOOT)) }
bts <- do.call(rbind, bl)
ci <- function(col) quantile(bts[, col], c(.025, .975), na.rm = TRUE)

# --- tables ---------------------------------------------------------------------------------
n_e <- prof %>% count(ers_grp) %>% mutate(ers_grp = as.character(ers_grp))
hte <- tibble(ers_grp = ERS_LEVELS, rd = pt_t[ERS_LEVELS],
              rd_lo = vapply(ERS_LEVELS, function(k) ci(k)[1], numeric(1)),
              rd_hi = vapply(ERS_LEVELS, function(k) ci(k)[2], numeric(1)),
              adjustment = "sofa + bmi", family = FAM, design = DESIGN, site = site_name) %>%
  left_join(n_e, by = "ers_grp") %>% mutate(ers_grp = factor(ers_grp, ERS_LEVELS))
gradient <- tibble(family = FAM, design = DESIGN,
                   statistic = c("RD(Stiff) - RD(Compliant)", "RD(p90 E_spec) - RD(p10 E_spec)",
                                 "RR(p90 E_spec) / RR(p10 E_spec)"),
                   estimate = c(grad_pt, slope_pt, rr_ratio_pt),
                   lo = c(ci("gradient")[1], ci("slope")[1], ci("rr_ratio")[1]),
                   hi = c(ci("gradient")[2], ci("slope")[2], ci("rr_ratio")[2]),
                   scale = c("risk difference", "risk difference", "ratio of risk ratios"),
                   site = site_name)
cmat <- bts[, paste0("c", seq_len(NC)), drop = FALSE]; rmat <- bts[, paste0("r", seq_len(NC)), drop = FALSE]
curve <- tibble(family = FAM, design = DESIGN, l_ers = c0$l_ers, spec_ers = exp(c0$l_ers),
                rd = c0$rd, rd_lo = apply(cmat, 2, quantile, .025, na.rm = TRUE),
                rd_hi = apply(cmat, 2, quantile, .975, na.rm = TRUE),
                rr = c0$rr, rr_lo = apply(rmat, 2, quantile, .025, na.rm = TRUE),
                rr_hi = apply(rmat, 2, quantile, .975, na.rm = TRUE), site = site_name) %>%
  filter(l_ers %in% GRID)
write_csv(hte,      file.path(final_dir, paste0(tag("hte"), ".csv")))
write_csv(gradient, file.path(final_dir, paste0(tag("gradient"), ".csv")))
write_csv(curve,    file.path(final_dir, paste0(tag("cate_curve"), ".csv")))

# positivity by E_spec tertile x arm, and the BMI-stratified gradient (chest-wall caveat)
pday <- bind_rows(des$bl$pday %>% mutate(arm = "pfvc_informed"), des$bh$pday %>% mutate(arm = "pbw_ceiling")) %>%
  inner_join(emod %>% select(hospitalization_id, ers_grp), by = "hospitalization_id")
overlap <- pday %>% group_by(arm, ers_grp) %>%
  summarise(n_eligible_days = n(), median_padhere = median(p_adhere),
            frac_padhere_lt05 = mean(p_adhere < 0.05), .groups = "drop") %>%
  mutate(family = FAM, design = DESIGN, site = site_name)
bmi_strat <- map_dfr(levels(emod$bmi_grp), function(b) {
  ids_b <- emod$hospitalization_id[emod$bmi_grp == b & !is.na(emod$bmi_grp)]
  lb <- long %>% filter(hospitalization_id %in% ids_b); pb <- prof %>% filter(hospitalization_id %in% ids_b)
  if (n_distinct(lb$hospitalization_id) < 200) return(tibble())
  fb <- suppressWarnings(tryCatch(glm(FORM_E, data = lb, family = binomial, weights = ipcw), error = function(e) NULL))
  if (is.null(fb)) return(tibble())
  r <- vapply(ERS_LEVELS, function(g) std_rd_e(fb, pb, g), numeric(1))
  tibble(bmi_grp = b, n = n_distinct(lb$hospitalization_id), rd_compliant = r["Compliant"],
         rd_stiff = r["Stiff"], gradient = unname(r["Stiff"] - r["Compliant"]))
}) %>% mutate(family = FAM, design = DESIGN, site = site_name)
write_csv(overlap,   file.path(final_dir, paste0(tag("overlap"), ".csv")))
if (nrow(bmi_strat)) write_csv(bmi_strat, file.path(final_dir, paste0(tag("bmi_strata"), ".csv")))

# --- figure: absolute over relative, with a rug of the modifier ------------------------------
rug <- tibble(x = exp(prof$l_ers))
p_rd <- ggplot(curve, aes(spec_ers, 100 * rd)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_ribbon(aes(ymin = 100 * rd_lo, ymax = 100 * rd_hi), alpha = .15, fill = "#0072B2") +
  geom_line(colour = "#0072B2", linewidth = 1) +
  geom_rug(data = rug %>% slice_sample(n = min(nrow(rug), 2000)), aes(x = x), inherit.aes = FALSE,
           alpha = .06, length = unit(.03, "npc")) +
  labs(x = NULL, y = "CATE: 28-d mortality RD (pp)",
       title = "A. Absolute scale", subtitle = sprintf("RD(p90)-RD(p10): %+.2f pp [%.2f, %.2f]",
                100 * slope_pt, 100 * ci("slope")[1], 100 * ci("slope")[2])) +
  theme_minimal(base_size = 11)
p_rr <- ggplot(curve, aes(spec_ers, rr)) +
  geom_hline(yintercept = 1, colour = "grey80") +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = .15, fill = "#D55E00") +
  geom_line(colour = "#D55E00", linewidth = 1) +
  geom_rug(data = rug %>% slice_sample(n = min(nrow(rug), 2000)), aes(x = x), inherit.aes = FALSE,
           alpha = .06, length = unit(.03, "npc")) +
  scale_y_log10() +
  labs(x = "Specific elastance  Ers x PFVC (cmH2O) - higher = aerated lung SMALLER than predicted",
       y = "Relative CATE: risk ratio", title = "B. Relative scale",
       subtitle = sprintf("RR(p90)/RR(p10): %.2f [%.2f, %.2f]", rr_ratio_pt, ci("rr_ratio")[1], ci("rr_ratio")[2])) +
  theme_minimal(base_size = 11)
fig <- p_rd / p_rr + patchwork::plot_annotation(
  title = paste0("PFVC-guided dosing by MEASURED mis-sizing (specific elastance) - ", site_name,
                 if (is_synthetic) " (SYNTHETIC - plumbing only)" else ""),
  subtitle = sprintf("Design '%s' (%s family). Unlike PBW/PFVC, E_spec is measured: R2 on demographics %.2f vs %.2f. Negative = the PFVC-informed arm lowers mortality more.",
                     DESIGN, FAM, idvar$r2_on_demographics[1], idvar$r2_on_demographics[2]))
ggsave(file.path(final_dir, paste0(tag("cate"), ".pdf")), fig, width = 9, height = 8)

# --- console ---------------------------------------------------------------------------------
cat(sprintf("\n=== 33_tte_ceiling_elastance [%s/%s]: %d patients ===\n", FAM, DESIGN, length(ids)))
cat("--- A. identifying variation (the point of the leaf) ---\n")
print(as.data.frame(idvar %>% transmute(modifier, r2_demo = round(r2_on_demographics, 3),
        r2_demo_bmi = round(r2_on_demographics_plus_bmi, 3), resid_share = round(residual_variance_share, 3))), row.names = FALSE)
cat("--- cross-tab: E_spec tertile x discordance tertile (row %) ---\n")
print(as.data.frame(crosstab %>% transmute(ers_grp, disc_grp3, n, pct = round(100 * frac_of_ers_tertile))), row.names = FALSE)
cat("--- B. RD by E_spec tertile (SOFA + BMI adjusted) ---\n")
print(as.data.frame(hte %>% transmute(ers_grp, n, rd_pp = round(100 * rd, 2),
        ci = sprintf("[%.2f, %.2f]", 100 * rd_lo, 100 * rd_hi))), row.names = FALSE)
print(as.data.frame(gradient %>% transmute(statistic, est = round(ifelse(scale == "risk difference", 100 * estimate, estimate), 3),
        ci = sprintf("[%.3f, %.3f]", ifelse(scale == "risk difference", 100 * lo, lo),
                     ifelse(scale == "risk difference", 100 * hi, hi)))), row.names = FALSE)
cat("--- D. positivity by E_spec tertile ---\n")
print(as.data.frame(overlap %>% transmute(arm, ers_grp, n_eligible_days, median_padhere = round(median_padhere, 3),
        frac_lt05 = round(frac_padhere_lt05, 3))), row.names = FALSE)
if (nrow(bmi_strat)) { cat("--- chest-wall check: gradient within BMI strata ---\n")
  print(as.data.frame(bmi_strat %>% transmute(bmi_grp, n, gradient_pp = round(100 * gradient, 2))), row.names = FALSE) }
message("Wrote ", tag("{identifying_variation,crosstab,hte,gradient,cate_curve,overlap,bmi_strata}"), ".csv + cate .pdf to ", final_dir)
