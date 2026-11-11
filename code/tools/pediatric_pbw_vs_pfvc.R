#!/usr/bin/env Rscript
# =============================================================================
# pediatric_pbw_vs_pfvc.R -- presentation figures: pediatric PBW formulas against GLI-2012
#                            predicted FVC, ages 3-18, by sex, height percentile and race
# =============================================================================
# Formulas only; no patient data. For boys and girls aged 3-18 at the 10th, 50th and 90th
# percentile of height for age (CDC 2000 stature-for-age), computes the predicted body weight
# each pediatric formula assigns and the race-specific GLI-2012 predicted FVC for each of the
# five GLI ethnicity groups, then the strain a 6 mL/kg PBW tidal volume delivers:
#
#   VT/PFVC at 6 mL/kg (%) = 6 mL/kg x PBW / PFVC
#
# For comparison, one adult heatmap applies the same measure to Devine PBW at the midpoint of each
# decade from the 20s to the 80s. Adult heights are held at the CDC age-20 10th, 50th and 90th
# percentiles at every age, so the decades differ only through GLI's age term; the adult and
# pediatric heatmaps share one colour scale.
#
# the same quantity as figure 1C (pooling/pooled_displays.R) for adults. The dashed line is
# that figure's reference, 11% of predicted FVC (about the 75th percentile of VT/PFVC in the
# ARMA low tidal volume arm).
#
# PBW formulas (H = height in cm):
#   Traub-Johnson (1980)  2.396 x exp(0.01863 H); derived in children 1-17 years, so drawn to 17
#   Traub-Kichen (1983)   H^2 x 1.65 / 1000 below 152.4 cm (5 ft); from 5 ft, 39 kg (boys) or
#                         42.2 kg (girls) + 2.27 kg per inch over 60
#   BMI method            50th-percentile BMI for age x (H / 100)^2
#   McLaren               50th-percentile weight at the height age, the age at which H is the
#                         median height. Drawn only while the height age is at most 18: the median
#                         height curve flattens in the late teens, so near the adult median a few
#                         cm move the height age by years and the weight kinks upward, and above
#                         the adult median (P90 teenagers) there is no height age at all
#   Moore                 weight for age at the child's own height percentile (P10 height -> P10
#                         weight)
#   Devine (ARDSNet)      50 kg (boys) or 45.5 kg (girls) + 2.3 kg per inch over 60; defined
#                         from 152.4 cm (5 ft). Below 5 ft it falls toward zero (girls under about
#                         102 cm get a negative weight), yet PICU attendings report applying it
#                         there. The line slides draw it only from 5 ft; the heatmaps and the
#                         Devine slide also show it extrapolated below 5 ft, where it is positive
#
# Sources: Traub SL, Johnson CE, Am J Hosp Pharm 1980;37:195-201 (exponential form); Traub SL,
# Kichen L, Am J Hosp Pharm 1983;40:107-110 (squared-height form; the from-5-ft constants as
# tabulated in clinPK::calc_ibw and Lexicomp); McLaren DS, Read WW, Lancet 1972;2:146-148
# (height age); Moore DJ et al., Nutr Res 1985;5:797-799 (same percentile); Le Grange D et al.,
# Pediatrics 2012;129:e438-46 (BMI method, and a comparison of the three chart-based methods).
# Traub-Kichen steps at 152.4 cm, where it switches equation; that step is the formula's own.
#
# At the 50th height percentile McLaren and Moore both reduce to the median weight for age and
# nearly coincide; they are drawn with different line types so both stay visible.
#
# The CDC height percentiles are not race-specific, so within a sex, age and height percentile
# every race has the same height and the same PBW: race enters only through the GLI coefficient.
#
# Reference data: reference-data/cdc_growth_2000/ (CDC 2000 growth charts; see its README).
# Writes figures/pediatric_pbw/: the grid as CSV, then each figure as PNG (300 dpi) and PDF, 16:9.
# Colours: Okabe-Ito for formulas and races.
#
# Usage: uvr run code/tools/pediatric_pbw_vs_pfvc.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here); library(rspiro); library(ragg); library(patchwork) })

out_dir <- here("figures", "pediatric_pbw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

VT_ML_PER_KG  <- 6
# The cut-off on every slide: the mean strain delivered in the ARMA low tidal volume arm, on each
# lung-volume scale, read from the per-arm summary that code/tools/calc_external_pfvc.R writes from
# the ARMA data (aggregates only). Heatmaps colour only cells above it; line slides draw it dashed.
ARMA_SUMMARY <- Sys.getenv("PBWPFVC_ARMA_SUMMARY", here("arma_pfvc_data_pfvc_by_arm.csv"))
ARMA_ARM <- "LTV"
arma_summary <- read_csv(ARMA_SUMMARY, show_col_types = FALSE) %>% filter(arm == ARMA_ARM)
arma_mean <- function(metric) {
  value <- arma_summary$mean[arma_summary$metric == metric]
  if (length(value) != 1)
    stop(metric, " is not in ", ARMA_SUMMARY, ". Rerun code/tools/calc_external_pfvc.R on the ARMA data",
         " (it now reports vttlc and vtva) and point PBWPFVC_ARMA_SUMMARY at its output.")
  round(value, 1)
}
ARMA_LTV_MEAN_VTPFVC <- arma_mean("vtpfvc")
DEVINE_MIN_HEIGHT_CM <- 152.4
TRAUB_JOHNSON_MAX_AGE <- 17
MCLAREN_MAX_HEIGHT_AGE <- 18
AGE_GRID <- seq(3, 18, by = 0.1)
HEIGHT_PERCENTILES <- c("P10", "P50", "P90")

GLI_ETHNICITY <- c("Caucasian" = 1, "African-American" = 2, "North East Asian" = 3,
                   "South East Asian" = 4, "Other/mixed" = 5)
FORMULA_LEVELS <- c("Traub-Johnson", "Traub-Kichen", "BMI method", "McLaren", "Moore", "Devine (ARDSNet)")
OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
FORMULA_COLOURS   <- set_names(OKABE_ITO[c(1, 2, 3, 5, 7, 6)], FORMULA_LEVELS)
FORMULA_LINETYPES <- set_names(c("solid", "solid", "solid", "solid", "22", "solid"), FORMULA_LEVELS)
RACE_COLOURS <- set_names(OKABE_ITO[c(8, 6, 2, 3, 1)], names(GLI_ETHNICITY))
# GLI's African-American and South East Asian female FVC nearly coincide; dashing one keeps both visible
RACE_LINETYPES <- set_names(c("solid", "solid", "solid", "22", "solid"), names(GLI_ETHNICITY))

# =============================================================================
# CDC 2000 reference curves, interpolated to the age grid
# =============================================================================
# bmiagerev.csv repeats its header row between the boys' and the girls' rows; it is dropped here
read_cdc <- function(file) read_csv(here("reference-data", "cdc_growth_2000", file),
                                    col_types = cols(.default = "c")) %>%
  filter(Sex != "Sex") %>%
  type_convert(col_types = cols(.default = "d")) %>%
  mutate(sex = if_else(Sex == 1, "Male", "Female"))
stature_for_age <- read_cdc("statage.csv")
weight_for_age  <- read_cdc("wtage.csv")
bmi_for_age     <- read_cdc("bmiagerev.csv")

# Value of one chart column at each age (years), by linear interpolation in months
chart_at_age <- function(chart, chart_sex, column, age_years) {
  sex_rows <- filter(chart, sex == chart_sex)
  approx(sex_rows$Agemos, sex_rows[[column]], xout = age_years * 12)$y
}

# Height age: the age (years) at which the median height equals height_cm. Median height
# rises monotonically over 2-20 years, so the inverse is an interpolation with axes swapped;
# outside the chart (taller than the adult median) it is NA.
height_age_years <- function(height_cm, chart_sex) {
  sex_rows <- filter(stature_for_age, sex == chart_sex)
  approx(sex_rows$P50, sex_rows$Agemos, xout = height_cm)$y / 12
}

# =============================================================================
# The grid: sex x age x height percentile, each formula's PBW, then x GLI ethnicity
# =============================================================================
inches_over_60 <- function(height_cm) height_cm / 2.54 - 60
devine_pbw <- function(height_cm, sex) if_else(sex == "Male", 50, 45.5) + 2.3 * inches_over_60(height_cm)

pbw_grid <- expand_grid(sex = c("Male", "Female"), height_percentile = HEIGHT_PERCENTILES, age = AGE_GRID) %>%
  group_by(sex, height_percentile) %>%
  mutate(
    height_cm   = chart_at_age(stature_for_age, first(sex), first(height_percentile), age),
    median_bmi  = chart_at_age(bmi_for_age, first(sex), "P50", age),
    moore_kg    = chart_at_age(weight_for_age, first(sex), first(height_percentile), age),
    height_age  = height_age_years(height_cm, first(sex)),
    mclaren_kg  = if_else(height_age <= MCLAREN_MAX_HEIGHT_AGE,
                          chart_at_age(weight_for_age, first(sex), "P50", height_age), NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    `Traub-Johnson` = if_else(age <= TRAUB_JOHNSON_MAX_AGE, 2.396 * exp(0.01863 * height_cm), NA_real_),
    `Traub-Kichen`  = if_else(height_cm < DEVINE_MIN_HEIGHT_CM, height_cm^2 * 1.65 / 1000,
                              if_else(sex == "Male", 39, 42.2) + 2.27 * inches_over_60(height_cm)),
    `BMI method`    = median_bmi * (height_cm / 100)^2,
    McLaren         = mclaren_kg,
    Moore           = moore_kg,
    `Devine (ARDSNet)` = if_else(height_cm >= DEVINE_MIN_HEIGHT_CM, devine_pbw(height_cm, sex), NA_real_)
  ) %>%
  select(sex, height_percentile, age, height_cm, height_age, all_of(FORMULA_LEVELS)) %>%
  pivot_longer(all_of(FORMULA_LEVELS), names_to = "formula", values_to = "pbw_kg")

# pred_GLI warns on heights under 1 m (3-year-olds); GLI-2012 covers ages 3-95 and these heights
gli_fvc_litres <- function(age, height_cm, sex, ethnicity) suppressWarnings(
  pred_GLI(age = age, height = height_cm / 100, gender = if_else(sex == "Male", 1, 2),
           ethnicity = ethnicity, param = "FVC"))
# GLI static lung volumes (TLC, ages 5-80) and alveolar volume (VA, from the TLCO equations, ages
# 5-90), as in script 03. Both are European-only, with no ethnicity term; NA outside the age range.
gli_tlc_litres <- function(age, height_cm, sex) suppressWarnings(
  pred_GLIvol(age = age, height = height_cm / 100, gender = if_else(sex == "Male", 1, 2), param = "TLC"))
gli_va_litres <- function(age, height_cm, sex) suppressWarnings(
  pred_GLIdiff(age = age, height = height_cm / 100, gender = if_else(sex == "Male", 1, 2), param = "VA"))

dosing_grid <- pbw_grid %>%
  expand_grid(race = names(GLI_ETHNICITY)) %>%
  mutate(pfvc_l = gli_fvc_litres(age, height_cm, sex, GLI_ETHNICITY[race]),
         vt_ml = VT_ML_PER_KG * pbw_kg,
         vtpfvc_at_6_ml_kg = vt_ml / (pfvc_l * 1000) * 100,
         # Devine at any height, as applied below 5 ft in practice; NA where it is not positive
         devine_below_5ft = height_cm < DEVINE_MIN_HEIGHT_CM,
         devine_any_height_kg = if_else(devine_pbw(height_cm, sex) > 0, devine_pbw(height_cm, sex), NA_real_),
         vtpfvc_devine_any_height = VT_ML_PER_KG * devine_any_height_kg / (pfvc_l * 1000) * 100,
         ptlc_l = gli_tlc_litres(age, height_cm, sex),
         pva_l = gli_va_litres(age, height_cm, sex),
         vttlc_at_6_ml_kg = vt_ml / (ptlc_l * 1000) * 100,
         vtva_at_6_ml_kg = vt_ml / (pva_l * 1000) * 100,
         vttlc_devine_any_height = VT_ML_PER_KG * devine_any_height_kg / (ptlc_l * 1000) * 100,
         vtva_devine_any_height = VT_ML_PER_KG * devine_any_height_kg / (pva_l * 1000) * 100,
         sex = factor(sex, c("Male", "Female"), c("Boys", "Girls")),
         formula = factor(formula, FORMULA_LEVELS),
         race = factor(race, names(GLI_ETHNICITY)),
         height_percentile = factor(height_percentile, HEIGHT_PERCENTILES,
                                    c("10th percentile", "50th percentile", "90th percentile")))
if (anyNA(dosing_grid$pfvc_l)) stop("pred_GLI returned NA for some grid rows")
write_csv(dosing_grid, file.path(out_dir, "pediatric_pbw_pfvc_grid.csv"))

# Adults: Devine PBW at each decade's midpoint, heights fixed at the CDC age-20 percentiles
ADULT_AGES <- seq(25, 85, by = 10)
adult_grid <- expand_grid(sex = c("Male", "Female"), height_percentile = HEIGHT_PERCENTILES,
                          age = ADULT_AGES, race = names(GLI_ETHNICITY)) %>%
  mutate(height_cm = map2_dbl(sex, height_percentile, ~ chart_at_age(stature_for_age, .x, .y, 20)),
         pbw_kg = devine_pbw(height_cm, sex),
         pfvc_l = gli_fvc_litres(age, height_cm, sex, GLI_ETHNICITY[race]),
         vtpfvc_at_6_ml_kg = VT_ML_PER_KG * pbw_kg / (pfvc_l * 1000) * 100,
         ptlc_l = gli_tlc_litres(age, height_cm, sex),
         pva_l = gli_va_litres(age, height_cm, sex),
         vttlc_at_6_ml_kg = VT_ML_PER_KG * pbw_kg / (ptlc_l * 1000) * 100,
         vtva_at_6_ml_kg = VT_ML_PER_KG * pbw_kg / (pva_l * 1000) * 100,
         decade = factor(sprintf("%d-%d", age - 5, age + 4), sprintf("%d-%d", ADULT_AGES - 5, ADULT_AGES + 4)),
         sex = factor(sex, c("Male", "Female"), c("Men", "Women")),
         race = factor(race, names(GLI_ETHNICITY)),
         height_percentile = factor(height_percentile, HEIGHT_PERCENTILES,
                                    c("10th percentile", "50th percentile", "90th percentile")))
if (anyNA(adult_grid$pfvc_l)) stop("pred_GLI returned NA for some adult grid rows")
write_csv(adult_grid, file.path(out_dir, "adult_devine_pfvc_grid.csv"))

# =============================================================================
# Figures: one axis of variation per slide
# =============================================================================
# Every slide varies one input (formula, height percentile, sex or GLI ethnicity) and holds
# the rest at the reference child below. Each has two panels: the strain itself across age,
# and its difference from the reference level in percentage points of predicted FVC. The
# difference removes the shared fall with age, which otherwise fills the y axis and hides
# the comparison. Height and sex act through the formula, so their difference panels are
# drawn for every formula. The last slide reads all four axes at selected ages.
REFERENCE <- list(formula = "BMI method", height_percentile = "50th percentile", sex = "Boys", race = "Caucasian")
SELECTED_AGES <- c(5, 10, 15)
HEATMAP_AGES <- c(3, 5, 8, 10, 15)   # one heatmap per age, beside the adult heatmap
# The children the compounding slides follow from the reference, one input at a time: the
# highest-strain children in the grid are short girls under Traub-Kichen before about age 6
# and tall girls under Moore after it, with the African-American GLI equation in both
EXAMPLE_CHILDREN <- list(
  list(formula = "Traub-Kichen", height_percentile = "10th percentile", sex = "Girls", race = "African-American"),
  list(formula = "Moore", height_percentile = "90th percentile", sex = "Girls", race = "African-American")
)
EXAMPLE_COLOURS <- c("#D55E00", "#CC79A7")

HEIGHT_COLOURS <- set_names(OKABE_ITO[c(6, 8, 5)], levels(dosing_grid$height_percentile))
SEX_COLOURS    <- c(Boys = "#0072B2", Girls = "#D55E00")
VTPFVC_LABEL   <- sprintf("VT/PFVC at %d mL/kg PBW (%%)", VT_ML_PER_KG)
DIFFERENCE_LABEL <- "Difference (percentage points of PFVC)"
ARMA_LINE_NOTE <- sprintf("Dashed: ARMA low-VT arm mean, %.1f%% of PFVC", ARMA_LTV_MEAN_VTPFVC)
arma_line <- geom_hline(yintercept = ARMA_LTV_MEAN_VTPFVC, linetype = 2, colour = "grey30")

theme_set(theme_minimal(base_size = 18) +
            theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold"), plot.title.position = "plot",
                  plot.caption = element_text(size = 11, colour = "grey30")))
SOURCE_NOTE <- paste("Race-specific GLI-2012 predicted FVC. Heights: CDC 2000 stature-for-age, which is not race-specific,",
                     "so race changes PFVC only.")

save_slide <- function(plot, name, width = 13.33, height = 7.5) {
  ggsave(file.path(out_dir, paste0(name, ".png")), plot, device = agg_png, width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot, width = width, height = height)
}
held_constant <- function(varied) {
  held <- REFERENCE[setdiff(names(REFERENCE), varied)]
  paste0("Held at: ", paste(unlist(held), collapse = ", "), if ("race" %in% names(held)) " GLI equation" else "")
}
# The grid with every input except `varied` (and `also_varied`) at its reference level
hold_reference <- function(varied, also_varied = character()) {
  held <- setdiff(names(REFERENCE), c(varied, also_varied))
  reduce(held, function(data, column) filter(data, .data[[column]] == REFERENCE[[column]]), .init = dosing_grid)
}
# Strain minus the reference level's strain, at each age within every other grouping
difference_from_reference <- function(data, varied) {
  data %>%
    group_by(across(c(age, all_of(setdiff(c("formula", "height_percentile", "sex", "race"), varied))))) %>%
    mutate(difference = vtpfvc_at_6_ml_kg - vtpfvc_at_6_ml_kg[.data[[varied]] == REFERENCE[[varied]]]) %>%
    ungroup()
}
colour_scales <- function(colours, linetypes = set_names(rep("solid", length(colours)), names(colours)),
                          legend_rows = 1) list(
  scale_colour_manual(values = colours, name = NULL),
  scale_linetype_manual(values = linetypes, name = NULL),
  guides(colour = guide_legend(nrow = legend_rows), linetype = guide_legend(nrow = legend_rows))
)
strain_panel <- function(data, varied, scales, title) {
  ggplot(data, aes(age, vtpfvc_at_6_ml_kg, colour = .data[[varied]], linetype = .data[[varied]])) +
    arma_line +
    geom_line(linewidth = 1.3, na.rm = TRUE) + scales +
    scale_x_continuous(breaks = seq(3, 18, by = 3)) +
    labs(title = title, x = "Age (years)", y = VTPFVC_LABEL)
}
difference_panel <- function(data, varied, scales, title) {
  ggplot(filter(data, .data[[varied]] != REFERENCE[[varied]]),
         aes(age, difference, colour = .data[[varied]], linetype = .data[[varied]])) +
    geom_hline(yintercept = 0, colour = "grey30") +
    geom_line(linewidth = 1.3, na.rm = TRUE) + scales +
    guides(colour = "none", linetype = "none") +   # the strain panel, which has every level, carries the legend
    scale_x_continuous(breaks = seq(3, 18, by = 3)) +
    labs(title = title, x = "Age (years)", y = DIFFERENCE_LABEL)
}
axis_slide <- function(left, right, title, varied, widths = c(1, 1)) {
  (left | right) + plot_layout(widths = widths, guides = "collect") +
    plot_annotation(title = title, subtitle = paste0(held_constant(varied), ". ", ARMA_LINE_NOTE),
                    caption = SOURCE_NOTE) &
    theme(legend.position = "bottom")
}

# 1. The ingredients: PBW by formula and predicted FVC by race, boys at the median height
ingredients <- filter(dosing_grid, height_percentile == REFERENCE$height_percentile, sex == REFERENCE$sex)
pbw_panel <- ingredients %>% filter(race == REFERENCE$race) %>%
  ggplot(aes(age, pbw_kg, colour = formula, linetype = formula)) +
  geom_line(linewidth = 1.3, na.rm = TRUE) +
  colour_scales(FORMULA_COLOURS, FORMULA_LINETYPES, legend_rows = 3) +
  labs(title = "PBW by formula", x = "Age (years)", y = "PBW (kg)")
pfvc_panel <- ingredients %>% distinct(age, race, pfvc_l) %>%
  ggplot(aes(age, pfvc_l, colour = race, linetype = race)) +
  geom_line(linewidth = 1.3) +
  colour_scales(RACE_COLOURS, RACE_LINETYPES, legend_rows = 3) +
  labs(title = "Predicted FVC by GLI-2012 ethnicity", x = "Age (years)", y = "Predicted FVC (L)")
save_slide((pbw_panel | pfvc_panel) &
             plot_annotation(title = "The two ingredients: boys at the 50th percentile of height for age", caption = SOURCE_NOTE) &
             theme(legend.position = "bottom"),
           "01_ingredients")

# 2. Formula
formula_data <- hold_reference("formula") %>% difference_from_reference("formula")
formula_scales <- colour_scales(FORMULA_COLOURS, FORMULA_LINETYPES, legend_rows = 1)
save_slide(axis_slide(strain_panel(formula_data, "formula", formula_scales, "Strain at 6 mL/kg"),
                      difference_panel(formula_data, "formula", formula_scales, paste("Minus the", REFERENCE$formula)),
                      "Formula: the same child, six PBWs", "formula"),
           "02_formula")

# 3. Height percentile: strain for the reference formula, difference from the median-height child for every formula
height_scales <- colour_scales(HEIGHT_COLOURS)
height_strain <- hold_reference("height_percentile")
height_difference <- hold_reference("height_percentile", "formula") %>% difference_from_reference("height_percentile")
save_slide(axis_slide(strain_panel(height_strain, "height_percentile", height_scales, paste("Strain at 6 mL/kg,", REFERENCE$formula)),
                      difference_panel(height_difference, "height_percentile", height_scales, "Minus the median-height child, by formula") +
                        facet_wrap(~ formula, nrow = 2),
                      "Height: short and tall children of the same age", "height_percentile", widths = c(1, 1.6)),
           "03_height")

# 4. Sex, at the same height percentile (so girls are the shorter at most ages)
sex_scales <- colour_scales(SEX_COLOURS)
sex_strain <- hold_reference("sex")
sex_difference <- hold_reference("sex", "formula") %>% difference_from_reference("sex")
save_slide(axis_slide(strain_panel(sex_strain, "sex", sex_scales, paste("Strain at 6 mL/kg,", REFERENCE$formula)),
                      difference_panel(sex_difference, "sex", sex_scales, "Girls minus boys, by formula") +
                        facet_wrap(~ formula, nrow = 2),
                      "Sex: girls and boys at the same height percentile", "sex", widths = c(1, 1.6)),
           "04_sex")

# 5. GLI ethnicity. Race changes PFVC only, so its proportional effect is the same under every formula
race_scales <- colour_scales(RACE_COLOURS, RACE_LINETYPES, legend_rows = 1)
race_data <- hold_reference("race") %>% difference_from_reference("race")
save_slide(axis_slide(strain_panel(race_data, "race", race_scales, "Strain at 6 mL/kg"),
                      difference_panel(race_data, "race", race_scales, paste("Minus the", REFERENCE$race, "equation")),
                      "GLI ethnicity: the same PBW, five predicted FVCs", "race"),
           "05_race")

# 5b. GLI ethnicity, boys and girls side by side: the strain panel alone, one panel per sex, one
# slide per formula. One y axis across the formulas'
# slides, so they compare directly. Devine is drawn from 5 ft only, as on the other line slides.
race_by_sex_data <- hold_reference(c("race", "sex", "formula"))
race_by_sex_limits <- range(race_by_sex_data$vtpfvc_at_6_ml_kg, ARMA_LTV_MEAN_VTPFVC, na.rm = TRUE)
walk(levels(dosing_grid$formula), function(chosen_formula) {
  race_by_sex <- race_by_sex_data %>% filter(formula == chosen_formula) %>%
    ggplot(aes(age, vtpfvc_at_6_ml_kg, colour = race, linetype = race)) +
    arma_line +
    geom_line(linewidth = 1.3, na.rm = TRUE) +
    facet_wrap(~ sex) + race_scales +
    scale_x_continuous(breaks = seq(3, 18, by = 3), limits = range(AGE_GRID)) +
    scale_y_continuous(limits = race_by_sex_limits) +
    labs(title = paste0(chosen_formula, ": the same PBW, five predicted FVCs by GLI ethnicity"),
         subtitle = paste0("Held at: ", REFERENCE$height_percentile, ", ", chosen_formula, ". ", ARMA_LINE_NOTE),
         x = "Age (years)", y = VTPFVC_LABEL,
         caption = paste(SOURCE_NOTE, if (chosen_formula == "Devine (ARDSNet)") "Devine drawn from 5 ft (152.4 cm)." else ""))
  slug <- str_replace_all(str_to_lower(str_remove(chosen_formula, " \\(.*")), "[^a-z]+", "_")
  save_slide(race_by_sex, paste0("05b_race_by_sex_", slug))
})

# 6. All four axes at selected ages: each row varies one input from the reference child
axis_rows <- bind_rows(
  hold_reference("formula") %>% mutate(axis = "Formula", level = as.character(formula)),
  hold_reference("height_percentile") %>% mutate(axis = "Height", level = as.character(height_percentile)),
  hold_reference("sex") %>% mutate(axis = "Sex", level = as.character(sex)),
  hold_reference("race") %>% mutate(axis = "GLI ethnicity", level = as.character(race))
) %>%
  filter(round(age, 1) %in% SELECTED_AGES, !is.na(vtpfvc_at_6_ml_kg)) %>%
  mutate(axis = factor(axis, c("Formula", "Height", "Sex", "GLI ethnicity")),
         level = fct_inorder(level),
         is_reference = formula == REFERENCE$formula & height_percentile == REFERENCE$height_percentile &
                        sex == REFERENCE$sex & race == REFERENCE$race,
         age_label = factor(paste("Age", age), paste("Age", SELECTED_AGES)))
write_csv(select(axis_rows, age, axis, level, is_reference, height_cm, pbw_kg, pfvc_l, vtpfvc_at_6_ml_kg),
          file.path(out_dir, "pediatric_pbw_axis_summary.csv"))
axis_summary <- ggplot(axis_rows, aes(vtpfvc_at_6_ml_kg, fct_rev(level))) +
  geom_vline(data = filter(axis_rows, is_reference) %>% distinct(age_label, vtpfvc_at_6_ml_kg),
             aes(xintercept = vtpfvc_at_6_ml_kg), colour = "grey60") +
  geom_vline(xintercept = ARMA_LTV_MEAN_VTPFVC, linetype = 2, colour = "grey30") +
  geom_point(aes(colour = axis, shape = is_reference), size = 4) +
  facet_grid(axis ~ age_label, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = set_names(OKABE_ITO[c(1, 3, 5, 6)], levels(axis_rows$axis)), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16), guide = "none") +
  labs(title = "Which input moves the strain most? One input at a time, at ages 5, 10 and 15",
       subtitle = paste0("Open circle and grey line: the reference child (", paste(unlist(REFERENCE), collapse = ", "), "). ", ARMA_LINE_NOTE),
       x = VTPFVC_LABEL, y = NULL, caption = SOURCE_NOTE) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing.y = unit(1, "lines"))
save_slide(axis_summary, "06_axis_summary")

# =============================================================================
# Compounding: the inputs' errors add up in one child
# =============================================================================
child_label <- function(child) paste0(child$race, " ", str_to_lower(child$sex), ", ",
                                      str_remove(child$height_percentile, " percentile"), " height percentile, ", child$formula)
child_rows <- function(child) filter(dosing_grid, formula == child$formula, sex == child$sex,
                                     height_percentile == child$height_percentile, race == child$race)
strain_of <- function(child, at_age) {
  row <- filter(child_rows(child), round(age, 1) == at_age)
  stopifnot(nrow(row) == 1)
  row$vtpfvc_at_6_ml_kg
}

# 7. Waterfall: from the reference child to each example child, changing one input per step in
# the order formula, sex, height, race. Each bar is the change in strain that input makes given
# the inputs already changed, so the bars sum to the child's total excess over the reference
# child; another order splits the total differently but reaches the same child.
STEP_ORDER <- c("formula", "sex", "height_percentile", "race")
STEP_NAMES <- c(formula = "Formula", sex = "Sex", height_percentile = "Height", race = "Race")
waterfall <- map_dfr(EXAMPLE_CHILDREN, function(child) map_dfr(SELECTED_AGES, function(at_age) {
  states <- accumulate(STEP_ORDER, function(state, input) { state[[input]] <- child[[input]]; state }, .init = REFERENCE)
  excess <- map_dbl(states, strain_of, at_age = at_age) - strain_of(REFERENCE, at_age)
  child_strain <- strain_of(child, at_age)   # before tibble(), whose `child` column would mask the list
  step_labels <- c(paste0(STEP_NAMES[STEP_ORDER], ":\n", str_replace(str_remove(unlist(child[STEP_ORDER]), " percentile"), "-", "-\n")), "Total")
  tibble(child = child_label(child), age = at_age,
         reference_strain = strain_of(REFERENCE, at_age), child_strain = child_strain,
         step = step_labels,
         position = seq_len(length(STEP_ORDER) + 1),
         start = c(head(excess, -1), 0), end = c(excess[-1], last(excess)),
         kind = c(rep("Step", length(STEP_ORDER)), "Total"))
})) %>%
  mutate(child = fct_inorder(child),
         panel = fct_reorder(sprintf("Age %d: %.1f%% to %.1f%%", age, reference_strain, child_strain), age),
         fill_group = case_when(kind == "Total" ~ "Total", end >= start ~ "Raises strain", TRUE ~ "Lowers strain"),
         bar_label = sprintf("%+.1f", round(end - start, 1) + 0))   # + 0 turns a rounded -0 into 0
write_csv(select(waterfall, child, age, step, start, end), file.path(out_dir, "pediatric_pbw_waterfall.csv"))
walk2(EXAMPLE_CHILDREN, seq_along(EXAMPLE_CHILDREN), function(child, child_index) {
  steps <- filter(waterfall, child == child_label(!!child))
  plot <- ggplot(steps) +
    geom_hline(yintercept = 0, colour = "grey30") +
    geom_rect(aes(xmin = position - 0.4, xmax = position + 0.4, ymin = start, ymax = end, fill = fill_group)) +
    geom_text(aes(x = position, y = if_else(end >= start, pmax(start, end), pmin(start, end)), label = bar_label,
                  vjust = if_else(end >= start, -0.4, 1.3)), size = 5.5) +
    facet_wrap(~ panel, nrow = 1) +
    scale_x_continuous(breaks = steps$position, labels = steps$step) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    scale_fill_manual(values = c(`Raises strain` = "#D55E00", `Lowers strain` = "#0072B2", Total = "grey40"), name = NULL) +
    labs(title = paste0("The errors compound: ", child_label(child)),
         subtitle = paste0("Change from the reference child (", child_label(REFERENCE), "), one input at a time"),
         x = NULL, y = "Excess strain (percentage points of PFVC)",
         caption = paste("Panel titles: the reference child's strain, then this child's.", SOURCE_NOTE)) +
    theme(axis.text.x = element_text(size = 12), panel.grid.major.x = element_blank(), legend.position = "top")
  save_slide(plot, sprintf("07%s_compounding_waterfall", letters[child_index]))
})

# 8. Across age: the example children and the reference child inside the range of all 180 children
all_children_range <- dosing_grid %>% filter(!is.na(vtpfvc_at_6_ml_kg)) %>%
  group_by(age) %>% summarise(lowest = min(vtpfvc_at_6_ml_kg), highest = max(vtpfvc_at_6_ml_kg), .groups = "drop")
named_children <- map_dfr(c(list(REFERENCE), EXAMPLE_CHILDREN), ~ child_rows(.x) %>% mutate(child = child_label(.x))) %>%
  mutate(child = fct_inorder(child))
trajectory_plot <- ggplot() +
  geom_ribbon(data = all_children_range, aes(age, ymin = lowest, ymax = highest), fill = "grey85") +
  arma_line +
  geom_line(data = named_children, aes(age, vtpfvc_at_6_ml_kg, colour = child), linewidth = 1.5, na.rm = TRUE) +
  scale_colour_manual(values = set_names(c("#000000", EXAMPLE_COLOURS), levels(named_children$child)), name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_x_continuous(breaks = seq(3, 18, by = 3)) +
  labs(title = "Which child is most exposed changes with age",
       subtitle = paste("Grey band: all 180 children (6 formulas x 2 sexes x 3 height percentiles x 5 GLI ethnicities).", ARMA_LINE_NOTE),
       x = "Age (years)", y = VTPFVC_LABEL, caption = SOURCE_NOTE)
save_slide(trajectory_plot, "08_compounding_across_age")

# Heatmaps: every child at ages 5, 10 and 15 (rows sex x GLI ethnicity, columns height percentile
# within formula) and every adult by decade (Devine only), on three lung-volume scales: predicted
# FVC (slides 9 children, 10 adults), predicted TLC (12, 13) and predicted alveolar volume (14, 15),
# children at each of HEATMAP_AGES the reference equation covers (TLC and VA start at 5). Cells at or
# below the ARMA low tidal volume arm's mean on that scale are grey; cells above it are coloured on
# viridis "inferno", brighter for more strain, values bold. Within a scale the children and the
# adults share one colour scale. The Devine column shows the formula at every height, italic below
# 5 ft. GLI's TLC and VA equations have no ethnicity term, so on those scales the rows are sex only.
LUNG_SCALES <- list(
  pfvc = list(name = "PFVC", long = "predicted FVC", column = "vtpfvc_at_6_ml_kg", devine = "vtpfvc_devine_any_height",
              arma_metric = "vtpfvc", by_race = TRUE, min_age = 3, child_slide = 9, adult_slide = 10,
              note = "Race-specific GLI-2012 predicted FVC."),
  tlc = list(name = "TLC", long = "predicted TLC", column = "vttlc_at_6_ml_kg", devine = "vttlc_devine_any_height",
             arma_metric = "vttlc", by_race = FALSE, min_age = 5, child_slide = 12, adult_slide = 13,
             note = "GLI 2021 predicted TLC (ages 5-80, European reference, no ethnicity term)."),
  va = list(name = "VA", long = "predicted alveolar volume", column = "vtva_at_6_ml_kg", devine = "vtva_devine_any_height",
            arma_metric = "vtva", by_race = FALSE, min_age = 5, child_slide = 14, adult_slide = 15,
            note = "GLI 2017 predicted alveolar volume, VA (ages 5-90, European reference, no ethnicity term).")
)
draw_heatmaps <- function(lung) {
  cutoff <- arma_mean(lung$arma_metric)
  value_label <- sprintf("VT/%s at\n6 mL/kg PBW (%%)", lung$name)
  child_ages <- HEATMAP_AGES[HEATMAP_AGES >= lung$min_age]
  child_cells <- dosing_grid %>% filter(round(age, 1) %in% child_ages) %>%
    mutate(extrapolated = formula == "Devine (ARDSNet)" & devine_below_5ft,
           strain = if_else(formula == "Devine (ARDSNet)", .data[[lung$devine]], .data[[lung$column]]),
           formula = fct_recode(formula, `Devine (ARDSNet)*` = "Devine (ARDSNet)"))
  adult_cells <- adult_grid %>% mutate(strain = .data[[lung$column]], extrapolated = FALSE)
  if (!lung$by_race) {   # the scale has no ethnicity term: one row per sex
    child_cells <- filter(child_cells, race == "Caucasian")
    adult_cells <- filter(adult_cells, race == "Caucasian")
  }
  limits <- range(c(child_cells$strain, adult_cells$strain), na.rm = TRUE)
  # Grey to the cut-off, inferno above it; compared at the printed precision, so every cell that
  # reads the cut-off is grey (the fill switches half a printed digit above it)
  above_cutoff <- function(strain) coalesce(round(strain, 1) > cutoff, FALSE)
  white_text <- function(strain) above_cutoff(strain) & strain < (cutoff + limits[2]) / 2
  above_colours <- viridisLite::inferno(64, begin = 0.3, end = 0.98)
  breakpoint <- scales::rescale(cutoff + 0.05, from = limits)
  fill_scale <- scale_fill_gradientn(
    colours = c("grey93", "grey78", above_colours),
    values = c(0, breakpoint, seq(breakpoint + 1e-6, 1, length.out = length(above_colours))),
    limits = limits, oob = scales::squish, na.value = "white",
    breaks = sort(c(cutoff, setdiff(scales::breaks_pretty(4)(limits), round(cutoff)))),
    labels = function(x) if_else(x == cutoff, paste0(x, " ARMA mean"), as.character(x)),
    name = value_label, guide = guide_colourbar(barheight = unit(14, "lines")))
  style_cells <- function(cells) cells %>%
    mutate(row_label = if (lung$by_race) fct_rev(interaction(race, sex, sep = ", ", lex.order = TRUE)) else fct_rev(sex),
           height_short = factor(str_remove(height_percentile, "th percentile"), c("10", "50", "90")),
           dark_cell = white_text(strain),
           label_face = case_when(above_cutoff(strain) & extrapolated ~ "bold.italic", above_cutoff(strain) ~ "bold",
                                  extrapolated ~ "italic", TRUE ~ "plain"))
  subtitle <- sprintf("Coloured, bold: above the mean delivered in the ARMA low tidal volume arm (%.1f%% of %s). Brighter: more strain",
                      cutoff, lung$long)
  heatmap <- function(cells, facet, x_label, title, caption, text_size) {
    ggplot(cells, aes(height_short, row_label, fill = strain)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = if_else(is.na(strain), "", sprintf("%.1f", strain)), colour = dark_cell, fontface = label_face),
                size = text_size, na.rm = TRUE) +
      facet_grid(cols = vars({{ facet }}), labeller = label_wrap_gen(12)) +
      fill_scale +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
      labs(title = title, subtitle = subtitle, x = x_label, y = NULL, caption = str_wrap(caption, 150)) +
      theme(panel.grid = element_blank(), legend.position = "right", legend.title = element_text(size = 12),
            plot.subtitle = element_text(size = 14), axis.text.y = element_text(size = 13), strip.text = element_text(size = 13))
  }
  walk(child_ages, function(at_age) {
    cells <- child_cells %>% filter(round(age, 1) == at_age) %>%
      filter(!all(is.na(strain)), .by = formula) %>%   # drop a formula undefined for every child at this age
      style_cells()
    plot <- heatmap(cells, formula, "Height percentile",
                    sprintf("Every child at age %d: VT/%s at %d mL/kg PBW", at_age, lung$name, VT_ML_PER_KG),
                    paste(lung$note, "Heights: CDC 2000 stature-for-age, not race-specific.",
                          "*Devine below 5 ft (italic) is extrapolated, as some PICUs apply it.",
                          "Same colour scale at ages", paste(child_ages, collapse = ", "), "and in adults.",
                          "Blank: formula or reference equation not defined, or Devine weight not positive."), 3.8)
    save_slide(plot, sprintf("%02d%s_%s_every_child_age_%02d", lung$child_slide, letters[match(at_age, child_ages)],
                             tolower(lung$name), at_age))
  })
  adult_plot <- heatmap(style_cells(adult_cells), decade, "Height percentile (CDC, age 20)",
                        sprintf("Every adult, by decade: VT/%s at %d mL/kg Devine PBW", lung$name, VT_ML_PER_KG),
                        paste(lung$note, "Each decade's midpoint age. Heights held at the CDC 2000 age-20 10th, 50th",
                              "and 90th percentiles at every age, so decades differ only through the reference",
                              "equation's age term. Same colour scale as the pediatric slides.",
                              "Blank: reference equation not defined at this age."), 3.3) +
    theme(axis.text.x = element_text(size = 11))
  save_slide(adult_plot, sprintf("%02d_%s_every_adult_by_decade", lung$adult_slide, tolower(lung$name)))
}
draw_heatmaps(LUNG_SCALES$pfvc)

# 11. Devine below 5 ft: the ARDSNet formula against McLaren (recommended for pediatric ARDS by
# Ward et al., Pediatr Crit Care Med 2018), Caucasian GLI equation. Points mark 5 ft.
devine_children <- dosing_grid %>%
  filter(race == REFERENCE$race, formula == "McLaren") %>%
  transmute(sex, height_percentile, age, height_cm, devine_below_5ft,
            `Devine (ARDSNet)` = devine_any_height_kg, McLaren = pbw_kg,
            pfvc_l) %>%
  pivot_longer(c(`Devine (ARDSNet)`, McLaren), names_to = "formula", values_to = "pbw_kg") %>%
  mutate(vtpfvc = VT_ML_PER_KG * pbw_kg / (pfvc_l * 1000) * 100)
five_feet <- devine_children %>% filter(formula == "Devine (ARDSNet)", !devine_below_5ft) %>%
  slice_min(age, by = c(sex, height_percentile))
devine_panel <- function(y, y_label, title) {
  ggplot(devine_children, aes(age, .data[[y]], colour = height_percentile, linetype = formula)) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(data = five_feet, size = 3, show.legend = FALSE) +
    facet_wrap(~ sex) +
    scale_colour_manual(values = HEIGHT_COLOURS, name = NULL) +
    scale_linetype_manual(values = c(`Devine (ARDSNet)` = "solid", McLaren = "22"), name = NULL) +
    scale_x_continuous(breaks = seq(3, 18, by = 3)) +
    labs(title = title, x = "Age (years)", y = y_label)
}
devine_slide <- (devine_panel("pbw_kg", "PBW (kg)", "Predicted body weight") |
                   (devine_panel("vtpfvc", VTPFVC_LABEL, "Strain at 6 mL/kg") +
                      arma_line)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Devine below 5 ft: too little volume in young children, too much approaching 5 ft",
                  subtitle = "Devine (ARDSNet) at any height against McLaren. Points: 5 ft (152.4 cm). Caucasian GLI equation",
                  caption = str_wrap(paste("Devine below 5 ft is extrapolated; girls under about 102 cm get a negative weight, not drawn.",
                                           SOURCE_NOTE), 160)) &
  theme(legend.position = "bottom")
save_slide(devine_slide, "11_devine_below_5ft")


# TLC and alveolar volume heatmaps last: they need the ARMA means on those scales
draw_heatmaps(LUNG_SCALES$tlc)
draw_heatmaps(LUNG_SCALES$va)

message("Wrote ", nrow(dosing_grid), " grid rows and the slides to ", out_dir)
