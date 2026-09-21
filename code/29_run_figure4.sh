#!/bin/bash
# =============================================================================
# 29_run_figure4: every analysis behind manuscript figure 4, then the figure
# =============================================================================
# Figure 4 asks whether a smaller predicted lung, at the same VT/PBW, goes with
# organ-injury markers that worsen faster over the first 7 days of ventilation. This
# script produces it for one site, from the CLIF tables to the PDF:
#
#   markers   platelets, bilirubin, creatinine (dialysis, continuous or intermittent,
#             as a third competing cause; ESRD censored at day 0), vasopressor dose on
#             pressor days, oxygen saturation index
#   arms      ventilated, all patients           (panels A, B and C)
#             no respiratory support, unmatched  (panel B, the negative control)
#             no respiratory support, matched to the ventilated cohort's severity
#             optional: ventilated by baseline SF class (SF_BANDS; off by default)
#   each fit adjusted and unadjusted for age, sex and race
#
# Steps, each logged to output/{site}_output/logs/figure4_{stamp}/:
#   1 build     scripts 01-03 for the ventilated cohort and the no-support cohort,
#               each only if its derived tables are missing (FORCE_BUILD=1 rebuilds)
#   2 panels    the 7-day panel of both cohorts
#   3 anchors   the severity-anchor distributions (final/injury/jm_severity_anchor_*)
#   4 floors    each control marker's severity floor = the ventilated cohort's median
#               anchor score (SEV_MIN overrides, e.g. "creatinine=2,platelets=3,...")
#   5 fits      22_biotrauma_fit.R and 23_biotrauma_report.R for every arm
#   6 figure    27_control_comparison.R, then 24_biotrauma_figures.R
# Fits already on disk with the same chain settings are reused, so a rerun after a
# failure costs only what failed. A failed step is reported and the run continues.
#
# The oxygen saturation index needs positive-pressure ventilation (mean airway
# pressure), so it has no control arm. A marker or arm with too few patients or
# deaths is skipped by the fit, with the reason in its manifest.
#
# Usage (from anywhere; the script moves to the repo root):
#   caffeinate -i nohup bash code/29_run_figure4.sh > figure4.out 2>&1 &
#   bash code/29_run_figure4.sh --dry-run
# Site default (2026-09-21): 26 fits at 2000 / 500 iterations, four at a time. At
# MIMIC the divergence terms the figure rests on converged at 2000 iterations; the
# hazard blocks did not converge at any length tried.
# Knobs (environment): ITER BURNIN CHAINS THIN (2000 / 500 / 3 / 5), PAR (fits at a
#   time, 4; an earlier estimate put a 7-day fit at a 7,000-patient site at 15-25 GB,
#   so four at once can need 60-100 GB: lower PAR on a smaller machine), MARKERS, CONTROL_MARKERS, CREATININE (1; 0 skips it),
#   SF_BANDS (baseline SF classes of the ventilated cohort, off by default; the lead
#   site runs SF_BANDS="235,315 115,235 0,115"), SEV_MIN, FORCE_BUILD, FORCE_PANEL.
#   Output: final/injury/biotrauma_fig_main_pfvc_7d_{site}.pdf.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

MARKERS=${MARKERS:-osi,pressor_dose,platelets,bilirubin}   # creatinine runs on its own, with RRT as a third cause
CONTROL_MARKERS=${CONTROL_MARKERS:-pressor_dose,platelets,bilirubin}
CREATININE=${CREATININE:-1}
SF_BANDS=${SF_BANDS:-}                   # e.g. "235,315 115,235 0,115"; off by default
SEV_MIN=${SEV_MIN:-}
ITER=${ITER:-2000}; BURNIN=${BURNIN:-500}; CHAINS=${CHAINS:-3}; THIN=${THIN:-5}; PAR=${PAR:-4}
FORCE_BUILD=${FORCE_BUILD:-0}
DRY=0; [[ "${1:-}" == "--dry-run" ]] && DRY=1

# the site's name from config.json (sed: Rscript's stdout carries renv notices)
BASE_SITE=${PBWPFVC_SITE_NAME:-$(sed -n 's/.*"site_name" *: *"\([^"]*\)".*/\1/p' config/config.json | head -n 1)}
BASE_SITE="${BASE_SITE%_niv}"; BASE_SITE="${BASE_SITE%_nosupport}"
[[ -n "$BASE_SITE" ]] || { echo "could not read site_name from config/config.json"; exit 1; }
ROOT="output/${BASE_SITE}_output"
LOG_DIR="$ROOT/logs/figure4_$(date +%Y%m%d_%H%M%S)"
[[ $DRY == 0 ]] && mkdir -p "$LOG_DIR"
unset PBWPFVC_COHORT
export PBWPFVC_SITE_NAME=$BASE_SITE
export PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 PBWPFVC_JM_MODIFIER=pfvc PBWPFVC_JM_MODELS=main
export PBWPFVC_JM_ITER=$ITER PBWPFVC_JM_BURNIN=$BURNIN PBWPFVC_JM_CHAINS=$CHAINS PBWPFVC_JM_THIN=$THIN PBWPFVC_JM_PAR=$PAR
echo "site ${BASE_SITE}; markers ${MARKERS}$([[ $CREATININE == 1 ]] && echo ",creatinine"); controls ${CONTROL_MARKERS}; SF bands '${SF_BANDS:-none}'; chains ${ITER}/${BURNIN} x ${CHAINS}"
[[ $DRY == 0 ]] && echo "logs -> $LOG_DIR"

FAILED=()
run_step () {     # name, then the command
  local step_name=$1; shift
  if [[ $DRY == 1 ]]; then echo "[dry] ${step_name}: $*"; return 0; fi
  echo "[$(date +%H:%M:%S)] ${step_name}"
  if "$@" > "${LOG_DIR}/${step_name}.log" 2>&1; then echo "    ok"
  else FAILED+=("$step_name"); echo "    FAILED, see ${LOG_DIR}/${step_name}.log"
       tail -5 "${LOG_DIR}/${step_name}.log" | sed 's/^/    | /'; fi
}
with_cohort () {  # cohort (imv | nosupport), then the command
  local cohort=$1; shift
  if [[ $cohort == imv ]]; then env -u PBWPFVC_COHORT "$@"; else env PBWPFVC_COHORT=$cohort "$@"; fi
}
# one arm: the markers, then creatinine with dialysis as a third competing cause
fit_arm () {      # arm name, cohort, marker list, then extra environment assignments
  local arm=$1 cohort=$2 markers=$3; shift 3
  # a cohort whose panel failed to build is not fitted: its old panel is out of date
  if [[ " ${FAILED[*]-} " == *" panel_${cohort} "* ]]; then
    echo "[$(date +%H:%M:%S)] ${arm}: skipped, the ${cohort} panel failed to build"; FAILED+=("${arm}_skipped"); return 0
  fi
  run_step "${arm}_fit"    with_cohort "$cohort" env "$@" PBWPFVC_JM_MARKERS="$markers" Rscript code/22_biotrauma_fit.R
  run_step "${arm}_report" with_cohort "$cohort" env "$@" PBWPFVC_JM_MARKERS="$markers" Rscript code/23_biotrauma_report.R
  if [[ $CREATININE == 1 ]]; then
    run_step "${arm}_creatinine_fit"    with_cohort "$cohort" env "$@" PBWPFVC_JM_MARKERS=creatinine PBWPFVC_JM_RRT_EVENT=1 Rscript code/22_biotrauma_fit.R
    run_step "${arm}_creatinine_report" with_cohort "$cohort" env "$@" PBWPFVC_JM_MARKERS=creatinine PBWPFVC_JM_RRT_EVENT=1 Rscript code/23_biotrauma_report.R
  fi
}

# ---- 1 build
build_cohort () { # cohort, folder holding its derived tables
  local cohort=$1 derived=$2
  # a cohort built before dialysis and ESRD entered the RRT definition (2026-09-21) lacks
  # rrt_sources_available.rds, and its panel cannot be built: rebuild it
  if [[ $FORCE_BUILD == 0 && -f "$derived/analysis_cross_sectional.parquet" && -f "$derived/rrt_sources_available.rds" ]]; then
    echo "[$(date +%H:%M:%S)] ${cohort}: cohort already built, scripts 01-03 skipped (FORCE_BUILD=1 rebuilds)"; return 0
  fi
  for script in 01_cohort_identification 02_quality_checks 03_variable_derivation; do
    run_step "build_${cohort}_${script}" with_cohort "$cohort" Rscript "code/${script}.R"
  done
}
build_cohort imv       "$ROOT/intermediate"
build_cohort nosupport "$ROOT/intermediate/controls/nosupport"

# ---- 2 panels, rebuilt only when something they are built from has changed: a fit made
#      on an older panel is refitted (22_biotrauma_fit.R compares the times), so an
#      unconditional rebuild would refit everything on every rerun. FORCE_PANEL=1 rebuilds.
FORCE_PANEL=${FORCE_PANEL:-0}
build_panel () {  # cohort, folder holding its derived tables
  local cohort=$1 derived=$2 panel="$2/jm_surv_7d.parquet" dep
  if [[ $FORCE_PANEL == 0 && $DRY == 0 && -f "$panel" ]]; then
    local stale=""
    for dep in code/21_biotrauma_panel.R code/10_panel_common.R code/20_biotrauma_grid.R utils/config.R \
               "$derived/analysis_cross_sectional.parquet" "$derived/cohort_dialysis.parquet" "$derived/cohort_esrd.parquet"; do
      [[ ! -e "$dep" || "$dep" -nt "$panel" ]] && stale="$dep" && break
    done
    if [[ -z "$stale" ]]; then echo "[$(date +%H:%M:%S)] panel_${cohort}: up to date, kept (FORCE_PANEL=1 rebuilds)"; return 0; fi
    echo "[$(date +%H:%M:%S)] panel_${cohort}: $(basename "$stale") is newer than the panel (or missing); rebuilding"
  fi
  run_step "panel_${cohort}" with_cohort "$cohort" Rscript code/21_biotrauma_panel.R
}
build_panel imv       "$ROOT/intermediate"
build_panel nosupport "$ROOT/intermediate/controls/nosupport"

# ---- 3 anchors
ANCHOR_MARKERS="creatinine,${CONTROL_MARKERS}"
run_step anchors_ventilated with_cohort imv       env PBWPFVC_JM_ANCHOR_ONLY=1 PBWPFVC_JM_MARKERS="$ANCHOR_MARKERS" Rscript code/22_biotrauma_fit.R
run_step anchors_nosupport  with_cohort nosupport env PBWPFVC_JM_ANCHOR_ONLY=1 PBWPFVC_JM_MARKERS="$ANCHOR_MARKERS" Rscript code/22_biotrauma_fit.R

# ---- 4 floors: the ventilated median anchor score per marker, unless SEV_MIN is given
FLOOR_FILE="${LOG_DIR}/severity_floors.txt"
if [[ $DRY == 1 ]]; then
  echo "[dry] floors: ventilated median of final/injury/jm_severity_anchor_7d_${BASE_SITE}.csv per marker (${ANCHOR_MARKERS}), unless SEV_MIN"
  SEV_SPEC="creatinine=M,platelets=M,..."; SEV_TAG="sev_..._"
else
  SEV_MIN="$SEV_MIN" ANCHOR_MARKERS="$ANCHOR_MARKERS" FLOOR_FILE="$FLOOR_FILE" \
    ANCHOR_FILE="$ROOT/final/injury/jm_severity_anchor_7d_${BASE_SITE}.csv" Rscript -e '
    suppressMessages(library(readr))
    markers <- strsplit(Sys.getenv("ANCHOR_MARKERS"), ",")[[1]]
    spec <- Sys.getenv("SEV_MIN")
    if (!nzchar(spec)) {
      a <- read_csv(Sys.getenv("ANCHOR_FILE"), show_col_types = FALSE)
      # the median band: the highest band at or above which half the ventilated patients sit
      floors <- vapply(markers, function(m) {
        d <- a[a$marker == m, ]
        if (!nrow(d)) stop("no anchor distribution for ", m, " in ", Sys.getenv("ANCHOR_FILE"))
        max(d$sev_anchor_from[d$pct_at_or_above_from >= 50])
      }, numeric(1))
      spec <- paste0(names(floors), "=", floors, collapse = ",")
    }
    pairs <- strsplit(strsplit(spec, ",")[[1]], "=")
    floors <- setNames(as.numeric(vapply(pairs, `[`, "", 2)), vapply(pairs, `[`, "", 1))
    floors <- floors[order(names(floors))]
    tag <- paste0("sev_", paste0(names(floors), floors, collapse = "_"), "_")
    writeLines(c(spec, tag), Sys.getenv("FLOOR_FILE"))' > "${LOG_DIR}/floors.log" 2>&1
  if [[ ! -s "$FLOOR_FILE" ]]; then
    echo "floors could not be set (see ${LOG_DIR}/floors.log); the matched control arm is skipped"; FAILED+=("floors")
    SEV_SPEC=""; SEV_TAG=""
  else
    SEV_SPEC=$(sed -n 1p "$FLOOR_FILE"); SEV_TAG=$(sed -n 2p "$FLOOR_FILE")
    echo "[$(date +%H:%M:%S)] severity floors (ventilated median anchor${SEV_MIN:+, from SEV_MIN}): ${SEV_SPEC}"
  fi
fi

# ---- 5 fits, arm by arm
fit_arm ventilated imv "$MARKERS"
for BAND in $SF_BANDS; do fit_arm "ventilated_sf${BAND/,/to}" imv "$MARKERS" PBWPFVC_JM_SF_BAND="$BAND"; done
fit_arm nosupport nosupport "$CONTROL_MARKERS"
[[ -n "$SEV_SPEC" ]] && fit_arm nosupport_matched nosupport "$CONTROL_MARKERS" PBWPFVC_JM_SEV_MIN="$SEV_SPEC"

# ---- 6 comparison table and the figure
FIG_MARKERS="platelets,bilirubin$([[ $CREATININE == 1 ]] && echo ",creatinine"),pressor_dose,osi"
run_step comparison with_cohort imv Rscript code/27_control_comparison.R
run_step figure with_cohort imv env PBWPFVC_JM_WITH_RRT=1 PBWPFVC_FIG_MARKERS="$FIG_MARKERS" PBWPFVC_FIG_SEV_TAG="$SEV_TAG" \
  Rscript code/24_biotrauma_figures.R

if [[ $DRY == 0 ]]; then
  echo
  echo "figure  -> $ROOT/final/injury/biotrauma_fig_main_pfvc_7d_${BASE_SITE}.pdf"
  echo "tables  -> $ROOT/final/injury/ and $ROOT/final/controls/"
  if (( ${#FAILED[@]} )); then echo "FAILED steps (${#FAILED[@]}): ${FAILED[*]}"; exit 1; else echo "all steps ok"; fi
fi
