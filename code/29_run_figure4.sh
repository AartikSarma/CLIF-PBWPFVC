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
#             pressor days, oxygen saturation index, and the SF ratio (panel C, the
#             positive control for mechanics: a larger VT/PFVC recruits lung, so SF
#             should be better early in the smaller predicted lung, then reverse)
#   arms      ventilated, all patients           (panels A, B and C)
#             no respiratory support             (panel B, the negative control): every
#                                                patient, the divergence read at the
#                                                ventilated cohort's mean severity
#             no support, hypoxemic on the index day (index-day SF <= 315, the ventilated
#                                                cohort's own gate), read at the ventilated severity:
#                                                the arms then differ in ventilation, not
#                                                hypoxemia (27 writes its DiD separately;
#                                                HYPOXEMIC_CONTROL_MARKERS, platelets by default)
#             optional: ventilated by baseline SF class (SF_BANDS; off by default)
#   The control is standardised, not matched (2026-09-21): severity cannot confound a
#   PFVC fixed by height, age, sex and race, but it could MODIFY the divergence, so the
#   control's divergence varies with its severity anchor and is read at the ventilated
#   mean (20_biotrauma_grid.R, PBWPFVC_JM_SEV_CENTER). No control patient is discarded,
#   and the severity x divergence term tests whether sicker controls diverge faster.
#   each fit adjusted and unadjusted for age, sex and race
#
# Steps, each logged to output/{site}_output/logs/figure4_{stamp}/:
#   1 build     scripts 01-03 for the ventilated cohort and the no-support cohort,
#               each only if its derived tables are missing (FORCE_BUILD=1 rebuilds)
#   2 panels    the 7-day panel of both cohorts
#   3 anchors   the severity-anchor distributions of both cohorts, and the ventilated
#               mean anchor per marker (final/injury/jm_severity_anchor_mean_*)
#   4 centres   each control marker's centre = that ventilated mean
#   5 fits      22_biotrauma_fit.R and 23_biotrauma_report.R for every arm
#   6 figure    27_control_comparison.R (with the difference-in-differences), then
#               24_biotrauma_figures.R: figure 4 and the checks figure
#               (biotrauma_fig_checks_*: the DiD)
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
# Site default (2026-09-24): 20 figure-4 fits at 2000 / 500 iterations, four at a time. At
# MIMIC the divergence terms the figure rests on converged at 2000 iterations; the
# hazard blocks did not converge at any length tried.
# Knobs (environment): ITER BURNIN CHAINS THIN (2000 / 500 / 3 / 5), PAR (fits at a
#   time, 4; an earlier estimate put a 7-day fit at a 7,000-patient site at 15-25 GB,
#   so four at once can need 60-100 GB: lower PAR on a smaller machine), MARKERS, CONTROL_MARKERS, CREATININE (1; 0 skips it),
#   SF_BANDS (baseline SF classes of the ventilated cohort, off by default; the lead
#   site runs SF_BANDS="235,315 115,235 0,115"), HYPOXEMIC_CONTROL_MARKERS (the hypoxemic
#   control arm, platelets by default, creatinine not fitted there; empty skips it),
#   FORCE_BUILD, FORCE_PANEL. The VT/PFVC companion (form vtpfvc) was dropped on
#   2026-09-24: at a fixed VT/PBW, VT/PFVC moves only with PBW/PFVC, whose variance
#   after age, sex and race is a few percent (Claim 5a), so the companion re-reads the
#   pfvc form on a scale the data cannot identify. 22-24 still accept the form.
#   Output: final/injury/biotrauma_fig_main_pfvc_7d_{site}.pdf.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

MARKERS=${MARKERS:-osi,sf,pressor_dose,platelets,bilirubin}   # creatinine runs on its own, with RRT as a third cause
CONTROL_MARKERS=${CONTROL_MARKERS:-pressor_dose,platelets,bilirubin}
HYPOXEMIC_CONTROL_MARKERS=${HYPOXEMIC_CONTROL_MARKERS-platelets}   # the hypoxemic control arm; set empty to skip
CREATININE=${CREATININE:-1}
SF_BANDS=${SF_BANDS:-}                   # e.g. "235,315 115,235 0,115"; off by default
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
  # rrt_sources_available.rds, and its panel cannot be built: rebuild it. A control built
  # before it was indexed at ICU admission (same day) lacks cohort_icu_stays.parquet.
  local icu_ok=1
  [[ $cohort == nosupport && ! -f "$derived/cohort_icu_stays.parquet" ]] && icu_ok=0
  if [[ $FORCE_BUILD == 0 && $icu_ok == 1 && -f "$derived/analysis_cross_sectional.parquet" && -f "$derived/rrt_sources_available.rds" ]]; then
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

# ---- 4 centres: the ventilated cohort's mean anchor per control marker, where each
#      control fit reads its divergence (PBWPFVC_JM_SEV_CENTER, 20_biotrauma_grid.R)
CENTER_FILE="$ROOT/final/injury/jm_severity_anchor_mean_7d_${BASE_SITE}.csv"
if [[ $DRY == 1 ]]; then
  echo "[dry] centres: ventilated mean anchor per marker (${ANCHOR_MARKERS}) from ${CENTER_FILE}"
  SEV_CENTER="creatinine=M,platelets=M,..."
else
  # the file is written by anchors_ventilated; one "marker=mean" pair per control marker
  SEV_CENTER=$(awk -F, -v want=",${ANCHOR_MARKERS}," 'NR > 1 && index(want, "," $1 ",") { printf "%s%s=%.4f", sep, $1, $2; sep = "," }' "$CENTER_FILE" 2>/dev/null)
  n_found=$(tr ',' '\n' <<< "$SEV_CENTER" | grep -c '=')
  n_want=$(tr ',' '\n' <<< "$ANCHOR_MARKERS" | grep -c .)
  if [[ $n_found != "$n_want" ]]; then
    echo "severity centres missing for some of ${ANCHOR_MARKERS} (found '${SEV_CENTER}' in ${CENTER_FILE}); the control arm is skipped"
    FAILED+=("centres"); SEV_CENTER=""
  else
    echo "[$(date +%H:%M:%S)] severity centres (ventilated mean anchor): ${SEV_CENTER}"
  fi
fi

# ---- 5 fits, arm by arm
fit_arm ventilated imv "$MARKERS"
for BAND in $SF_BANDS; do fit_arm "ventilated_sf${BAND/,/to}" imv "$MARKERS" PBWPFVC_JM_SF_BAND="$BAND"; done
[[ -n "$SEV_CENTER" ]] && fit_arm nosupport nosupport "$CONTROL_MARKERS" PBWPFVC_JM_SEV_CENTER="$SEV_CENTER"
# the hypoxemic control: the same control, index-day SF <= 315, without the creatinine fits
if [[ -n "$SEV_CENTER" && -n "$HYPOXEMIC_CONTROL_MARKERS" ]]; then
  CREATININE_ALL_ARMS=$CREATININE; CREATININE=0
  fit_arm nosupport_hypoxemic nosupport "$HYPOXEMIC_CONTROL_MARKERS" PBWPFVC_JM_SEV_CENTER="$SEV_CENTER" PBWPFVC_JM_SF_BAND="0,315"
  CREATININE=$CREATININE_ALL_ARMS
fi

# ---- 6 comparison table and the figure
FIG_MARKERS="platelets,bilirubin$([[ $CREATININE == 1 ]] && echo ",creatinine"),pressor_dose,osi,sf"
run_step comparison with_cohort imv Rscript code/27_control_comparison.R
run_step figure with_cohort imv env PBWPFVC_JM_WITH_RRT=1 PBWPFVC_FIG_MARKERS="$FIG_MARKERS" \
  Rscript code/24_biotrauma_figures.R

if [[ $DRY == 0 ]]; then
  echo
  echo "figure  -> $ROOT/final/injury/biotrauma_fig_main_pfvc_7d_${BASE_SITE}.pdf"
  echo "tables  -> $ROOT/final/injury/ and $ROOT/final/controls/"
  if (( ${#FAILED[@]} )); then echo "FAILED steps (${#FAILED[@]}): ${FAILED[*]}"; exit 1; else echo "all steps ok"; fi
fi
