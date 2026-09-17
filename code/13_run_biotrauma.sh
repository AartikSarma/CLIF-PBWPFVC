#!/usr/bin/env bash
# =============================================================================
# Script 13 (runner): the biotrauma analysis for one site
# =============================================================================
# One command per site. Produces, for creatinine, platelets, any vasopressor
# and the SF ratio over the first 48 hours of ventilation:
#   * the joint model of each marker's trajectory with death and extubation
#     (13_biotrauma_fit.R, pfvc form: log PFVC as a level and a divergence over
#     time beside the clinician's dose), adjusted and unadjusted, its report
#     (jm_level_contrast_*, jm_estimates_*, jm_association_hr_*) and figures
#     (biotrauma_fig_trajectory_* is the headline: the predicted marker
#     difference from the median-PFVC patient at +/- 1 SD over the window)
#   * the fixed-horizon comparator per marker (13_injury_at_horizon.R) and the
#     longitudinal model alone (13_quick_lme.R), with the channel decomposition
#   * one summary table (13_overnight_summary.R)
#
# Every stage is cached: a rerun skips panels, comparators and quick LMEs whose
# outputs exist, and the fit script reuses every finished fit from its result
# file, so a crash costs only the fit that was running. FRESH=1 redoes all.
#
# Memory: fits run PAR at a time (default 1), each with CHAINS processes, and the
# stored draws are thinned; a 7,000-patient site needs roughly 15-25 GB per fit.
# Time at 10,000 iterations: about 40-80 min per lab/SF fit, longer for the
# vasopressor part; eight fits two at a time is roughly 5-7 hours.
#
# Usage (from anywhere; the script moves to the repo root):
#   caffeinate -i nohup bash code/13_run_biotrauma.sh > biotrauma.out 2>&1 &
#   bash code/13_run_biotrauma.sh --dry-run
# Knobs (environment): HORIZONS ("48"; add 72 24 for the sensitivities), FORMS
#   ("pfvc"; vtpfvc = the same contrast told as VT/PFVC at a given VT/PBW;
#   channels, disc_level available), MARKERS_INJ, MARKERS_JM, ITER,
#   BURNIN, CHAINS, THIN, PAR, HEARTBEAT, FRESH (1: ignore every cache),
#   SKIP_JM (1: no joint models), plus PBWPFVC_SITE_NAME / PBWPFVC_TABLES_PATH.
# Logs: output/{site}_output/logs/biotrauma_{stamp}/{stage}.log and status.tsv.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

MARKERS_INJ=${MARKERS_INJ:-creatinine,platelets,ne_equiv,sf}     # comparator / quick-LME names
MARKERS_JM=${MARKERS_JM:-creatinine,platelets,any_pressor,sf}     # fit-script names (any_pressor = the hurdle's binary part)
HORIZONS=${HORIZONS:-48}
FORMS=${FORMS:-pfvc}
ITER=${ITER:-10000}; BURNIN=${BURNIN:-2000}; CHAINS=${CHAINS:-3}; THIN=${THIN:-5}
PAR=${PAR:-1}                  # fits at a time; raise to 2 once one fit has been watched to fit in memory
HEARTBEAT=${HEARTBEAT:-300}
FRESH=${FRESH:-0}              # 1: ignore every cache and redo everything
SKIP_JM=${SKIP_JM:-0}          # 1: panels, comparators, quick LMEs and summary only
COHORT=${COHORT:-imv}          # niv (HFNC/NIV first, the middle arm) or nosupport (room air / cannula, the control); under {site}_{cohort}, no joint models
DRY=0; [[ "${1:-}" == "--dry-run" ]] && DRY=1

# site name from config.json by sed: Rscript's stdout carries renv's start-up notices
SITE=${PBWPFVC_SITE_NAME:-$(sed -n 's/.*"site_name" *: *"\([^"]*\)".*/\1/p' config/config.json | head -n 1)}
[[ -n "$SITE" ]] || { echo "could not read site_name from config/config.json"; exit 1; }
case $COHORT in
  imv) ;;
  niv|nosupport) SITE="${SITE%_niv}"; SITE="${SITE%_nosupport}_${COHORT}"; export PBWPFVC_SITE_NAME=$SITE PBWPFVC_COHORT=$COHORT; SKIP_JM=1 ;;
  *) echo "COHORT must be imv, niv or nosupport"; exit 1 ;;
esac
FINAL="output/${SITE}_output/final"
INTER="output/${SITE}_output/intermediate"
STAMP=$(date +%Y%m%d_%H%M%S)
LOGDIR="output/${SITE}_output/logs/biotrauma_${STAMP}"
STATUS="$LOGDIR/status.tsv"
set_rc() { eval "RC_${1//[^A-Za-z0-9]/_}=$2"; }
rc_of()  { eval "echo \${RC_${1//[^A-Za-z0-9]/_}:-1}"; }

echo "site $SITE (cohort $COHORT); markers $MARKERS_INJ / $MARKERS_JM; horizons $HORIZONS; forms $FORMS; chains $ITER/$BURNIN x $CHAINS, thin $THIN, $PAR fits at a time; fresh $FRESH"
# preflight: the site's script-03 outputs must exist. The control cohort builds
# its own (scripts 01-03 under {site}_niv); the analytic cohort must be built first.
need="$INTER/ne_equiv_admin.parquet"
if [[ $DRY == 0 && ! -f "$need" ]]; then
  if [[ $COHORT != imv ]]; then
    echo "[$(date +%T)] building the $COHORT cohort for $SITE: scripts 01-03"
    mkdir -p "output/${SITE}_output"
    for s in 01_cohort_identification 02_quality_checks 03_variable_derivation; do
      Rscript "code/$s.R" > "output/${SITE}_output/${s}.log" 2>&1 || { echo "ABORT: $s failed (output/${SITE}_output/${s}.log)"; exit 1; }
      echo "[$(date +%T)] $s done"
    done
  else
    echo "ABORT: $need is missing. config/config.json names site '$SITE'; run scripts 01-03 for it, or point the config at the site you meant."
    exit 1
  fi
fi
if [[ $DRY == 0 ]]; then
  mkdir -p "$LOGDIR"
  printf 'stage\tstart\tend\texit\n' > "$STATUS"
  stale=$(ls "$FINAL"/quick_lme_*_[0-9]*h_"$SITE".csv 2>/dev/null || true)
  [[ -n "$stale" ]] && echo "WARNING: stale horizon-tagged quick_lme files present (older script); excluded by the summary and pooling:" && echo "$stale"
  echo "logs -> $LOGDIR"
fi

# run_stage <name> [VAR=value ...] -- <command ...>
run_stage() {
  local name=$1; shift
  local envs=()
  while [[ "$1" != "--" ]]; do envs+=("$1"); shift; done; shift
  if [[ $DRY == 1 ]]; then echo "[dry] $name: ${envs[@]+"${envs[@]}"} $*"; set_rc "$name" 0; return 0; fi
  local t0; t0=$(date +%FT%T)
  echo "[$(date +%H:%M:%S)] START $name"
  env ${envs[@]+"${envs[@]}"} "$@" > "$LOGDIR/$name.log" 2>&1
  local rc=$?
  set_rc "$name" $rc
  printf '%s\t%s\t%s\t%s\n' "$name" "$t0" "$(date +%FT%T)" "$rc" >> "$STATUS"
  echo "[$(date +%H:%M:%S)] END   $name (exit $rc)"
  return 0
}
# skip_if <name> <file>: a cached stage (its output exists and FRESH is off)
skip_if() {
  local name=$1 f=$2
  if [[ $FRESH == 0 && $DRY == 0 && -f "$f" ]]; then
    echo "[$(date +%H:%M:%S)] CACHED $name ($(basename "$f") exists)"
    set_rc "$name" 0; printf '%s\t%s\t%s\t%s\n' "$name" cached cached 0 >> "$STATUS"; return 0
  fi
  return 1
}
# number of usable fits (converged or rhat_fail) in a manifest; 0 when absent
usable_fits() {
  local f=$1
  [[ -f "$f" ]] || { echo 0; return; }
  Rscript -e "m <- read.csv('$f'); cat(sum(m\$status %in% c('converged', 'rhat_fail')), '\n')" 2>/dev/null | tail -n 1 | tr -dc '0-9'
}
tag_of() { local form=$1 h=$2; echo "${form}_${h}h_${SITE}"; }

# ---- 1. panels (the 72-hour one also serves the quick LMEs)
PANELS="$HORIZONS"; [[ " $HORIZONS " == *" 72 "* ]] || PANELS="$HORIZONS 72"
for H in $PANELS; do
  skip_if "panel_${H}h" "$INTER/jm_long_${H}h.parquet" || \
    run_stage "panel_${H}h" PBWPFVC_JM_HORIZON_H=$H -- Rscript code/13_biotrauma_panel.R
done

# ---- 2. comparator (all horizons inside the script)
for m in ${MARKERS_INJ//,/ }; do
  skip_if "injury_${m}" "$FINAL/injury_negctrl_${m}_${SITE}.csv" || \
    run_stage "injury_${m}" PBWPFVC_INJ_MARKER=$m -- Rscript code/13_injury_at_horizon.R
done

# ---- 3. quick LME (72-hour panel; both exposures inside the script)
if [[ $(rc_of panel_72h) -eq 0 ]]; then
  for m in ${MARKERS_INJ//,/ }; do
    skip_if "quick_${m}" "$FINAL/quick_sf_channels_${m}_${SITE}.csv" || \
      run_stage "quick_${m}" PBWPFVC_INJ_MARKER=$m PBWPFVC_JM_HORIZON_H=72 -- Rscript code/13_quick_lme.R
  done
else
  echo "quick LME skipped: the 72-hour panel failed"
fi

# ---- 4. joint models: per horizon, per form; finished fits are reused by the fit
#         script itself; report + figures when the fit left usable models
[[ $SKIP_JM == 1 ]] && echo "joint models skipped (SKIP_JM=1)"
for H in $HORIZONS; do
  [[ $SKIP_JM == 1 ]] && break
  if [[ $(rc_of panel_${H}h) -ne 0 ]]; then echo "fits at ${H}h skipped: panel failed"; continue; fi
  for FORM in $FORMS; do
    run_stage "fit_${FORM}_${H}h" PBWPFVC_JM_HORIZON_H=$H PBWPFVC_JM_MODIFIER=$FORM PBWPFVC_JM_MARKERS=$MARKERS_JM \
      PBWPFVC_JM_MODELS=main PBWPFVC_JM_ITER=$ITER PBWPFVC_JM_BURNIN=$BURNIN PBWPFVC_JM_CHAINS=$CHAINS \
      PBWPFVC_JM_THIN=$THIN PBWPFVC_JM_PAR=$PAR PBWPFVC_JM_FRESH=$FRESH PBWPFVC_JM_HEARTBEAT=$HEARTBEAT \
      -- Rscript code/13_biotrauma_fit.R
    n_ok=$(usable_fits "$FINAL/jm_manifest_$(tag_of $FORM $H).csv"); [[ $DRY == 1 ]] && n_ok=1
    if [[ $(rc_of fit_${FORM}_${H}h) -eq 0 && ${n_ok:-0} -gt 0 ]]; then
      run_stage "report_${FORM}_${H}h"  PBWPFVC_JM_HORIZON_H=$H PBWPFVC_JM_MODIFIER=$FORM -- Rscript code/13_biotrauma_report.R
      run_stage "figures_${FORM}_${H}h" PBWPFVC_JM_HORIZON_H=$H PBWPFVC_JM_MODIFIER=$FORM -- Rscript code/13_biotrauma_figures.R
    else
      echo "report/figures for ${FORM} ${H}h skipped: fit exit $(rc_of fit_${FORM}_${H}h), usable fits ${n_ok:-0}"
    fi
  done
done

# ---- 5. summary
run_stage "summary" -- Rscript code/13_overnight_summary.R
[[ $DRY == 0 ]] && { echo; echo "---- status.tsv"; cat "$STATUS"; echo "summary -> $FINAL/overnight_summary_${SITE}.csv"; echo "headline figure -> $FINAL/biotrauma_fig_trajectory_pfvc_48h_${SITE}.pdf"; }
exit 0
