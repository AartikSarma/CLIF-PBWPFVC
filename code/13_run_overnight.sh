#!/usr/bin/env bash
# =============================================================================
# Script 13 (runner): the overnight biotrauma run, one site
# =============================================================================
# Runs, in order, each as a separate Rscript subprocess with its own log:
#   1. the period panels for every horizon (13_biotrauma_panel.R)
#   2. the fixed-horizon comparator per marker (13_injury_at_horizon.R): the
#      survivor-only contrast, the channel decomposition and its supports
#   3. the quick LME per marker (13_quick_lme.R, on the 72-hour panel)
#   4. the joint models per horizon and form (13_biotrauma_fit.R), each followed
#      by the report and the figures when the fit produced usable models
#   5. the summary (13_overnight_summary.R)
# Independent stages continue after a failure; status.tsv records every stage.
# The fit script exits 0 even when every fit fails (per-fit errors become
# manifest rows), so the report/figures stages are gated on the manifest.
#
# Usage (from anywhere; the script moves to the repo root):
#   caffeinate -i nohup bash code/13_run_overnight.sh > overnight.out 2>&1 &
#   bash code/13_run_overnight.sh --dry-run          # print the stages only
# Knobs (environment): MARKERS_INJ, MARKERS_JM, HORIZONS, FORMS, ITER, BURNIN,
#   CHAINS, HEARTBEAT (fit progress interval, s), RESUME (1: reuse fit bundles
#   already on disk, fit only what is missing), plus PBWPFVC_SITE_NAME /
#   PBWPFVC_TABLES_PATH as in utils/config.R.
# Logs: output/{site}_output/logs/overnight_{stamp}/{stage}.log and status.tsv;
# the existing jm_*/injury_*/quick_* tables are copied there first, because the
# fits run with PBWPFVC_JM_FRESH=1 and replace their tables.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

MARKERS_INJ=${MARKERS_INJ:-creatinine,platelets,ne_equiv,sf}     # comparator / quick-LME names
MARKERS_JM=${MARKERS_JM:-creatinine,platelets,any_pressor,sf}     # fit-script names (any_pressor = the hurdle's binary part)
HORIZONS=${HORIZONS:-48 72 24}                                    # fit order: the primary window first
FORMS=${FORMS:-pfvc channels disc_level}
ITER=${ITER:-25000}; BURNIN=${BURNIN:-5000}; CHAINS=${CHAINS:-3}
HEARTBEAT=${HEARTBEAT:-300}
RESUME=${RESUME:-0}            # 1: rebuild tables from fit bundles already on disk, fit only what is missing
DRY=0; [[ "${1:-}" == "--dry-run" ]] && DRY=1

# site name from config.json by sed: Rscript's stdout carries renv's start-up notices
SITE=${PBWPFVC_SITE_NAME:-$(sed -n 's/.*"site_name" *: *"\([^"]*\)".*/\1/p' config/config.json | head -n 1)}
[[ -n "$SITE" ]] || { echo "could not read site_name from config/config.json"; exit 1; }
FINAL="output/${SITE}_output/final"
STAMP=$(date +%Y%m%d_%H%M)
LOGDIR="output/${SITE}_output/logs/overnight_${STAMP}"
STATUS="$LOGDIR/status.tsv"
# per-stage exit codes in plain variables (macOS bash 3.2 has no associative arrays)
set_rc() { eval "RC_${1//[^A-Za-z0-9]/_}=$2"; }
rc_of()  { eval "echo \${RC_${1//[^A-Za-z0-9]/_}:-1}"; }

echo "site $SITE; markers $MARKERS_INJ / $MARKERS_JM; horizons $HORIZONS; forms $FORMS; chains $ITER/$BURNIN x $CHAINS"
if [[ $DRY == 0 ]]; then
  mkdir -p "$LOGDIR/previous_tables"
  printf 'stage\tstart\tend\texit\n' > "$STATUS"
  # keep what the fresh fits will replace
  find "$FINAL" -maxdepth 1 \( -name 'jm_*' -o -name 'injury_*' -o -name 'quick_*' -o -name 'biotrauma_fig_*' \) \
    -exec cp {} "$LOGDIR/previous_tables/" \; 2>/dev/null
  stale=$(ls "$FINAL"/quick_lme_*_[0-9]*h_"$SITE".csv 2>/dev/null || true)
  [[ -n "$stale" ]] && echo "WARNING: stale horizon-tagged quick_lme files present (older script); they are excluded by the summary and pooling:" && echo "$stale"
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
# number of usable fits (converged or rhat_fail) in a manifest; 0 when absent
usable_fits() {
  local f=$1
  [[ -f "$f" ]] || { echo 0; return; }
  Rscript -e "m <- read.csv('$f'); cat(sum(m\$status %in% c('converged', 'rhat_fail')), '\n')" 2>/dev/null | tail -n 1 | tr -dc '0-9'
}
tag_of() { local form=$1 h=$2; echo "${form}_${h}h_${SITE}"; }

# ---- 1. panels
for H in $HORIZONS; do
  run_stage "panel_${H}h" PBWPFVC_JM_HORIZON_H=$H -- Rscript code/13_biotrauma_panel.R
done
[[ " $HORIZONS " == *" 72 "* ]] || run_stage "panel_72h" PBWPFVC_JM_HORIZON_H=72 -- Rscript code/13_biotrauma_panel.R

# ---- 2. comparator (all horizons inside the script)
for m in ${MARKERS_INJ//,/ }; do
  run_stage "injury_${m}" PBWPFVC_INJ_MARKER=$m -- Rscript code/13_injury_at_horizon.R
done

# ---- 3. quick LME (72-hour panel; both exposures inside the script)
if [[ $(rc_of panel_72h) -eq 0 ]]; then
  for m in ${MARKERS_INJ//,/ }; do
    run_stage "quick_${m}" PBWPFVC_INJ_MARKER=$m PBWPFVC_JM_HORIZON_H=72 -- Rscript code/13_quick_lme.R
  done
else
  echo "quick LME skipped: the 72-hour panel failed"
fi

# ---- 4. joint models: per horizon, per form; report + figures when the fit left usable models
for H in $HORIZONS; do
  if [[ $(rc_of panel_${H}h) -ne 0 ]]; then echo "fits at ${H}h skipped: panel failed"; continue; fi
  for FORM in $FORMS; do
    run_stage "fit_${FORM}_${H}h" PBWPFVC_JM_HORIZON_H=$H PBWPFVC_JM_MODIFIER=$FORM PBWPFVC_JM_MARKERS=$MARKERS_JM \
      PBWPFVC_JM_MODELS=main PBWPFVC_JM_ITER=$ITER PBWPFVC_JM_BURNIN=$BURNIN PBWPFVC_JM_CHAINS=$CHAINS \
      PBWPFVC_JM_FRESH=1 PBWPFVC_JM_HEARTBEAT=$HEARTBEAT PBWPFVC_JM_RESUME=$RESUME -- Rscript code/13_biotrauma_fit.R
    n_ok=$(usable_fits "$FINAL/jm_manifest_$(tag_of $FORM $H).csv"); [[ $DRY == 1 ]] && n_ok=1
    if [[ $(rc_of fit_${FORM}_${H}h) -eq 0 && $n_ok -gt 0 ]]; then
      run_stage "report_${FORM}_${H}h"  PBWPFVC_JM_HORIZON_H=$H PBWPFVC_JM_MODIFIER=$FORM -- Rscript code/13_biotrauma_report.R
      run_stage "figures_${FORM}_${H}h" PBWPFVC_JM_HORIZON_H=$H PBWPFVC_JM_MODIFIER=$FORM -- Rscript code/13_biotrauma_figures.R
    else
      echo "report/figures for ${FORM} ${H}h skipped: fit exit $(rc_of fit_${FORM}_${H}h), usable fits $n_ok"
    fi
  done
done

# ---- 5. summary
run_stage "summary" -- Rscript code/13_overnight_summary.R
[[ $DRY == 0 ]] && { echo; echo "---- status.tsv"; cat "$STATUS"; echo "summary -> $FINAL/overnight_summary_${SITE}.csv"; }
exit 0
