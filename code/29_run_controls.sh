#!/bin/bash
# =============================================================================
# the biotrauma suite (2x) (controls runner): every arm of the divergence-by-lung-size comparison
# =============================================================================
# Two stages, because the severity floor is chosen from the ventilated cohort's
# anchor distribution and that has to be read first.
#
#   bash code/29_run_controls.sh anchors
#       rebuilds the 7-day panel of both cohorts (they must carry the SOFA
#       components) and writes final/jm_severity_anchor_* for each. No fits.
#
#   SEV_MIN="platelets=2,bilirubin=1" bash code/29_run_controls.sh fits
#       fits and reports, in order: ventilated; ventilated in each baseline SF
#       class (235-315, 115-235, <=115); ventilated above the floor; no support;
#       no support above the floor. Then 27_control_comparison.R. Fits already on
#       disk with the same chain settings are reused, so the two unrestricted
#       arms cost nothing if they have been run.
#
# The no-support cohort needs its own scripts 01-03 outputs ({site}_nosupport).
# Knobs (environment): MARKERS, SEV_MIN, SF_BANDS, ITER, BURNIN, CHAINS, THIN, HORIZON.
# A failed arm is reported and the run continues: the arms are independent.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

STAGE=${1:-}
MARKERS=${MARKERS:-platelets,bilirubin}
SF_BANDS=${SF_BANDS:-"235,315 115,235 0,115"}
SEV_MIN=${SEV_MIN:-}
ITER=${ITER:-2000}; BURNIN=${BURNIN:-500}; CHAINS=${CHAINS:-3}; THIN=${THIN:-5}
HORIZON=${HORIZON:-7}

BASE_SITE=$(Rscript -e 'cat(jsonlite::fromJSON("config/config.json")$site_name)' 2>/dev/null)
[ -n "$BASE_SITE" ] || { echo "could not read site_name from config/config.json"; exit 1; }
LOG_DIR="output/${BASE_SITE}_output/logs/controls_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"
echo "site ${BASE_SITE}; markers ${MARKERS}; logs in ${LOG_DIR}"

export PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=$HORIZON PBWPFVC_JM_MODIFIER=pfvc
export PBWPFVC_JM_MARKERS=$MARKERS PBWPFVC_JM_MODELS=main
export PBWPFVC_JM_ITER=$ITER PBWPFVC_JM_BURNIN=$BURNIN PBWPFVC_JM_CHAINS=$CHAINS PBWPFVC_JM_THIN=$THIN

use_cohort () {   # imv | nosupport
  if [ "$1" = imv ]; then unset PBWPFVC_COHORT PBWPFVC_SITE_NAME
  else export PBWPFVC_COHORT=$1 PBWPFVC_SITE_NAME="${BASE_SITE}_$1"; fi
}
run_step () {     # name, then the command
  local step_name=$1; shift
  echo "[$(date +%H:%M:%S)] ${step_name}"
  if "$@" > "${LOG_DIR}/${step_name}.log" 2>&1; then echo "    ok"
  else echo "    FAILED, see ${LOG_DIR}/${step_name}.log"; tail -5 "${LOG_DIR}/${step_name}.log" | sed 's/^/    | /'; fi
}
run_arm () {      # name, then extra environment assignments for this arm
  local arm_name=$1; shift
  run_step "${arm_name}_fit"    env "$@" Rscript code/22_biotrauma_fit.R
  run_step "${arm_name}_report" env "$@" Rscript code/23_biotrauma_report.R
}

case "$STAGE" in
  anchors)
    for COHORT in imv nosupport; do
      use_cohort $COHORT
      run_step "${COHORT}_panel"   Rscript code/21_biotrauma_panel.R
      run_step "${COHORT}_anchors" env PBWPFVC_JM_ANCHOR_ONLY=1 Rscript code/22_biotrauma_fit.R
    done
    echo "anchor distributions:"
    echo "  output/${BASE_SITE}_output/final/jm_severity_anchor_${HORIZON}d_${BASE_SITE}.csv"
    echo "  output/${BASE_SITE}_nosupport_output/final/jm_severity_anchor_${HORIZON}d_${BASE_SITE}_nosupport.csv"
    ;;
  fits)
    [ -n "$SEV_MIN" ] || { echo "set SEV_MIN, e.g. SEV_MIN=\"platelets=2,bilirubin=1\" (run the anchors stage first)"; exit 1; }
    use_cohort imv
    run_arm ventilated PBWPFVC_UNRESTRICTED=1
    for BAND in $SF_BANDS; do run_arm "ventilated_sf${BAND/,/to}" PBWPFVC_JM_SF_BAND=$BAND; done
    run_arm ventilated_matched PBWPFVC_JM_SEV_MIN=$SEV_MIN
    use_cohort nosupport
    run_arm nosupport PBWPFVC_UNRESTRICTED=1
    run_arm nosupport_matched PBWPFVC_JM_SEV_MIN=$SEV_MIN
    use_cohort imv
    run_step comparison Rscript code/27_control_comparison.R
    grep -A40 "divergence per day" "${LOG_DIR}/comparison.log" | grep -v "^Warning"
    ;;
  *) echo "usage: bash code/29_run_controls.sh anchors | SEV_MIN=... bash code/29_run_controls.sh fits"; exit 1 ;;
esac
