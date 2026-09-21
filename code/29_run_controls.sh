#!/bin/bash
# =============================================================================
# 29_run_controls: the control cohorts, built together, and every arm of the
# divergence-by-lung-size comparison
# =============================================================================
# Everything lands in the site's one output folder: the controls' aggregates in
# output/{site}_output/final/controls/, the comparison in final/injury/. Three stages,
# because the severity floor is chosen from the ventilated cohort's anchor
# distribution and that has to be read before the matched fits run.
#
#   bash code/29_run_controls.sh build
#       scripts 01-03 for every control cohort (CONTROL_COHORTS, default
#       "nosupport"). The ventilated cohort must already be built (00_run_pipeline.R).
#       The noninvasive cohort (niv) is not a control: NIPPV delivers large, unlimited
#       positive-pressure volumes. Build it only on request, CONTROL_COHORTS="nosupport niv".
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
# Knobs (environment): CONTROL_COHORTS, MARKERS, SEV_MIN, SF_BANDS (set it empty to skip
# the ventilated SF strata), ITER, BURNIN, CHAINS, THIN, HORIZON.
# A failed arm is reported and the run continues: the arms are independent.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

STAGE=${1:-}
MARKERS=${MARKERS:-platelets,bilirubin}
SF_BANDS=${SF_BANDS-"235,315 115,235 0,115"}   # SF_BANDS= (empty) skips the ventilated SF strata
SEV_MIN=${SEV_MIN:-}
CONTROL_COHORTS=${CONTROL_COHORTS:-nosupport}
ITER=${ITER:-2000}; BURNIN=${BURNIN:-500}; CHAINS=${CHAINS:-3}; THIN=${THIN:-5}
HORIZON=${HORIZON:-7}

BASE_SITE=${PBWPFVC_SITE_NAME:-$(sed -n 's/.*"site_name" *: *"\([^"]*\)".*/\1/p' config/config.json | head -n 1)}
[ -n "$BASE_SITE" ] || { echo "could not read site_name from config/config.json"; exit 1; }
# an older shell may still export PBWPFVC_SITE_NAME={site}_{cohort}; utils/config.R strips it, so must this
BASE_SITE="${BASE_SITE%_niv}"; BASE_SITE="${BASE_SITE%_nosupport}"
LOG_DIR="output/${BASE_SITE}_output/logs/controls_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"
echo "site ${BASE_SITE}; markers ${MARKERS}; logs in ${LOG_DIR}"

export PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=$HORIZON PBWPFVC_JM_MODIFIER=pfvc
export PBWPFVC_JM_MARKERS=$MARKERS PBWPFVC_JM_MODELS=main
export PBWPFVC_JM_ITER=$ITER PBWPFVC_JM_BURNIN=$BURNIN PBWPFVC_JM_CHAINS=$CHAINS PBWPFVC_JM_THIN=$THIN

use_cohort () {   # imv | nosupport | niv; utils/config.R routes a control into the site's own folder
  if [ "$1" = imv ]; then unset PBWPFVC_COHORT; else export PBWPFVC_COHORT=$1; fi
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
  build)
    for COHORT in $CONTROL_COHORTS; do
      use_cohort $COHORT
      for SCRIPT in 01_cohort_identification 02_quality_checks 03_variable_derivation; do
        run_step "${COHORT}_${SCRIPT}" Rscript code/${SCRIPT}.R
      done
    done
    use_cohort imv
    ;;
  anchors)
    for COHORT in imv nosupport; do
      use_cohort $COHORT
      run_step "${COHORT}_panel"   Rscript code/21_biotrauma_panel.R
      run_step "${COHORT}_anchors" env PBWPFVC_JM_ANCHOR_ONLY=1 Rscript code/22_biotrauma_fit.R
    done
    echo "anchor distributions:"
    echo "  output/${BASE_SITE}_output/final/injury/jm_severity_anchor_${HORIZON}d_${BASE_SITE}.csv"
    echo "  output/${BASE_SITE}_output/final/controls/jm_severity_anchor_${HORIZON}d_${BASE_SITE}_nosupport.csv"
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
  *) echo "usage: bash code/29_run_controls.sh build | anchors | fits   (fits needs SEV_MIN=...)"; exit 1 ;;
esac
