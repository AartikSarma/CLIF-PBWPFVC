#!/bin/bash
# =============================================================================
# One-time migration: sort a flat final/ folder into its block subfolders
# =============================================================================
# Before 2026-09-19 every script wrote to output/{site}_output/final/. Each now
# writes to final/<block>/ (utils/config.R, final_dir_for()):
#   cross_sectional/  figures 1-3     injury/   figure 4     causal/  figure 5
#   supplement/       what code/supplement/ scripts write    controls/  as before
# This moves the files already on disk, by file-name prefix, so nothing is rerun.
# File names do not change, with one exception: overnight_summary_* becomes
# biotrauma_summary_*, the name 28_biotrauma_summary.R now writes.
#
# Usage:  bash code/tools/migrate_final_subfolders.sh SITE            (dry run)
#         bash code/tools/migrate_final_subfolders.sh SITE --apply
# Only files at the top of final/ are touched; subfolders (controls/, the TTE's
# tagged sensitivity folders) are left alone. Nothing is overwritten. A file whose
# prefix no rule knows is listed and left where it is: no current script writes it,
# so keep it or delete it by hand.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
SITE=${1:-}; APPLY=${2:-}
[ -n "$SITE" ] || { echo "usage: bash code/tools/migrate_final_subfolders.sh SITE [--apply]"; exit 1; }
FINAL="output/${SITE}_output/final"
[ -d "$FINAL" ] || { echo "no $FINAL"; exit 1; }

block_of () {   # file name -> block; the first matching rule wins, so the narrow prefixes come first
  case "$1" in
    norm_bracket_*|harm_*|pos_*|dose_*|cbias_export_*|ivpolicy_*|sevenday_*|oi_*|dp_by_*|dp_vs_*|vtpbw_band_scan*|sawtooth_*) echo supplement ;;
    jm_*|injury_*|quick_*|biotrauma_fig_*|biotrauma_summary_*|overnight_summary_*) echo injury ;;
    tte_*|vtpbw_titration_*|mppbw_additive_*|iv_preference*) echo causal ;;
    regression_*|norm_*|bias_*|table*|negative_control*|evalues_*|cox_*|km_*|aic_*|evidence_ratio_*|size_*|distribution_*|consort_*|attrition_log_*|dist_*|lab_summary_*|vital_summary_*) echo cross_sectional ;;
    *) echo "" ;;
  esac
}

n_moved=0; n_kept=0; unknown=()
for path in "$FINAL"/*; do
  [ -f "$path" ] || continue
  name=$(basename "$path"); block=$(block_of "$name")
  if [ -z "$block" ]; then unknown+=("$name"); continue; fi
  target="$FINAL/$block/${name/#overnight_summary_/biotrauma_summary_}"
  if [ -e "$target" ]; then n_kept=$((n_kept + 1)); echo "  exists, left in place: $target"; continue; fi
  n_moved=$((n_moved + 1))
  if [ "$APPLY" = "--apply" ]; then mkdir -p "$FINAL/$block" && mv "$path" "$target"; fi
done
echo "$FINAL: $n_moved files to move, $n_kept already at their destination"
for block in cross_sectional injury causal supplement; do
  [ -d "$FINAL/$block" ] && echo "  $block/: $(find "$FINAL/$block" -maxdepth 1 -type f | wc -l | tr -d ' ') files"
done
if [ ${#unknown[@]} -gt 0 ]; then
  echo "no rule for these ${#unknown[@]} files (no current script writes them); left in $FINAL:"
  printf '  %s\n' "${unknown[@]}"
fi
[ "$APPLY" = "--apply" ] || echo "dry run: nothing moved. Re-run with --apply."
