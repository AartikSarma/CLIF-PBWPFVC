#!/bin/bash
# =============================================================================
# One-time migration: fold the old sibling control folders into the site's folder
# =============================================================================
# Before 2026-09-18 a control cohort ran under a site name of its own and wrote to
#   output/{site}_{cohort}_output/{intermediate,final,logs}
# It now lives inside the site's one folder (utils/config.R):
#   output/{site}_output/intermediate/controls/{cohort}/
#   output/{site}_output/final/controls/
#   output/{site}_output/logs/{cohort}/
# File names do not change, so nothing needs refitting: cached fits are found again.
#
# Usage:  bash code/tools/migrate_control_folders.sh SITE            (dry run: prints the moves)
#         bash code/tools/migrate_control_folders.sh SITE --apply
# Nothing is overwritten: a file already present at the destination is left alone
# and reported. An old folder is removed only once it is empty.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
SITE=${1:-}; APPLY=${2:-}
[ -n "$SITE" ] || { echo "usage: bash code/tools/migrate_control_folders.sh SITE [--apply]"; exit 1; }
ROOT="output/${SITE}_output"
[ -d "$ROOT" ] || { echo "no $ROOT: run from a checkout that holds the site's outputs"; exit 1; }

move_tree () {   # source dir, destination dir
  local src=$1 dst=$2 n_moved=0 n_kept=0
  [ -d "$src" ] || return 0
  while IFS= read -r -d '' path; do
    rel=${path#"$src"/}
    if [ -e "$dst/$rel" ]; then n_kept=$((n_kept + 1)); echo "  exists, left in place: $dst/$rel"; continue; fi
    n_moved=$((n_moved + 1))
    if [ "$APPLY" = "--apply" ]; then mkdir -p "$(dirname "$dst/$rel")" && mv "$path" "$dst/$rel"; fi
  done < <(find "$src" -type f -print0)
  echo "  $src -> $dst: $n_moved files to move, $n_kept already there"
  if [ "$APPLY" = "--apply" ]; then find "$src" -type d -empty -delete 2>/dev/null; fi
}

for COHORT in nosupport niv; do
  OLD="output/${SITE}_${COHORT}_output"
  [ -d "$OLD" ] || { echo "$OLD: not present, nothing to do"; continue; }
  echo "$OLD:"
  move_tree "$OLD/intermediate" "$ROOT/intermediate/controls/$COHORT"
  move_tree "$OLD/final"        "$ROOT/final/controls"
  move_tree "$OLD/logs"         "$ROOT/logs/$COHORT"
  move_tree "$OLD"              "$ROOT/logs/$COHORT"      # stray top-level logs (01-03, biotrauma.out)
  if [ "$APPLY" = "--apply" ]; then rmdir "$OLD" 2>/dev/null && echo "  removed empty $OLD"; fi
done
[ "$APPLY" = "--apply" ] || echo "dry run: nothing moved. Re-run with --apply."
