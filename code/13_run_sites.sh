#!/usr/bin/env bash
# =============================================================================
# Script 13 (sites): the biotrauma analysis for several sites, one after another
# =============================================================================
# For each site given as NAME=TABLES_PATH[=FILE_TYPE]: sets the config overrides
# utils/config.R reads (PBWPFVC_SITE_NAME, _TABLES_PATH, _FILE_TYPE; the config
# file is not edited), builds the cohort with scripts 01-03 when the site's
# script-03 outputs are missing, then runs code/13_run_biotrauma.sh. Sites run
# sequentially, so memory is one site's. Every knob of 13_run_biotrauma.sh
# (HORIZONS, FORMS, ITER, PAR, FRESH, ...) passes through the environment.
#
# Usage:
#   caffeinate -i nohup bash code/13_run_sites.sh \
#       UCSF=/path/to/ucsf_clif MIMIC=/path/to/mimic_clif > sites.out 2>&1 &
# CONTROL (default "nosupport niv") names the control arms run after each site:
# nosupport = room air / cannula only (the negative control), niv = HFNC/NIV first
# (the middle arm of the strain gradient); each built and analysed under
# output/{site}_{arm}_output/. CONTROL="" skips them.
#   bash code/13_run_sites.sh --dry-run UCSF=/path/to/ucsf_clif
# Each site's own log is output/{site}_output/biotrauma.out; the per-stage logs
# are under output/{site}_output/logs/biotrauma_{stamp}/.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
CONTROL=${CONTROL:-nosupport niv}   # control arms run after each site (COHORT values; empty string = none)
DRY=0; [[ "${1:-}" == "--dry-run" ]] && { DRY=1; shift; }
[[ $# -ge 1 ]] || { echo "usage: bash code/13_run_sites.sh [--dry-run] NAME=TABLES_PATH[=FILE_TYPE] ..."; exit 1; }

for spec in "$@"; do
  name=${spec%%=*}; rest=${spec#*=}; path=${rest%%=*}; ftype=parquet
  [[ "$rest" == *=* ]] && ftype=${rest#*=}
  path=${path/#\~/$HOME}
  echo "================================================================"
  echo "[$(date +%FT%T)] site $name  tables $path  ($ftype)"
  if [[ ! -d "$path" ]]; then echo "ABORT: tables path does not exist: $path"; continue; fi
  export PBWPFVC_SITE_NAME=$name PBWPFVC_TABLES_PATH=$path PBWPFVC_FILE_TYPE=$ftype
  out="output/${name}_output"; mkdir -p "$out"
  if [[ ! -f "$out/intermediate/ne_equiv_admin.parquet" ]]; then
    echo "[$(date +%T)] cohort not built for $name: running scripts 01-03"
    for s in 01_cohort_identification 02_quality_checks 03_variable_derivation; do
      if [[ $DRY == 1 ]]; then echo "[dry] Rscript code/$s.R"; continue; fi
      Rscript "code/$s.R" > "$out/${s}.log" 2>&1; rc=$?
      echo "[$(date +%T)] $s exit $rc (log $out/${s}.log)"
      if [[ $rc -ne 0 ]]; then echo "ABORT $name: $s failed"; continue 2; fi
    done
  else
    echo "[$(date +%T)] cohort present for $name; scripts 01-03 skipped"
  fi
  if [[ $DRY == 1 ]]; then bash code/13_run_biotrauma.sh --dry-run | head -3; continue; fi
  bash code/13_run_biotrauma.sh > "$out/biotrauma.out" 2>&1
  echo "[$(date +%FT%T)] $name done; $(grep -c 'exit [1-9]' "$out/biotrauma.out") failed stages; log $out/biotrauma.out"
  tail -n 4 "$out/biotrauma.out"
  for arm in $CONTROL; do
    mkdir -p "output/${name}_${arm}_output"
    COHORT=$arm bash code/13_run_biotrauma.sh > "output/${name}_${arm}_output/biotrauma.out" 2>&1
    echo "[$(date +%FT%T)] ${name}_${arm} done; $(grep -c 'exit [1-9]' "output/${name}_${arm}_output/biotrauma.out") failed stages"
    tail -n 3 "output/${name}_${arm}_output/biotrauma.out"
  done
done
unset PBWPFVC_SITE_NAME PBWPFVC_TABLES_PATH PBWPFVC_FILE_TYPE
echo "[$(date +%FT%T)] all sites done"
