#!/bin/bash
# =====================================================================
# Resolution scan, run from the WORKING FOLDER that holds the settings
# files, normally with the ./run_scan launcher that lives there:
#
#     cd ~/Dropbox/Documents/MuonColliderSimulation/Analysis
#     (edit __scan_config.txt, and if needed __smearing_config.txt,
#      __cuts_config.txt, __input_files_config.txt)
#     ./run_scan
#
# What one scan does:
#   1. Checks all four settings files and the input files first; stops
#      before running anything if something is wrong.
#   2. Creates runs/<date>_<time>_scan_<time|angle>/ and copies the four
#      settings files into it at the START.
#   3. Runs resolution_scan.py (see templates/__scan_config.txt for the
#      method), saving everything printed to run_log.txt there, plus
#      scan_results.csv, the plots and scan_<time|angle>.pdf.
#   4. If it succeeded, puts a copy of the PDF in the working folder as
#      _scan_<time|angle>_<date>_<time>.pdf (replacing the previous one of
#      the same kind). If it failed, the run folder is renamed ..._FAILED.
# The step* and _highlights folders of ./run_all are never touched.
# =====================================================================
set -u
set -o pipefail
CODE_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PYTHON="${PYTHON:-python3}"
say() { printf '%s\n' "$*"; }
stop() { printf '\nERROR: %s\n\n' "$*" >&2; exit 1; }

case "${1:-}" in -h|--help) sed -n '2,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;; esac
WORKDIR="${1:-$PWD}"
[ -d "$WORKDIR" ] || stop "working folder not found: $WORKDIR"
WORKDIR="$(cd "$WORKDIR" && pwd -P)"
SIM_DIR="${SIM_DIR:-$(dirname "$WORKDIR")}"

command -v "$PYTHON" >/dev/null 2>&1 || stop "'$PYTHON' not found"
for f in __scan_config.txt __cuts_config.txt __smearing_config.txt __input_files_config.txt; do
    if [ ! -f "$WORKDIR/$f" ]; then
        [ -f "$CODE_DIR/templates/$f" ] && cp "$CODE_DIR/templates/$f" "$WORKDIR/$f" \
            && stop "$f was missing, so a fresh copy was made from the template:
       $WORKDIR/$f
       Check its values, then run again."
        stop "$f not found in $WORKDIR"
    fi
done
PARAM="$("$PYTHON" "$CODE_DIR/resolution_scan.py" --check "$WORKDIR/__scan_config.txt")" \
    || stop "fix __scan_config.txt as listed above, then run again. Nothing was run."
"$PYTHON" "$CODE_DIR/check_configs.py" --quiet "$WORKDIR/__cuts_config.txt" "$WORKDIR/__smearing_config.txt" \
    "$WORKDIR/__input_files_config.txt" --sim-dir "$SIM_DIR" \
    || stop "fix the settings file(s) as listed above, then run again. Nothing was run."
assignments="$("$PYTHON" "$CODE_DIR/input_files.py" shell "$WORKDIR/__input_files_config.txt" "$SIM_DIR")" \
    || stop "could not read __input_files_config.txt"
eval "$assignments"
[ -f "$COMBINED" ] || stop "the combined BIB file $COMBINED does not exist yet - run ./run_all once to build it"
missing=""
for pkg in numpy uproot awkward matplotlib; do
    "$PYTHON" -c "import $pkg" >/dev/null 2>&1 || missing="$missing $pkg"
done
[ -z "$missing" ] || stop "Python ($PYTHON) is missing:$missing   Install with: $PYTHON -m pip install$missing"

STAMP="$(date +%Y-%m-%d_%H%M%S)"
RUN_ID="${STAMP}_scan_$PARAM"
RUN_DIR="$WORKDIR/runs/$RUN_ID"
mkdir -p "$RUN_DIR" || stop "cannot create $RUN_DIR"
cp "$WORKDIR"/__scan_config.txt "$WORKDIR"/__cuts_config.txt "$WORKDIR"/__smearing_config.txt \
   "$WORKDIR"/__input_files_config.txt "$RUN_DIR/" || stop "cannot copy the settings files"
{
    say "Run:      $RUN_ID"
    say "Started:  $(date)"
    say "Code:     $CODE_DIR"
    if git -C "$CODE_DIR" --no-optional-locks rev-parse --short HEAD >/dev/null 2>&1; then
        say "Code version: git commit $(git -C "$CODE_DIR" --no-optional-locks log -1 --format='%h (%ad) %s' --date=short)"
        ch="$(git -C "$CODE_DIR" --no-optional-locks status --porcelain --untracked-files=no)"
        [ -z "$ch" ] && say "No uncommitted code changes." || { say "Uncommitted code changes:"; printf '%s\n' "$ch" | sed 's/^/    /'; }
    fi
    say "Signal:   $SIGNAL"
    say "BIB:      $COMBINED"
} > "$RUN_DIR/code_version.txt"

export MPLBACKEND=Agg
T0=$(date +%s)
{
    say "======================================================================"
    say " Resolution scan $RUN_ID"
    say " Settings (copies kept in runs/$RUN_ID/):"
    "$PYTHON" "$CODE_DIR/check_configs.py" "$RUN_DIR/__cuts_config.txt" "$RUN_DIR/__smearing_config.txt" \
        "$RUN_DIR/__input_files_config.txt" --sim-dir "$SIM_DIR" | sed 's/^/   /'
    grep -v '^#' "$RUN_DIR/__scan_config.txt" | grep -v '^ *$' | sed 's/^/   /'
    say "======================================================================"
    cd "$RUN_DIR" && "$PYTHON" "$CODE_DIR/resolution_scan.py" "$RUN_DIR/__scan_config.txt" \
        "$RUN_DIR/__cuts_config.txt" "$RUN_DIR/__smearing_config.txt" "$SIGNAL" "$COMBINED" "$RUN_DIR"
} 2>&1 | tee -a "$RUN_DIR/run_log.txt"
status=${PIPESTATUS[0]}
el=$(( $(date +%s) - T0 ))

if [ "$status" -ne 0 ] || [ ! -f "$RUN_DIR/scan_$PARAM.pdf" ]; then
    mv "$RUN_DIR" "${RUN_DIR}_FAILED" 2>/dev/null
    say ""
    say "######################################################################"
    say " SCAN DID NOT COMPLETE. Details: runs/${RUN_ID}_FAILED/run_log.txt"
    say "######################################################################"
    exit 1
fi
TOP="_scan_${PARAM}_${STAMP}.pdf"
rm -f -- "$WORKDIR"/_scan_${PARAM}_*.pdf 2>/dev/null
cp "$RUN_DIR/scan_$PARAM.pdf" "$WORKDIR/$TOP"
{
    say ""
    say "======================================================================"
    say " SCAN COMPLETE: runs/$RUN_ID   ($((el / 60))m $((el % 60))s)"
    say " Results: $TOP (also scan_results.csv and the plots in runs/$RUN_ID/)"
    say "======================================================================"
} | tee -a "$RUN_DIR/run_log.txt"
