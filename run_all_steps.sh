#!/bin/bash
# =====================================================================
# Full BIB analysis, steps 1-4, run from a WORKING FOLDER that holds the
# three editable settings files. Normally started with the ./run_all
# launcher that lives in that folder:
#
#     cd ~/Dropbox/Documents/MuonColliderSimulation/Analysis
#     (edit __cuts_config.txt, __smearing_config.txt, __input_files_config.txt)
#     ./run_all
#
# What one run does:
#   1. Checks everything first - settings files (typos, bad values),
#      input files, Python packages. If anything is wrong it stops
#      before running anything.
#   2. Creates runs/<date>_<time>/ in the working folder and copies the
#      settings files into it at the START, so the archive records
#      exactly what was used (editing the files during a run has no
#      effect on that run).
#   3. Runs steps 1-4, writing every plot and table into that run folder,
#      and saves everything printed to runs/<date>_<time>/run_log.txt.
#   4. Only if ALL steps succeed: writes summary.txt (settings + BIB
#      rejection / signal efficiency) and replaces the step* folders in
#      the working folder with this run's copy, so they always show the
#      latest SUCCESSFUL run (latest_run.txt says which one).
#      If anything fails, the run folder is renamed <date>_<time>_FAILED
#      and the step* folders are left as they were.
#
# Input files: named in __input_files_config.txt, and read from the Data/
# (ROOT) and Geometry/ (XML) folders inside the folder ABOVE the working
# folder (override that folder with SIM_DIR=...). If the combined BIB file
# is missing, or doesn't contain exactly the plus/minus/ipp files listed,
# it is (re)built first. Python is python3 (override: PYTHON=...).
#
# Splitting a run into pieces (only needed under a time limit):
#   ./run_all --steps 1,2                   start a run with steps 1-2
#   ./run_all --run <date>_<time> --steps 3,4   add steps 3-4 to it
# The run is completed automatically once all four steps are done.
# =====================================================================
set -u
set -o pipefail

CODE_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PYTHON="${PYTHON:-python3}"
STEP_DIRS="step1_basic_plots step1_basic_plots_minus step1_basic_plots_combined
step2_incidence_angles step2_incidence_angles_minus step2_incidence_angles_combined
step2_time_of_flight step2_time_of_flight_minus step2_time_of_flight_combined
step2_time_of_flight_signal step3_signal_overlay
step4_cuts step4_track_efficiency step4_n1_cuts"

# The key plots, gathered into _highlights/ after every successful run
# (in the run folder, and in the working folder next to the step folders).
HIGHLIGHTS="step4_cuts/density_before_after_cuts.png
step4_track_efficiency/track_efficiency_vs_pt.png
step3_signal_overlay/inv_radius_per_subsystem_with_signal.png
step3_signal_overlay/inv_radius_per_subsystem_zoom_with_signal.png
step3_signal_overlay/time_corrected_per_subsystem_with_signal.png
step3_signal_overlay/time_corrected_per_subsystem_zoom_with_signal.png
step3_signal_overlay/z_axis_intercept_per_subsystem_with_signal.png
step3_signal_overlay/z_axis_intercept_per_subsystem_zoom_with_signal.png
step4_n1_cuts/inv_radius_per_subsystem_n1.png
step4_n1_cuts/inv_radius_per_subsystem_zoom_n1.png
step4_n1_cuts/time_corrected_per_subsystem_n1.png
step4_n1_cuts/time_corrected_per_subsystem_zoom_n1.png
step4_n1_cuts/z_axis_intercept_per_subsystem_n1.png
step4_n1_cuts/z_axis_intercept_per_subsystem_zoom_n1.png"

tilde() {   # show a path with the home folder written as ~
    case "$1" in "$HOME"/*) printf '~/%s' "${1#"$HOME"/}" ;; *) printf '%s' "$1" ;; esac
}
say() { printf '%s\n' "$*"; }
stop() { printf '\nERROR: %s\n\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- arguments
WORKDIR=""
STEPS="1 2 3 4"
RUN_ID=""
while [ $# -gt 0 ]; do
    case "$1" in
        --steps) [ $# -ge 2 ] || stop "--steps needs a value, e.g. --steps 1,2"
                 STEPS="$(printf '%s' "$2" | tr ',' ' ')"; shift 2 ;;
        --run)   [ $# -ge 2 ] || stop "--run needs a run name (a folder name under runs/)"
                 RUN_ID="$2"; shift 2 ;;
        -h|--help) sed -n '2,/^# =====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*) stop "unknown option '$1' (try --help)" ;;
        *)  [ -z "$WORKDIR" ] || stop "unexpected argument '$1' (try --help)"
            WORKDIR="$1"; shift ;;
    esac
done
for s in $STEPS; do
    case "$s" in 1|2|3|4) ;; *) stop "--steps: '$s' is not a step number (1 to 4)" ;; esac
done
[ -n "$WORKDIR" ] || WORKDIR="$PWD"
[ -d "$WORKDIR" ] || stop "working folder not found: $WORKDIR"
WORKDIR="$(cd "$WORKDIR" && pwd -P)"

SIM_DIR="${SIM_DIR:-${DATA_DIR:-$(dirname "$WORKDIR")}}"   # holds Data/ and Geometry/

# ---------------------------------------------------------------- checks
command -v "$PYTHON" >/dev/null 2>&1 || stop "'$PYTHON' not found. Install Python 3, or run with PYTHON=/path/to/python3 ./run_all"

if [ -z "$RUN_ID" ]; then
    for f in __cuts_config.txt __smearing_config.txt __input_files_config.txt; do
        if [ ! -f "$WORKDIR/$f" ]; then
            if [ -f "$CODE_DIR/templates/$f" ]; then
                cp "$CODE_DIR/templates/$f" "$WORKDIR/$f"
                stop "$f was missing, so a fresh copy was made from the template:
       $WORKDIR/$f
       Check its values, then run again."
            fi
            stop "$f not found in $WORKDIR"
        fi
    done
    CFG_DIR="$WORKDIR"
else
    case "$RUN_ID" in
        *_FAILED|*_INTERRUPTED) stop "run $RUN_ID did not complete - start a new run instead" ;;
    esac
    [ -d "$WORKDIR/runs/$RUN_ID" ] || stop "no run named $RUN_ID in $WORKDIR/runs/"
    [ ! -f "$WORKDIR/runs/$RUN_ID/summary.txt" ] || stop "run $RUN_ID is already complete"
    CFG_DIR="$WORKDIR/runs/$RUN_ID"
fi

[ -f "$CFG_DIR/__input_files_config.txt" ] || stop "__input_files_config.txt not found in $CFG_DIR"
"$PYTHON" "$CODE_DIR/check_configs.py" --quiet "$CFG_DIR/__cuts_config.txt" "$CFG_DIR/__smearing_config.txt" \
    "$CFG_DIR/__input_files_config.txt" --sim-dir "$SIM_DIR" \
    || stop "fix the settings file(s) as listed above, then run again. Nothing was run."
assignments="$("$PYTHON" "$CODE_DIR/input_files.py" shell "$CFG_DIR/__input_files_config.txt" "$SIM_DIR")" \
    || stop "could not read __input_files_config.txt"
eval "$assignments"   # PLUS MINUS IPP COMBINED SIGNAL GEOM_MAIN GEOM_VERTEX GEOM_IT GEOM_OT

missing=""
for pkg in numpy uproot awkward matplotlib; do
    "$PYTHON" -c "import $pkg" >/dev/null 2>&1 || missing="$missing $pkg"
done
[ -z "$missing" ] || stop "Python ($PYTHON) is missing these packages:$missing
       Install them with:   $PYTHON -m pip install$missing"

# ---------------------------------------------------------------- run folder
if [ -z "$RUN_ID" ]; then
    RUN_ID="$(date +%Y-%m-%d_%H%M%S)"
    while [ -e "$WORKDIR/runs/$RUN_ID" ] || [ -e "$WORKDIR/runs/${RUN_ID}_FAILED" ]; do
        sleep 1; RUN_ID="$(date +%Y-%m-%d_%H%M%S)"
    done
    RUN_DIR="$WORKDIR/runs/$RUN_ID"
    mkdir -p "$RUN_DIR" || stop "cannot create $RUN_DIR"
    cp "$WORKDIR/__cuts_config.txt" "$WORKDIR/__smearing_config.txt" "$WORKDIR/__input_files_config.txt" "$RUN_DIR/" \
        || stop "cannot copy the settings files into $RUN_DIR"
    {
        say "Run:            $RUN_ID"
        say "Started:        $(date)"
        say "Working folder: $WORKDIR"
        say "Code:           $CODE_DIR"
        if git -C "$CODE_DIR" --no-optional-locks rev-parse --short HEAD >/dev/null 2>&1; then
            say "Code version:   git commit $(git -C "$CODE_DIR" --no-optional-locks log -1 --format='%h (%ad) %s' --date=short)"
            changes="$(git -C "$CODE_DIR" --no-optional-locks status --porcelain --untracked-files=no)"
            if [ -n "$changes" ]; then
                say "Uncommitted code changes at start of run:"
                printf '%s\n' "$changes" | sed 's/^/    /'
            else
                say "No uncommitted code changes."
            fi
        fi
        say "Python:         $("$PYTHON" -c 'import sys; print(sys.version.split()[0])') ($(command -v "$PYTHON"))"
    } > "$RUN_DIR/code_version.txt"
else
    RUN_DIR="$WORKDIR/runs/$RUN_ID"
fi

CUTS="$RUN_DIR/__cuts_config.txt"
export SMEARING_CONFIG="$RUN_DIR/__smearing_config.txt"
export BIB_GEOMETRY_FILES="$GEOM_MAIN:$GEOM_VERTEX:$GEOM_IT:$GEOM_OT"
export BIB_RUN_ALL=1
export MPLBACKEND=Agg

# ---------------------------------------------------------------- the steps
run_py() {   # run_py "<label>" <script.py> <args...>
    local label="$1" script="$2" t0
    shift 2
    say ""
    say "---- $label"
    t0=$(date +%s)
    if ! "$PYTHON" "$CODE_DIR/$script" "$@"; then
        say "!!!! FAILED: $label"
        printf '%s\n' "$label" > "$RUN_DIR/FAILED_STEP.txt"
        return 1
    fi
    say "     ($(( $(date +%s) - t0 ))s)"
}

step1() {
    run_py "Step 1: basic plots, BIB plus"     make_basic_plots.py "$PLUS"     "$RUN_DIR/step1_basic_plots" &&
    run_py "Step 1: basic plots, BIB minus"    make_basic_plots.py "$MINUS"    "$RUN_DIR/step1_basic_plots_minus" &&
    run_py "Step 1: basic plots, BIB combined" make_basic_plots.py "$COMBINED" "$RUN_DIR/step1_basic_plots_combined"
}
step2() {
    run_py "Step 2: incidence angles, BIB plus"     incidence_angle_plots.py "$PLUS"     "$RUN_DIR/step2_incidence_angles" &&
    run_py "Step 2: incidence angles, BIB minus"    incidence_angle_plots.py "$MINUS"    "$RUN_DIR/step2_incidence_angles_minus" &&
    run_py "Step 2: incidence angles, BIB combined" incidence_angle_plots.py "$COMBINED" "$RUN_DIR/step2_incidence_angles_combined" &&
    run_py "Step 2: time of flight, BIB plus"       time_of_flight_plots.py  "$PLUS"     "$RUN_DIR/step2_time_of_flight" &&
    run_py "Step 2: time of flight, BIB minus"      time_of_flight_plots.py  "$MINUS"    "$RUN_DIR/step2_time_of_flight_minus" &&
    run_py "Step 2: time of flight, BIB combined"   time_of_flight_plots.py  "$COMBINED" "$RUN_DIR/step2_time_of_flight_combined" &&
    run_py "Step 2: time of flight, signal"         time_of_flight_plots.py  "$SIGNAL"   "$RUN_DIR/step2_time_of_flight_signal"
}
step3() {
    run_py "Step 3: signal overlay - r-z map and densities" \
        signal_overlay_plots.py "$COMBINED" "$SIGNAL" "$RUN_DIR/step3_signal_overlay" &&
    run_py "Step 3: signal overlay - z-intercept, pT, time" \
        signal_overlay_angle_plots.py "$COMBINED" "$SIGNAL" "$RUN_DIR/step3_signal_overlay" "$CUTS"
}
step4() {
    run_py "Step 4: cuts - BIB rejection and signal efficiency" \
        apply_cuts.py "$COMBINED" "$SIGNAL" "$CUTS" "$RUN_DIR/step4_cuts" &&
    run_py "Step 4: track-finding efficiency vs pT" \
        track_efficiency.py "$SIGNAL" --cuts "$CUTS" --out "$RUN_DIR/step4_track_efficiency" &&
    run_py "Step 4: N-1 cut plots" \
        n1_cut_plots.py "$COMBINED" "$SIGNAL" "$CUTS" "$RUN_DIR/step4_n1_cuts"
}

main() {
    say "======================================================================"
    say " BIB analysis run $RUN_ID   (steps: $STEPS)"
    say " Working folder: $(tilde "$WORKDIR")"
    say " Settings (copies kept in runs/$RUN_ID/):"
    "$PYTHON" "$CODE_DIR/check_configs.py" "$CUTS" "$SMEARING_CONFIG" \
        "$RUN_DIR/__input_files_config.txt" --sim-dir "$SIM_DIR" | sed 's/^/   /'
    say "======================================================================"
    # run the scripts from inside the run folder, so the output folders
    # they report show up as short relative paths (step1_basic_plots/ ...)
    cd "$RUN_DIR" || return 1
    run_py "Combined BIB file" input_files.py prepare "$RUN_DIR/__input_files_config.txt" "$SIM_DIR" \
        || return 1
    if ! grep -q "^Input files" "$RUN_DIR/code_version.txt" 2>/dev/null; then
        {
            say "Input files (bytes, modified):"
            for f in "$PLUS" "$MINUS" ${IPP:+"$IPP"} "$COMBINED" "$SIGNAL" \
                     "$GEOM_MAIN" "$GEOM_VERTEX" "$GEOM_IT" "$GEOM_OT"; do
                say "    $(ls -l "$f" | awk '{print $5", "$6" "$7" "$8}')  $(tilde "$f")"
            done
        } >> "$RUN_DIR/code_version.txt"
    fi
    local s t0
    for s in $STEPS; do
        if grep -qx "step $s" "$RUN_DIR/steps_done.txt" 2>/dev/null; then
            say ""
            say "(step $s was already done in this run - skipping it)"
            continue
        fi
        t0=$(date +%s)
        "step$s" || return 1
        say "step $s" >> "$RUN_DIR/steps_done.txt"
        say ""
        say "==== step $s done ($(( $(date +%s) - t0 ))s)"
    done
}

T_START=$(date +%s)
interrupted=0
trap 'interrupted=1' INT TERM
main 2>&1 | tee -i -a "$RUN_DIR/run_log.txt"
status=${PIPESTATUS[0]}
trap - INT TERM
if [ "$interrupted" = 1 ] && [ "$status" = 0 ]; then status=130; fi

latest_note() {
    if [ -f "$WORKDIR/latest_run.txt" ]; then
        say " The step folders here still show the previous successful run ($(sed -n 's/^The step folders in this folder show run: //p' "$WORKDIR/latest_run.txt"))."
    else
        say " The step folders here were not changed."
    fi
}

# ---------------------------------------------------------------- failure
if [ "$status" -ne 0 ]; then
    suffix="_FAILED"
    [ "$interrupted" = 1 ] && suffix="_INTERRUPTED"
    what="$(cat "$RUN_DIR/FAILED_STEP.txt" 2>/dev/null || echo "stopped")"
    mv "$RUN_DIR" "$RUN_DIR$suffix" 2>/dev/null && RUN_DIR="$RUN_DIR$suffix"
    {
        say ""
        say "######################################################################"
        say " RUN DID NOT COMPLETE: $what"
        say " Details:  runs/$(basename "$RUN_DIR")/run_log.txt"
        latest_note
        say "######################################################################"
    } | tee -a "$RUN_DIR/run_log.txt"
    exit 1
fi

# ---------------------------------------------------------------- partial run
todo=""
for s in 1 2 3 4; do
    grep -qx "step $s" "$RUN_DIR/steps_done.txt" 2>/dev/null || todo="$todo $s"
done
if [ -n "$todo" ]; then
    say ""
    say "Run $RUN_ID is not finished yet - steps still to do:$todo"
    say "Continue it with:   ./run_all --run $RUN_ID --steps $(echo $todo | tr ' ' ',')"
    exit 0
fi

# ---------------------------------------------------------------- finish
{
    say "Run $RUN_ID"
    "$PYTHON" "$CODE_DIR/check_configs.py" --brief "$RUN_DIR/__cuts_config.txt" "$RUN_DIR/__smearing_config.txt" \
        "$RUN_DIR/__input_files_config.txt" --sim-dir "$SIM_DIR"
    say ""
    awk '/Summary: BIB rejection vs. signal efficiency/ {print; f=1; next} f && /^  [^ ]/ {print; next} f {exit}' "$RUN_DIR/run_log.txt"
} > "$RUN_DIR/summary.txt"

missing_hl=""
mkdir -p "$RUN_DIR/_highlights"
for f in $HIGHLIGHTS; do
    if [ -f "$RUN_DIR/$f" ]; then
        cp "$RUN_DIR/$f" "$RUN_DIR/_highlights/"
    else
        missing_hl="$missing_hl $f"
    fi
done

promote() {
    local d
    for d in $STEP_DIRS _highlights; do
        [ -d "$RUN_DIR/$d" ] || { say "missing output folder: runs/$RUN_ID/$d"; return 1; }
    done
    say "The step folders are being updated from run $RUN_ID - if you can read this, the update did not finish; the complete results are in runs/$RUN_ID/" > "$WORKDIR/latest_run.txt"
    for d in $STEP_DIRS _highlights; do
        rm -rf -- "${WORKDIR:?}/${d:?}" 2>/dev/null
        [ ! -e "$WORKDIR/$d" ] || { say "could not remove the old $d folder"; return 1; }
        cp -R "$RUN_DIR/$d" "$WORKDIR/$d" || return 1
    done
    {
        say "The step folders in this folder show run: $RUN_ID"
        say "(so does _highlights/; the complete archive of that run, with its settings and log, is runs/$RUN_ID/)"
        say ""
        cat "$RUN_DIR/summary.txt"
    } > "$WORKDIR/latest_run.txt"
}

elapsed=$(( $(date +%s) - T_START ))
{
    say ""
    say "======================================================================"
    if promote; then
        say " RUN COMPLETE: runs/$RUN_ID   ($((elapsed / 60))m $((elapsed % 60))s)"
        say " The step folders here now show this run; key plots are in _highlights/."
    else
        say " RUN COMPLETE: runs/$RUN_ID   ($((elapsed / 60))m $((elapsed % 60))s)"
        say " BUT the step folders here could not be updated (see message above);"
        say " the complete results are in runs/$RUN_ID/."
    fi
    [ -z "$missing_hl" ] || say " WARNING - these plots were not produced, so they are not in _highlights:$missing_hl"
    say ""
    sed 's/^/ /' "$RUN_DIR/summary.txt"
    say "======================================================================"
} 2>&1 | tee -a "$RUN_DIR/run_log.txt"
