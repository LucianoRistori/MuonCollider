#!/bin/bash
# Re-run every script that depends on cuts_config.txt (step 4: BIB
# rejection/cutflow, tracking efficiency vs pT, and N-1 cut-validation
# plots), so after editing the cut values a single command regenerates
# everything they affect.
#
# Every run is archived as one self-contained, timestamped folder under
# results/step4_runs/ - together with a copy of the cuts_config.txt that
# produced it, so old hypotheses are never lost or overwritten and each
# run's plots/tables/cuts stay grouped together as a unit. The familiar
# fixed-name folders (results/step4_cuts/, results/step4_track_efficiency/,
# results/step4_n1_cuts/) are also refreshed each time as an always-
# "latest" convenience copy, for anything (e.g. the project doc) that
# just wants the current results without digging into results/step4_runs/.
#
# Everything printed to the terminal during the run is also saved as
# run_log.txt alongside the other archived files (and mirrored to
# results/step4_run_log.txt as the "latest" copy). Most of it duplicates
# the CSVs (the cutflow and track-efficiency tables are printed in full,
# not just summarized), but it also captures a couple of things that
# aren't written to any CSV - notably apply_cuts.py's final "BIB
# rejection vs. signal efficiency" summary - and it's a lot faster to
# skim than opening several CSVs.
#
# Usage:
#   ./run_step4.sh
#
# Edit cuts_config.txt first, then run this - no arguments needed for
# the normal case. Paths below assume the standard layout: this script
# lives next to the other analysis scripts (e.g. ~/code/MuonCollider),
# and the ROOT files / results live under DATA_DIR (e.g.
# ~/Dropbox/Documents/MuonColliderSimulation). Override DATA_DIR,
# BIB_FILE or SIGNAL_FILE as environment variables if a run ever needs
# different inputs, e.g.:
#   SIGNAL_FILE=/path/to/other_signal.root ./run_step4.sh

set -e
set -o pipefail
cd "$(dirname "$0")"

DATA_DIR="${DATA_DIR:-$HOME/Dropbox/Documents/MuonColliderSimulation}"
BIB_FILE="${BIB_FILE:-$DATA_DIR/ntu_bib_2evt.root}"
SIGNAL_FILE="${SIGNAL_FILE:-$DATA_DIR/ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root}"
CUTS_CONFIG="${CUTS_CONFIG:-cuts_config.txt}"

TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
RUN_DIR="$DATA_DIR/results/step4_runs/$TIMESTAMP"
mkdir -p "$RUN_DIR"
cp "$CUTS_CONFIG" "$RUN_DIR/cuts_config.txt"

# Everything from here on is both shown on the terminal and saved to
# run_log.txt (via the `tee` at the very end of this script), by running
# the rest of the work inside this function.
run_all() {
    echo "Using cuts: $CUTS_CONFIG"
    echo "  BIB file:    $BIB_FILE"
    echo "  Signal file: $SIGNAL_FILE"
    echo "Archiving this run to: $RUN_DIR"
    echo

    echo "=== Step 4a: apply_cuts.py (BIB rejection / cutflow, with vs. without cuts) ==="
    python3 apply_cuts.py "$BIB_FILE" "$SIGNAL_FILE" "$CUTS_CONFIG" \
        "$RUN_DIR/cuts"

    echo
    echo "=== Step 4b: track_efficiency.py (tracking efficiency vs. generated pT) ==="
    python3 track_efficiency.py "$SIGNAL_FILE" --cuts "$CUTS_CONFIG" \
        --out "$RUN_DIR/track_efficiency"

    echo
    echo "=== Step 4c: n1_cut_plots.py (N-1 cut-validation plots) ==="
    python3 n1_cut_plots.py "$BIB_FILE" "$SIGNAL_FILE" "$CUTS_CONFIG" \
        "$RUN_DIR/n1_cuts"

    echo
    echo "Done."
    echo "Archived (this run, permanent): $RUN_DIR/"
    echo "  cuts_config.txt   (the cuts used for this run)"
    echo "  cuts/             (cutflow_bib.csv, cutflow_signal.csv,"
    echo "                     density_before_after_cuts.csv/.png)"
    echo "  track_efficiency/ (track_efficiency_vs_pt.csv/.png)"
    echo "  n1_cuts/          (N-1 z-intercept/pT/time plots, full+zoom, 6 PNGs)"
    echo "  run_log.txt       (everything printed during this run)"
}

run_all 2>&1 | tee "$RUN_DIR/run_log.txt"

# Refresh the "latest" convenience copies at the familiar fixed paths.
# The archive above (RUN_DIR) is the permanent record of this run; these
# are just a mirror of it for anything that wants "the current results"
# without a timestamp in the path.
LATEST_CUTS="$DATA_DIR/results/step4_cuts"
LATEST_TRACK="$DATA_DIR/results/step4_track_efficiency"
LATEST_N1="$DATA_DIR/results/step4_n1_cuts"
rm -rf "$LATEST_CUTS" "$LATEST_TRACK" "$LATEST_N1"
cp -r "$RUN_DIR/cuts" "$LATEST_CUTS"
cp -r "$RUN_DIR/track_efficiency" "$LATEST_TRACK"
cp -r "$RUN_DIR/n1_cuts" "$LATEST_N1"
cp "$RUN_DIR/run_log.txt" "$DATA_DIR/results/step4_run_log.txt"

echo
echo "Updated (latest convenience copy):"
echo "  $LATEST_CUTS/"
echo "  $LATEST_TRACK/"
echo "  $LATEST_N1/"
echo "  $DATA_DIR/results/step4_run_log.txt"
