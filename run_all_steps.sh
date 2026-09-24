#!/bin/bash
# Runs steps 1-4 of the BIB analysis pipeline end to end, against the
# current data files in DATA_DIR. Steps 1-3 don't have their own driver
# (unlike step 4, which already has run_step4.sh), so this covers them
# with the standard per-file / per-output-dir invocations described in
# the project README, then calls run_step4.sh for step 4.
set -e
set -o pipefail
cd "$(dirname "$0")"

DATA_DIR="${DATA_DIR:-$HOME/Dropbox/Documents/MuonColliderSimulation}"
PLUS="$DATA_DIR/ntu_bib_plus_1evt.root"
MINUS="$DATA_DIR/ntu_bib_minus_1evt.root"
COMBINED="$DATA_DIR/ntu_bib_ipp_3evt.root"
SIGNAL="$DATA_DIR/ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root"
RESULTS="$DATA_DIR/results"

echo "DATA_DIR:  $DATA_DIR"
echo "PLUS:      $PLUS"
echo "MINUS:     $MINUS"
echo "COMBINED:  $COMBINED  (plus+minus+ipp)"
echo "SIGNAL:    $SIGNAL"
echo

echo "=== Step 1: make_basic_plots.py ==="
python3 make_basic_plots.py "$PLUS"     "$RESULTS/step1_basic_plots"
python3 make_basic_plots.py "$MINUS"    "$RESULTS/step1_basic_plots_minus"
python3 make_basic_plots.py "$COMBINED" "$RESULTS/step1_basic_plots_combined"

echo
echo "=== Step 2a: incidence_angle_plots.py ==="
python3 incidence_angle_plots.py "$PLUS"     "$RESULTS/step2_incidence_angles"
python3 incidence_angle_plots.py "$MINUS"    "$RESULTS/step2_incidence_angles_minus"
python3 incidence_angle_plots.py "$COMBINED" "$RESULTS/step2_incidence_angles_combined"

echo
echo "=== Step 2b: time_of_flight_plots.py ==="
python3 time_of_flight_plots.py "$PLUS"     "$RESULTS/step2_time_of_flight"
python3 time_of_flight_plots.py "$MINUS"    "$RESULTS/step2_time_of_flight_minus"
python3 time_of_flight_plots.py "$COMBINED" "$RESULTS/step2_time_of_flight_combined"
python3 time_of_flight_plots.py "$SIGNAL"   "$RESULTS/step2_time_of_flight_signal"

echo
echo "=== Step 3: signal_overlay_plots.py + signal_overlay_angle_plots.py ==="
python3 signal_overlay_plots.py "$COMBINED" "$SIGNAL" "$RESULTS/step3_signal_overlay"
python3 signal_overlay_angle_plots.py "$COMBINED" "$SIGNAL" "$RESULTS/step3_signal_overlay"

echo
echo "=== Step 4: run_step4.sh (apply_cuts.py, track_efficiency.py, n1_cut_plots.py) ==="
./run_step4.sh

echo
echo "ALL STEPS COMPLETE"
