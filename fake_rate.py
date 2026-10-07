"""
fake_rate.py - expected number of fake tracks, E[#fakes], in the 6-layer
IT + OT barrel tower of the fake-rate framework (Sections 14-16 of the
paper; calibration in fake_rate_framework/real_tower_calibration.py),
from the BIB hit density of each IT/OT barrel layer after the cuts.

E[#fakes] = P_true x prod_i n_i,  n_i = (BIB hits after cuts / area)_i x W_i^2
with W_i^2 the plane areas of the projective tower (1 m^2 outermost plane)
and P_true the calibrated probability for 3 track models.

Used by resolution_scan.py (./run_scan), and by ./run_all at the end of
its cuts step:
    python3 fake_rate.py <step4_cuts folder>
reads density_per_layer_before_after_cuts.csv there, prints E[#fakes] and
writes fake_rate.csv next to it.
"""
import csv
import sys
from pathlib import Path

IT_BARREL, OT_BARREL = 3, 5
TOWER = (IT_BARREL, OT_BARREL)
LAYERS = (0, 1, 2)

# Real per-layer sensitive areas (mm^2) of the IT/OT barrel layers
AREA_MM2 = {
    (IT_BARREL, 0): 1043723.52, (IT_BARREL, 1): 2203416.32, (IT_BARREL, 2): 5167881.04,
    (OT_BARREL, 0): 14003290.56, (OT_BARREL, 1): 19482839.04, (OT_BARREL, 2): 28615419.84,
}
# Plane areas W_i^2 (mm^2) of the projective tower (1 m^2 outermost plane)
Y_TRUE = {(IT_BARREL, 0): 164.0, (IT_BARREL, 1): 354.0, (IT_BARREL, 2): 554.0,
          (OT_BARREL, 0): 819.0, (OT_BARREL, 1): 1153.0, (OT_BARREL, 2): 1486.0}
W2_SYNTH = {k: (1000.0 * y / 1486.0) ** 2 for k, y in Y_TRUE.items()}
# Calibrated fake probability P_true for this tower, per track model
P_TRUE = {"exact_helix": 7.3054e-26, "conservative_quad": 6.2186e-26, "line_B0": 3.6895e-28}
MODEL_LABEL = {"exact_helix": "exact helix", "conservative_quad": "conservative quad",
               "line_B0": "line (B=0)"}


def efakes(n_after, area=AREA_MM2):
    """n_after: {(system, layer): BIB hits after the cuts} for the 6 tower
    layers. Returns (prod_n, {model: E[#fakes]})."""
    prod_n = 1.0
    for key in AREA_MM2:
        prod_n *= n_after[key] / area[key] * W2_SYNTH[key]
    return prod_n, {m: p * prod_n for m, p in P_TRUE.items()}


def summary_line(e):
    return ("Expected fake tracks E[#fakes], 6-layer IT+OT barrel tower: "
            + ",  ".join(f"{MODEL_LABEL[m]} {e[m]:.2e}" for m in P_TRUE))


def main(argv):
    if len(argv) != 1:
        print(__doc__.split("Used by")[1], file=sys.stderr)
        return 2
    folder = Path(argv[0])
    n_after, area = {}, {}
    with open(folder / "density_per_layer_before_after_cuts.csv") as f:
        for r in csv.DictReader(f):
            key = (int(r["system"]), int(r["layer"]))
            if key in AREA_MM2:
                n_after[key] = int(r["bib_n_hits_after"])
                area[key] = float(r["area_mm2"])
    missing = [k for k in AREA_MM2 if k not in n_after]
    if missing:
        print(f"ERROR: no IT/OT barrel layer(s) {missing} in the density table", file=sys.stderr)
        return 1
    for k in AREA_MM2:      # the calibration assumes these exact layer areas
        if abs(area[k] / AREA_MM2[k] - 1) > 1e-3:
            print(f"WARNING: layer {k} area {area[k]:.0f} mm^2 differs from the calibration's "
                  f"{AREA_MM2[k]:.0f} mm^2 (different geometry?)", file=sys.stderr)
    prod_n, e = efakes(n_after, area)
    with open(folder / "fake_rate.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["prod_n"] + [f"efakes_{m}" for m in P_TRUE]
                   + [f"bib_hits_after_{'IT' if s == 3 else 'OT'}_barrel_L{l}" for s, l in AREA_MM2])
        w.writerow([prod_n] + [e[m] for m in P_TRUE] + [n_after[k] for k in AREA_MM2])
    print("BIB hits after the cuts, IT/OT barrel layers:  "
          + "  ".join(f"{'IT' if s == 3 else 'OT'} L{l} {n_after[(s, l)]:,}" for s, l in AREA_MM2))
    print(summary_line(e))
    print("  (density after cuts of each layer x tower plane area, times the calibrated "
          "fake probability; Sections 14-16)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
