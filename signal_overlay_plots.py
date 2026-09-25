"""
Step 3: overlay the physics "signal" sample (a muon-gun particle sample,
e.g. ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root - 100k
independent single-muon events) on top of the BIB hit picture, to see
where genuine signal tracks land relative to where BIB background hits
concentrate.

This script only produces the spatial r-z hit map (BIB colored by
subsystem, signal overlaid in black) plus a numeric comparison CSV. An
earlier version also produced three "hit density per subsystem/layer/disk"
bar charts with a signal bar added alongside the BIB mean/peak bars, but
those were dropped: BIB density spans ~5-9 orders of magnitude across
subsystems on its own, and cramming a per-event signal density onto the
same axis stretched the plots further without adding useful information
(the density-only BIB comparison is already in make_basic_plots.py's
own output). See "Normalization note" below for the numbers that remain
useful, still reported in the CSV and the console printout.

Normalization note: the BIB file represents one aggregate background
sample (however many real bunch crossings / whatever time window that
corresponds to - not specified here), while the signal file is 100k
independent single-particle events. There's no known common luminosity
or time normalization between the two, so this script does NOT attempt
to convert either into a physical rate. Instead:
  - BIB density (mean, peak) is reported exactly as in make_basic_plots.py:
    hits/mm^2 for the whole BIB sample.
  - Signal density is reported PER EVENT: total signal hits in a region,
    divided by both area and the number of signal events (n_events) - i.e.
    the average hits/mm^2 a single physics muon leaves in that region.
(For a plotted, hits/collision/bin comparison of the two samples' shapes
in incidence-angle/curvature space, see signal_overlay_angle_plots.py.)

Usage:
    python3 signal_overlay_plots.py <bib.root> <signal.root> [output_dir] [geometry_dir]
"""
import csv
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bib_common import (
    load_hits, region_table, add_peak_density, subsystem_density_table,
    prepare_output_dir, _display_path, SYSTEM_NAMES,
    load_smearing_config, smearing_rng, describe_smearing,
    under_run_all, short_path, loaded_line, print_table,
    apply_position_time_smearing,
    apply_angle_smearing,
)
import geometry as geom_mod

PEAK_PERCENTILE = 99.0
PEAK_KEY = f"peak_density_p{PEAK_PERCENTILE:.0f}_hits_per_mm2"
PEAK_TARGET_HITS_PER_BIN = 20.0

GEOMETRY_FILES = (
    "MuSIC_v2.xml", "Vertex_o2_v06_01.xml",
    "InnerTracker_o2_v07_01.xml", "OuterTracker_o2_v07_01.xml",
)


def region_key(row):
    return (row["system"], row["side"], row["layer"])


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <bib.root> <signal.root> "
              f"[output_dir] [geometry_dir]")
        sys.exit(1)
    bib_file = sys.argv[1]
    signal_file = sys.argv[2]
    outdir = Path(sys.argv[3] if len(sys.argv) > 3 else "../output_signal_overlay")
    geom_dir = Path(sys.argv[4]) if len(sys.argv) > 4 else Path(bib_file).resolve().parent
    outdir = prepare_output_dir(outdir)

    smearing_config = os.environ.get("SMEARING_CONFIG", "").strip()
    smear_cfg = None
    smear_rng = None
    if smearing_config:
        smear_cfg = load_smearing_config(smearing_config)
        smear_rng = smearing_rng(smear_cfg)
        joined = describe_smearing(smear_cfg)
        if not under_run_all():
            print(f"Smearing: {joined}  ({short_path(smearing_config)})")
    bib_hits = load_hits(bib_file)
    if smear_cfg is not None:
        apply_position_time_smearing(bib_hits, smear_cfg, smear_rng)
    n_bib_hit = len(bib_hits["x"])
    print(loaded_line(bib_file, bib_hits, "BIB"))
    sig_hits = load_hits(signal_file)
    if smear_cfg is not None:
        apply_position_time_smearing(sig_hits, smear_cfg, smear_rng)
    n_sig_hit = len(sig_hits["x"])
    n_sig_events = sig_hits["_n_events"]
    print(loaded_line(signal_file, sig_hits, "signal"))

    bib_rows = region_table(bib_hits)
    sig_rows = region_table(sig_hits)

    geometry_used = False
    if all((geom_dir / f).exists() for f in GEOMETRY_FILES):
        try:
            area_lookup = geom_mod.build_area_lookup(geom_dir)
            geom_mod.annotate_rows_with_geometry_area(bib_rows, area_lookup)
            geom_mod.annotate_rows_with_geometry_area(sig_rows, area_lookup)
            geometry_used = True
            print("Sensitive areas from the detector geometry")
        except Exception as e:
            print(f"WARNING: geometry parsing failed ({e}); using hit-inferred area.")
    else:
        print(f"NOTE: geometry files not found in {short_path(geom_dir)}; "
              f"using hit-inferred area.")

    add_peak_density(bib_hits, bib_rows, bin_size_mm=None, percentile=PEAK_PERCENTILE,
                      target_hits_per_bin=PEAK_TARGET_HITS_PER_BIN)
    bib_sys_rows = subsystem_density_table(bib_rows)
    sig_sys_rows = subsystem_density_table(sig_rows)  # mean density only, no peak
    sys_ids = sorted(SYSTEM_NAMES.keys())

    bib_by_sys = {r["system"]: r for r in bib_sys_rows}
    sig_by_sys = {r["system"]: r for r in sig_sys_rows}

    # ---- r-z map: BIB hits (by subsystem) + signal hits overlaid --------
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = {1: "#e6194B", 2: "#f58231", 3: "#3cb44b",
              4: "#4363d8", 5: "#911eb4", 6: "#42d4f4"}
    rng = np.random.default_rng(0)
    bib_system = bib_hits["system"]
    for s in sys_ids:
        mask = bib_system == s
        idx = np.where(mask)[0]
        if len(idx) > 40000:
            idx = rng.choice(idx, 40000, replace=False)
        ax.scatter(bib_hits["z"][idx], bib_hits["r"][idx], s=0.5,
                   color=colors[s], label=f"BIB {SYSTEM_NAMES[s]}", alpha=0.5)
    # signal hits on top, subsampled, distinct marker
    n_sig_plot = min(len(sig_hits["x"]), 60000)
    sig_idx = rng.choice(len(sig_hits["x"]), n_sig_plot, replace=False) \
        if len(sig_hits["x"]) > n_sig_plot else np.arange(len(sig_hits["x"]))
    ax.scatter(sig_hits["z"][sig_idx], sig_hits["r"][sig_idx], s=1.5,
               color="black", label="signal hits", alpha=0.6, marker=".")
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("r (mm)")
    ax.set_title("Hit map (r vs z): BIB (colored by subsystem, subsampled) + "
                 "signal hits overlaid (black)")
    ax.legend(markerscale=15, loc="upper right", fontsize=7)
    plt.tight_layout()
    plt.savefig(outdir / "rz_map_with_signal.png", dpi=140)
    plt.close(fig)

    # ---- comparison CSV -------------------------------------------------
    csv_path = outdir / "signal_vs_bib_summary.csv"
    fields = ["system", "system_name", "bib_n_hits", "bib_mean_density_hits_per_mm2",
              f"bib_{PEAK_KEY}", "signal_n_hits", "signal_n_events",
              "signal_mean_density_hits_per_mm2_per_event"]
    with open(csv_path, "w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=fields, restval="")
        w_.writeheader()
        for s in sys_ids:
            b = bib_by_sys[s]
            sg = sig_by_sys.get(s, {"n_hits": 0, "density_hits_per_mm2": 0.0})
            w_.writerow({
                "system": s, "system_name": SYSTEM_NAMES[s],
                "bib_n_hits": b["n_hits"],
                "bib_mean_density_hits_per_mm2": b["density_hits_per_mm2"],
                f"bib_{PEAK_KEY}": b.get(PEAK_KEY, ""),
                "signal_n_hits": sg["n_hits"],
                "signal_n_events": n_sig_events,
                "signal_mean_density_hits_per_mm2_per_event":
                    sg["density_hits_per_mm2"] / n_sig_events,
            })

    print(f"Wrote rz_map_with_signal.png and {csv_path.name} to {short_path(outdir)}/")
    density_rows = []
    for s in sys_ids:
        b = bib_by_sys[s]
        sg = sig_by_sys.get(s, {"density_hits_per_mm2": 0.0})
        sd = sg["density_hits_per_mm2"] / n_sig_events
        density_rows.append([SYSTEM_NAMES[s], f"{b['density_hits_per_mm2']:.4g}",
                             f"{b[PEAK_KEY]:.4g}", f"{sd:.4g}"])
    print()
    print_table("Hit density (hits/mm^2): BIB (whole sample) vs. signal (per event):",
                ["subsystem", "BIB mean", "BIB peak", "signal per event"], density_rows)


if __name__ == "__main__":
    main()
