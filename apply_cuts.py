"""
Step 4: apply the selection cuts defined in cuts_config.txt (time-of-
flight-corrected hit time, z-axis intercept, curvature-based momentum -
see bib_common.load_cuts / bib_common.apply_cuts) to both the BIB and
signal samples, and compare - the main point of this script - hit
density and hit counts WITH cuts vs. WITHOUT cuts, per subsystem, so a
cut hypothesis can be judged by how much it suppresses BIB relative to
how much signal it keeps.

Two kinds of comparison are produced:
  1. A cutflow (hit counts before / after each individual cut / after
     all three combined, per subsystem) for BIB and for signal
     separately - this shows each cut's own selectivity.
  2. A before/after hit-density comparison (mean density, hits/mm^2,
     using the true detector geometry area when available - same
     definition as step 1/3) for BIB only. (Signal hit density isn't
     computed here - it isn't a meaningful quantity for a single-track
     sample; see track_efficiency.py instead for the signal-side
     question this cut set is meant to answer: what fraction of muon
     tracks are still found after the cuts.)

Usage:
    python3 apply_cuts.py <bib.root> <signal.root> [cuts_config] [output_dir] [geometry_dir]
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bib_common import (
    load_hits, add_incidence_angles, add_time_of_flight, load_cuts, apply_cuts,
    mask_hits, region_table, add_peak_density, subsystem_density_table,
    prepare_output_dir, _display_path, SYSTEM_NAMES,
)
import geometry as geom_mod

PEAK_PERCENTILE = 99.0
PEAK_KEY = f"peak_density_p{PEAK_PERCENTILE:.0f}_hits_per_mm2"
PEAK_TARGET_HITS_PER_BIN = 20.0
GEOMETRY_FILES = (
    "MuSIC_v2.xml", "Vertex_o2_v06_01.xml",
    "InnerTracker_o2_v07_01.xml", "OuterTracker_o2_v07_01.xml",
)
CUT_ORDER = ("t_corrected_ns", "z_axis_intercept_mm", "momentum_gev")
BIB_COLOR = "#3b7dd8"
BIB_COLOR_AFTER = "#0b2e63"


# Only these per-hit fields are needed downstream of add_time_of_flight
# (for the cuts themselves, and for region_table/add_peak_density). This
# script is the first to combine add_incidence_angles + add_time_of_flight
# (many intermediate per-hit arrays) with region_table/add_peak_density
# (which iterates per-region over the full hit arrays) in one process, on
# the full 16M-hit BIB sample - together that's enough memory to get
# OOM-killed on a laptop, so intermediate fields no longer needed (the
# sensor-normal/incidence-angle components, TOF intermediates, raw
# momentum, etc.) are dropped as soon as the fields we actually need
# (t_corrected_ns, z_axis_intercept_mm, inv_radius_per_mm) are computed.
FIELDS_NEEDED_FOR_CUTS_AND_DENSITY = {
    "system", "side", "layer", "module", "x", "y", "z", "r", "u", "v",
    "t_corrected_ns", "z_axis_intercept_mm", "inv_radius_per_mm",
}


def slim_hits(hits):
    """Drop per-hit fields not needed for cuts/density (see
    FIELDS_NEEDED_FOR_CUTS_AND_DENSITY), freeing their memory immediately.
    Mutates hits in place and returns it."""
    for k in list(hits.keys()):
        if k.startswith("_"):
            continue
        if k not in FIELDS_NEEDED_FOR_CUTS_AND_DENSITY:
            del hits[k]
    return hits


def density_table(hits, area_lookup, with_peak):
    rows = region_table(hits)
    if area_lookup is not None:
        geom_mod.annotate_rows_with_geometry_area(rows, area_lookup)
    if with_peak:
        add_peak_density(hits, rows, bin_size_mm=None, percentile=PEAK_PERCENTILE,
                          target_hits_per_bin=PEAK_TARGET_HITS_PER_BIN)
    return subsystem_density_table(rows)


def cutflow_rows(hits, combined_mask, per_cut_masks):
    system = hits["system"]
    sys_ids = sorted(SYSTEM_NAMES.keys())
    rows = []
    for s in sys_ids:
        smask = system == s
        n_total = int(smask.sum())
        row = {"system": s, "system_name": SYSTEM_NAMES[s], "n_total": n_total}
        for name in CUT_ORDER:
            row[f"n_after_{name}"] = int((smask & per_cut_masks[name]).sum())
        row["n_after_all"] = int((smask & combined_mask).sum())
        row["frac_after_all"] = row["n_after_all"] / n_total if n_total else float("nan")
        rows.append(row)
    n_total_all = int(len(system))
    overall = {"system": "ALL", "system_name": "ALL", "n_total": n_total_all}
    for name in CUT_ORDER:
        overall[f"n_after_{name}"] = int(per_cut_masks[name].sum())
    overall["n_after_all"] = int(combined_mask.sum())
    overall["frac_after_all"] = overall["n_after_all"] / n_total_all if n_total_all else float("nan")
    rows.append(overall)
    return rows


def write_cutflow_csv(path, rows):
    fields = ["system", "system_name", "n_total"] + \
             [f"n_after_{n}" for n in CUT_ORDER] + ["n_after_all", "frac_after_all"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def print_cutflow(label, rows):
    print(f"\n-- {label}: cutflow (n_hits after each individual cut, and after all combined) --")
    header = f"  {'region':12s}  {'n_total':>10s}  " + \
             "  ".join(f"{n.split('_')[0]:>10s}" for n in CUT_ORDER) + \
             f"  {'combined':>10s}  {'frac_pass':>9s}"
    print(header)
    for r in rows:
        print(f"  {r['system_name']:12s}  {r['n_total']:>10,}  " +
              "  ".join(f"{r[f'n_after_{n}']:>10,}" for n in CUT_ORDER) +
              f"  {r['n_after_all']:>10,}  {r['frac_after_all']*100:8.3f}%")


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <bib.root> <signal.root> "
              f"[cuts_config] [output_dir] [geometry_dir]")
        sys.exit(1)
    bib_file = sys.argv[1]
    signal_file = sys.argv[2]
    cuts_config = sys.argv[3] if len(sys.argv) > 3 else \
        str(Path(__file__).resolve().parent / "cuts_config.txt")
    outdir = Path(sys.argv[4] if len(sys.argv) > 4 else "../output_cuts")
    geom_dir = Path(sys.argv[5]) if len(sys.argv) > 5 else Path(bib_file).resolve().parent
    outdir = prepare_output_dir(outdir)

    import gc

    print(f"Loading BIB file: {_display_path(bib_file)}")
    bib_hits = load_hits(bib_file)
    add_incidence_angles(bib_hits)
    add_time_of_flight(bib_hits)
    n_bib_hit = len(bib_hits["x"])
    slim_hits(bib_hits)
    gc.collect()
    print(f"  {bib_hits['_tree_name']}: {bib_hits['_n_events']} event(s), n_hit={n_bib_hit:,}")

    print(f"Loading signal file: {_display_path(signal_file)}")
    sig_hits = load_hits(signal_file)
    add_incidence_angles(sig_hits)
    add_time_of_flight(sig_hits)
    n_sig_events = sig_hits["_n_events"]
    n_sig_hit = len(sig_hits["x"])
    slim_hits(sig_hits)
    gc.collect()
    print(f"  {sig_hits['_tree_name']}: {n_sig_events:,} event(s), n_hit={n_sig_hit:,}")

    print(f"\nLoading cuts: {_display_path(cuts_config)}")
    cuts = load_cuts(cuts_config)
    for name in CUT_ORDER:
        c = cuts[name]
        state = "enabled" if c["enabled"] else "DISABLED"
        print(f"  {name:20s}  center={c['center']:+.4g}  halfwidth={c['halfwidth']:.4g}  [{state}]")

    bib_mask, bib_per_cut = apply_cuts(bib_hits, cuts)
    sig_mask, sig_per_cut = apply_cuts(sig_hits, cuts)

    bib_cutflow = cutflow_rows(bib_hits, bib_mask, bib_per_cut)
    sig_cutflow = cutflow_rows(sig_hits, sig_mask, sig_per_cut)
    write_cutflow_csv(outdir / "cutflow_bib.csv", bib_cutflow)
    write_cutflow_csv(outdir / "cutflow_signal.csv", sig_cutflow)
    print_cutflow("BIB", bib_cutflow)
    print_cutflow("Signal", sig_cutflow)

    # The cut fields themselves aren't needed for the density tables below
    # (only system/side/layer/module/x/y/z/r/u/v are) - drop them now to
    # free more memory before the (heaviest) full-16M-hit BIB peak-density
    # computation.
    for hits in (bib_hits, sig_hits):
        for k in ("t_corrected_ns", "z_axis_intercept_mm", "inv_radius_per_mm"):
            hits.pop(k, None)
    gc.collect()

    # ---- geometry-based area, if available -------------------------------
    area_lookup = None
    if all((geom_dir / f).exists() for f in GEOMETRY_FILES):
        try:
            area_lookup = geom_mod.build_area_lookup(geom_dir)
            print(f"\nUsing true geometry-based area from {_display_path(geom_dir)}.")
        except Exception as e:
            print(f"\nWARNING: geometry parsing failed ({e}); using hit-inferred area.")
    else:
        print(f"\nNOTE: geometry files not found in {_display_path(geom_dir)}; "
              f"using hit-inferred area.")

    # ---- BIB density before/after cuts -------------------------------------
    bib_before = density_table(bib_hits, area_lookup, with_peak=True)
    bib_after = density_table(mask_hits(bib_hits, bib_mask), area_lookup, with_peak=True)

    bib_before_by_sys = {r["system"]: r for r in bib_before}
    bib_after_by_sys = {r["system"]: r for r in bib_after}
    sys_ids = sorted(SYSTEM_NAMES.keys())

    csv_path = outdir / "density_before_after_cuts.csv"
    fields = ["system", "system_name",
              "bib_n_hits_before", "bib_mean_density_before", f"bib_{PEAK_KEY}_before",
              "bib_n_hits_after", "bib_mean_density_after", f"bib_{PEAK_KEY}_after",
              "bib_rejection_frac"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, restval="")
        w.writeheader()
        for s in sys_ids:
            bb = bib_before_by_sys.get(s, {"n_hits": 0, "density_hits_per_mm2": 0.0})
            ba = bib_after_by_sys.get(s, {"n_hits": 0, "density_hits_per_mm2": 0.0})
            bib_n_before = bb["n_hits"]
            w.writerow({
                "system": s, "system_name": SYSTEM_NAMES[s],
                "bib_n_hits_before": bib_n_before,
                "bib_mean_density_before": bb["density_hits_per_mm2"],
                f"bib_{PEAK_KEY}_before": bb.get(PEAK_KEY, ""),
                "bib_n_hits_after": ba["n_hits"],
                "bib_mean_density_after": ba["density_hits_per_mm2"],
                f"bib_{PEAK_KEY}_after": ba.get(PEAK_KEY, ""),
                "bib_rejection_frac": 1.0 - ba["n_hits"] / bib_n_before if bib_n_before else "",
            })

    # ---- plot: BIB density before/after ------------------------------------
    labels = [SYSTEM_NAMES[s] for s in sys_ids]
    x = np.arange(len(sys_ids))
    width = 0.35

    fig, ax1 = plt.subplots(figsize=(11, 5))

    bib_mean_before = [bib_before_by_sys.get(s, {}).get("density_hits_per_mm2", 0.0) for s in sys_ids]
    bib_mean_after = [bib_after_by_sys.get(s, {}).get("density_hits_per_mm2", 0.0) for s in sys_ids]
    b1 = ax1.bar(x - width / 2, bib_mean_before, width, color=BIB_COLOR, label="BIB, no cuts")
    b2 = ax1.bar(x + width / 2, bib_mean_after, width, color=BIB_COLOR_AFTER, label="BIB, with cuts")
    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=20)
    ax1.set_ylabel("mean hit density (hits/mm$^2$)")
    ax1.set_title("BIB hit density: with vs. without cuts")
    ax1.legend()
    for bars in (b1, b2):
        for b in bars:
            ax1.text(b.get_x() + b.get_width() / 2, b.get_height(),
                     f"{b.get_height():.3g}", ha="center", va="bottom", fontsize=7)

    plt.tight_layout()
    plt.savefig(outdir / "density_before_after_cuts.png", dpi=140)
    plt.close(fig)

    print(f"\nWrote cutflow_bib.csv, cutflow_signal.csv, density_before_after_cuts.csv, "
          f"and density_before_after_cuts.png to {_display_path(outdir.resolve())}")

    print("\n-- Summary: BIB rejection vs. signal efficiency (all cuts combined) --")
    for r_b, r_s in zip(bib_cutflow, sig_cutflow):
        rej = 1.0 - r_b["frac_after_all"] if r_b["n_total"] else float("nan")
        eff = r_s["frac_after_all"]
        print(f"  {r_b['system_name']:12s}  BIB rejection={rej*100:6.2f}%   "
              f"signal efficiency={eff*100:6.2f}%")


if __name__ == "__main__":
    main()
