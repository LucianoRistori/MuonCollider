"""
Step 4: apply the selection cuts defined in __cuts_config.txt (time-of-
flight-corrected hit time, z-axis intercept, curvature-based momentum -
see bib_common.load_cuts / bib_common.apply_cuts) to both the BIB and
signal samples, and compare - the main point of this script - hit
density and hit counts WITH cuts vs. WITHOUT cuts, per subsystem, so a
cut hypothesis can be judged by how much it suppresses BIB relative to
how much signal it keeps.

The cuts are set per subsystem (__cuts_config.txt): every hit is judged
by its own subsystem's cuts.

Two kinds of comparison are produced:
  1. A cutflow (hit counts before / after each individual cut / after
     all three combined, per subsystem) for BIB and for signal
     separately - this shows each cut's own selectivity.
     The run's summary table is built from it: per subsystem, the BIB
     rejection factor (hits before / hits after = 1/(1 - R), R being
     the fraction of BIB hits removed) and the signal hit efficiency in
     the limit pT -> infinity (bib_common.efficiency_at_infinite_pt -
     not averaged over the sample, which is flat in 1/pT from 1.5
     GeV/c and so mostly made of muons below the pT cut). Also written
     to summary_by_region.csv.
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
import os
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
    load_smearing_config, smearing_rng, describe_smearing,
    under_run_all, short_path, loaded_line, print_table, resolve_geometry, default_cuts_config,
    apply_position_time_smearing,
    apply_angle_smearing,
    rejection_factor, format_rejection_factor, pt_inf_fit_min, efficiency_at_infinite_pt,
    PT_INF_FIT_MIN_GEV,
)
import cuts_table
import geometry as geom_mod

PEAK_PERCENTILE = 99.0
PEAK_KEY = f"peak_density_p{PEAK_PERCENTILE:.0f}_hits_per_mm2"
PEAK_TARGET_HITS_PER_BIN = 20.0
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


def slim_hits(hits, keep=()):
    """Drop per-hit fields not needed for cuts/density (see
    FIELDS_NEEDED_FOR_CUTS_AND_DENSITY) or listed in `keep`, freeing their
    memory immediately. Mutates hits in place and returns it."""
    for k in list(hits.keys()):
        if k.startswith("_"):
            continue
        if k not in FIELDS_NEEDED_FOR_CUTS_AND_DENSITY and k not in keep:
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
    print()
    print_table(
        f"{label} cutflow - hits passing each cut on its own, and all cuts combined:",
        ["region", "total"] + [n.split("_")[0] for n in CUT_ORDER] + ["combined", "passing"],
        [[r["system_name"], f"{r['n_total']:,}"]
         + [f"{r[f'n_after_{n}']:,}" for n in CUT_ORDER]
         + [f"{r['n_after_all']:,}", f"{r['frac_after_all']*100:.3f}%"] for r in rows])


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <bib.root> <signal.root> "
              f"[cuts_config] [output_dir] [geometry_dir]")
        sys.exit(1)
    bib_file = sys.argv[1]
    signal_file = sys.argv[2]
    cuts_config = sys.argv[3] if len(sys.argv) > 3 else \
        default_cuts_config()
    outdir = Path(sys.argv[4] if len(sys.argv) > 4 else "../output_cuts")
    geom_dir = Path(sys.argv[5]) if len(sys.argv) > 5 else Path(bib_file).resolve().parent
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

    import gc
    bib_hits = load_hits(bib_file)
    if smear_cfg is not None:
        apply_position_time_smearing(bib_hits, smear_cfg, smear_rng)
    add_incidence_angles(bib_hits)
    if smear_cfg is not None:
        apply_angle_smearing(bib_hits, smear_cfg, smear_rng)
    add_time_of_flight(bib_hits)
    n_bib_hit = len(bib_hits["x"])
    slim_hits(bib_hits)
    gc.collect()
    print(loaded_line(bib_file, bib_hits, "BIB"))
    sig_hits = load_hits(signal_file)
    if smear_cfg is not None:
        apply_position_time_smearing(sig_hits, smear_cfg, smear_rng)
    add_incidence_angles(sig_hits)
    if smear_cfg is not None:
        apply_angle_smearing(sig_hits, smear_cfg, smear_rng)
    add_time_of_flight(sig_hits)
    n_sig_events = sig_hits["_n_events"]
    n_sig_hit = len(sig_hits["x"])
    slim_hits(sig_hits, keep={"event_id"})   # event_id: generated pT of each hit
    gc.collect()
    print(loaded_line(signal_file, sig_hits, "signal"))

    cuts = load_cuts(cuts_config)
    if not under_run_all():
        print(f"Cuts: {short_path(cuts_config)}")

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
    geometry = resolve_geometry(geom_dir)
    if all(p.exists() for p in geometry.values()):
        try:
            area_lookup = geom_mod.build_area_lookup(geometry)
            print("Sensitive areas from the detector geometry")
        except Exception as e:
            print(f"WARNING: geometry parsing failed ({e}); using hit-inferred area.")
    else:
        print(f"NOTE: geometry files not found in {short_path(geom_dir)}; "
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
              "bib_rejection_factor"]
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
                "bib_rejection_factor": rejection_factor(bib_n_before, ba["n_hits"])
                if bib_n_before else "",
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

    print(f"Wrote cutflow_bib.csv, cutflow_signal.csv, density_before_after_cuts.csv "
          f"and .png, summary_by_region.csv to {short_path(outdir)}/")

    # Signal hit efficiency for pT -> infinity, per subsystem and overall:
    # each hit gets the generated pT of its event's muon.
    px, py = sig_hits["_primary"]["px"], sig_hits["_primary"]["py"]
    sig_pt = np.sqrt(px ** 2 + py ** 2)[sig_hits["event_id"]]
    sig_sys = sig_hits["system"]
    eff_inf = {}
    for s in sorted(SYSTEM_NAMES.keys()) + ["ALL"]:
        sel = (sig_sys == s) if s != "ALL" else np.ones(len(sig_sys), dtype=bool)
        pt_min = pt_inf_fit_min(cuts, None if s == "ALL" else [s])
        eff_inf[s] = efficiency_at_infinite_pt(sig_pt[sel], sig_mask[sel], pt_min,
                                               groups=sig_hits["event_id"][sel]) + (pt_min,)

    # When the cuts differ between subsystems, each region's own cut values
    # are shown next to its results.
    per_system = any(len(set(c["per_system"].values())) > 1 for c in cuts.values())
    summary_rows, csv_rows = [], []
    for r_b, r_s in zip(bib_cutflow, sig_cutflow):
        s = r_b["system"]
        eff, unc, n_fit, pt_min = eff_inf[s]
        cut_cells = []
        if per_system:
            cut_cells = ([cuts_table.fmt(cuts[n]["per_system"][s]) for n in CUT_ORDER]
                         if s in SYSTEM_NAMES else ["", "", ""])
        summary_rows.append([r_b["system_name"]] + cut_cells
                            + [format_rejection_factor(r_b["n_total"], r_b["n_after_all"]),
                               f"{eff*100:.1f} +- {unc*100:.1f}%"])
        csv_rows.append({
            "system": s, "system_name": r_b["system_name"],
            "time_cut_ns": cuts_table.fmt(cuts["t_corrected_ns"]["per_system"].get(s)) if s in SYSTEM_NAMES else "",
            "z0_cut_mm": cuts_table.fmt(cuts["z_axis_intercept_mm"]["per_system"].get(s)) if s in SYSTEM_NAMES else "",
            "pt_cut_gev": cuts_table.fmt(cuts["momentum_gev"]["per_system"].get(s)) if s in SYSTEM_NAMES else "",
            "bib_n_hits": r_b["n_total"], "bib_n_hits_after_cuts": r_b["n_after_all"],
            "bib_rejection_factor": rejection_factor(r_b["n_total"], r_b["n_after_all"]),
            "signal_efficiency_pt_inf": eff, "signal_efficiency_pt_inf_unc": unc,
            "fit_pt_min_gev": pt_min, "fit_n_hits": n_fit,
        })
    with open(outdir / "summary_by_region.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
        w.writeheader()
        w.writerows(csv_rows)
    header = (["region"] + (["time (ns)", "z0 (mm)", "pT (GeV/c)"] if per_system else [])
              + ["BIB rejection factor", "signal efficiency, pT -> inf"])
    print()
    print_table("Summary: BIB rejection factor and signal hit efficiency for pT -> inf "
                "(all cuts combined)", header, summary_rows)
    pt_mins = sorted({v[3] for v in eff_inf.values()})
    window = (f"pT > {pt_mins[0]:g} GeV/c" if len(pt_mins) == 1
              else "pT above twice their subsystem's pT cut (at least "
                   f"{PT_INF_FIT_MIN_GEV:g} GeV/c)")
    print("  Rejection factor: BIB hits before / after the cuts = 1/(1 - R), "
          "R = fraction removed.")
    print(f"  Efficiency for pT -> inf: fit eff + c/pT^2 to the hits of muons with {window}.")
    print()


if __name__ == "__main__":
    main()
