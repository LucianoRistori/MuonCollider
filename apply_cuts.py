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
     to summary_by_region.csv. Two more tables break this down by cut:
     each cut alone, and each cut on top of the other two (as in the N-1
     plots: BIB hits passing the other two / passing all three, and the
     signal efficiency among the hits passing the other two - dropping
     that cut divides the combined numbers by these).
  2. A before/after hit-density comparison (mean density, hits/mm^2,
     using the true detector geometry area when available - same
     definition as step 1/3) for BIB only: per subsystem
     (density_before_after_cuts.csv / .png), and per layer or disk
     (density_per_layer_before_after_cuts.csv and a table in the log:
     mean and peak density before and after the cuts, and the layer's
     rejection factor). (Signal hit density isn't
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
    signal_muon_hits_only,
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
    """Returns (one row per layer/disk, one row per subsystem)."""
    rows = region_table(hits)
    if area_lookup is not None:
        geom_mod.annotate_rows_with_geometry_area(rows, area_lookup)
    if with_peak:
        add_peak_density(hits, rows, bin_size_mm=None, percentile=PEAK_PERCENTILE,
                          target_hits_per_bin=PEAK_TARGET_HITS_PER_BIN)
    return rows, subsystem_density_table(rows)


SIDE_TAG = {0: "", 1: " +z", 3: " -z"}
LAYER_CSV_FIELDS = [
    "system", "system_name", "side", "layer", "label", "mean_r_mm", "mean_z_mm", "area_mm2",
    "bib_n_hits_before", "bib_mean_density_before", f"bib_{PEAK_KEY}_before",
    "bib_n_hits_after", "bib_mean_density_after", f"bib_{PEAK_KEY}_after",
    "bib_rejection_factor"]


def density_per_layer(before_rows, after_rows):
    """
    One row per layer (barrel) or disk (endcap: +z and -z listed
    separately, -z first), in the order barrel layers by radius, disks by
    |z|: BIB hits, mean density (hits/mm^2 over the layer's sensitive
    area, from the geometry) and peak density (PEAK_PERCENTILE-th
    percentile over bins of ~PEAK_TARGET_HITS_PER_BIN hits, as in step 1)
    before and after all cuts. The density after the cuts uses the same
    area as before, so before/after is exactly the layer's rejection
    factor.
    """
    after = {(r["system"], r["side"], r["layer"]): r for r in after_rows}
    rows = []
    for b in sorted(before_rows, key=lambda r: (r["system"], r["layer"], r["side"] != 3)):
        a = after.get((b["system"], b["side"], b["layer"]))
        n_after = a["n_hits"] if a else 0
        area = b["area_mm2"]
        rows.append({
            "system": b["system"], "system_name": b["system_name"], "side": b["side"],
            "layer": b["layer"], "label": f"L{b['layer']}{SIDE_TAG.get(b['side'], '')}",
            "mean_r_mm": b["mean_r"], "mean_z_mm": b["mean_z"], "area_mm2": area,
            "bib_n_hits_before": b["n_hits"],
            "bib_mean_density_before": b["density_hits_per_mm2"],
            f"bib_{PEAK_KEY}_before": b.get(PEAK_KEY, float("nan")),
            "bib_n_hits_after": n_after,
            "bib_mean_density_after": n_after / area if area > 0 else float("nan"),
            f"bib_{PEAK_KEY}_after": a.get(PEAK_KEY, float("nan")) if a else float("nan"),
            "bib_rejection_factor": rejection_factor(b["n_hits"], n_after),
        })
    return rows


def print_density_per_layer(rows):
    def g(x):
        return f"{x:.3g}" if x == x else "-"          # x == x: not NaN
    table, last = [], None
    for r in rows:
        barrel = r["side"] == 0
        table.append([r["system_name"] if r["system_name"] != last else "", r["label"],
                      f"r {r['mean_r_mm']:.0f} mm" if barrel else f"z {r['mean_z_mm']:+.0f} mm",
                      g(r["bib_mean_density_before"]), g(r["bib_mean_density_after"]),
                      g(r[f"bib_{PEAK_KEY}_before"]), g(r[f"bib_{PEAK_KEY}_after"]),
                      format_rejection_factor(r["bib_n_hits_before"], r["bib_n_hits_after"])])
        last = r["system_name"]
    print()
    print_table(f"BIB hit density per layer, before and after all cuts (hits/mm^2 per collision; "
                f"peak = p{PEAK_PERCENTILE:.0f} over bins of ~{PEAK_TARGET_HITS_PER_BIN:.0f} hits):",
                ["subsystem", "layer", "position", "mean before", "mean after", "peak before",
                 "peak after", "rejection factor"], table, align="lllrrrrr")


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


CUT_SHORT = {"t_corrected_ns": "time", "z_axis_intercept_mm": "z0", "momentum_gev": "pt"}
CUT_SHORT_TITLE = {"t_corrected_ns": "time", "z_axis_intercept_mm": "z0", "momentum_gev": "pT"}


def cut_is_off(cuts, name, s):
    """True if cut `name` removes nothing in region s (a system id, or "ALL"
    for all of them): off, or - for pT - a cut at 0."""
    systems = SYSTEM_NAMES.keys() if s == "ALL" else [s]
    return all(v is None or (name == "momentum_gev" and v == 0)
               for v in (cuts[name]["per_system"].get(x) for x in systems))


def per_cut_breakdown(cuts, bib_sys, bib_mask, bib_per_cut,
                      sig_sys, sig_mask, sig_per_cut, sig_pt, sig_event):
    """
    For each region (system id, and "ALL") and each cut: what the cut does
    on its own ("alone": BIB hits before / after it, signal efficiency of
    it alone) and on top of the other two ("n1": BIB hits passing the other
    two / passing all three, signal efficiency among the hits passing the
    other two) - see main(). Returns {region: {cut: {"alone": (rejection
    factor as text, rejection factor, eff, unc), "n1": (...)} or None}}.
    """
    regions = sorted(SYSTEM_NAMES.keys()) + ["ALL"]

    def bib_counts(mask):
        c = np.bincount(bib_sys[mask], minlength=32)
        return {**{s: int(c[s]) for s in SYSTEM_NAMES}, "ALL": int(mask.sum())}
    n_total = bib_counts(np.ones(len(bib_sys), dtype=bool))
    n_all = bib_counts(bib_mask)
    results = {s: {} for s in regions}
    for name in CUT_ORDER:
        o1, o2 = [o for o in CUT_ORDER if o != name]
        n_alone = bib_counts(bib_per_cut[name])
        n_others = bib_counts(bib_per_cut[o1] & bib_per_cut[o2])
        sig_others = sig_per_cut[o1] & sig_per_cut[o2]
        for s in regions:
            if cut_is_off(cuts, name, s):
                results[s][name] = None
                continue
            sel = (sig_sys == s) if s != "ALL" else np.ones(len(sig_sys), dtype=bool)
            pt_min = pt_inf_fit_min(cuts, None if s == "ALL" else [s])
            e_alone = efficiency_at_infinite_pt(sig_pt[sel], sig_per_cut[name][sel], pt_min,
                                                groups=sig_event[sel])
            sel_n1 = sel & sig_others
            e_n1 = efficiency_at_infinite_pt(sig_pt[sel_n1], sig_mask[sel_n1], pt_min,
                                             groups=sig_event[sel_n1])
            results[s][name] = {
                "alone": (format_rejection_factor(n_total[s], n_alone[s]),
                          rejection_factor(n_total[s], n_alone[s]), e_alone[0], e_alone[1]),
                "n1": (format_rejection_factor(n_others[s], n_all[s]),
                       rejection_factor(n_others[s], n_all[s]), e_n1[0], e_n1[1]),
            }
    return results


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
    sig_hits = load_hits(signal_file, muon_hits_only=signal_muon_hits_only(cuts_config))
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
    bib_before_layers, bib_before = density_table(bib_hits, area_lookup, with_peak=True)
    bib_after_layers, bib_after = density_table(mask_hits(bib_hits, bib_mask), area_lookup,
                                                with_peak=True)

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

    # ---- BIB density per layer, before/after cuts --------------------------
    layer_rows = density_per_layer(bib_before_layers, bib_after_layers)
    with open(outdir / "density_per_layer_before_after_cuts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=LAYER_CSV_FIELDS)
        w.writeheader()
        w.writerows(layer_rows)

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

    print(f"Wrote cutflow_bib.csv, cutflow_signal.csv, density_before_after_cuts.csv and .png, "
          f"density_per_layer_before_after_cuts.csv, summary_by_region.csv to {short_path(outdir)}/")
    print_density_per_layer(layer_rows)

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

    # Each cut alone, and each cut on top of the other two (as in the N-1
    # plots), per subsystem and overall: BIB rejection factor from the hit
    # counts, signal efficiency for pT -> inf from the same fit as above.
    # per_cut_results[region][cut] = {"alone": (rejection text, rejection,
    # eff, unc), "n1": (...)}, or None where that cut is off.
    per_cut_results = per_cut_breakdown(cuts, bib_hits["system"], bib_mask, bib_per_cut,
                                        sig_sys, sig_mask, sig_per_cut, sig_pt,
                                        sig_hits["event_id"])

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
        for name in CUT_ORDER:
            for key in ("alone", "n1"):
                res = per_cut_results[s][name]
                tag = f"{CUT_SHORT[name]}_{key}"
                csv_rows[-1].update({
                    f"bib_rejection_factor_{tag}": res[key][1] if res else "",
                    f"signal_efficiency_pt_inf_{tag}": res[key][2] if res else "",
                    f"signal_efficiency_pt_inf_{tag}_unc": res[key][3] if res else "",
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
    whose = ("the muon's own hits" if "_n_hits_all" in sig_hits
             else "all hits (secondaries included)")
    print(f"  Efficiency for pT -> inf: fit eff + c/pT^2 to {whose}, muons with {window}.")

    header = ["region"] + sum(([f"{CUT_SHORT_TITLE[n]}: rej.", "eff."] for n in CUT_ORDER), [])
    for key, title in (("alone", "Each cut alone - BIB rejection factor and signal efficiency "
                                 "for pT -> inf:"),
                       ("n1", "Each cut on top of the other two, as in the N-1 plots:")):
        rows = []
        for r_b in bib_cutflow:
            s = r_b["system"]
            row = [r_b["system_name"]]
            for name in CUT_ORDER:
                res = per_cut_results[s][name]
                row += ([res[key][0], f"{res[key][2]*100:.1f} +- {res[key][3]*100:.1f}%"] if res
                        else ["off", ""])
            rows.append(row)
        print()
        print_table(title, header, rows)
    print("  Dropping a cut divides the combined rejection factor and efficiency by its "
          "values in this last table.")
    print()


if __name__ == "__main__":
    main()
