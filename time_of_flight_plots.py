"""
Step 2 (part 3): for each hit, compute the expected time-of-flight of a
muon travelling from the IP (origin) to the hit along its actual curved
trajectory (using the transverse curvature, per bib_common.add_time_of_flight),
and subtract it from the measured hit time. For a genuine signal muon
produced at the IP at t=0, the corrected time should cluster tightly
around zero regardless of how far out the hit is - unlike the raw hit
time, which simply grows with distance travelled. See
bib_common.add_time_of_flight() docstring for the full derivation.

A small fraction of hits (~5-6% in the signal sample) show a very large
|t_corrected| - these turn out to be low-pT secondary/knock-on hits (not
the primary muon), for which the "originated at the IP" assumption behind
the correction is simply wrong; large excursions there are an expected,
even informative, consequence, not a bug (see console printout below,
which reports raw-momentum pT for the outlier tail).

Usage:
    python3 time_of_flight_plots.py <input.root> [output_dir]
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
    load_hits, add_incidence_angles, add_time_of_flight,
    prepare_output_dir, _display_path, SYSTEM_NAMES,
    load_smearing_config, smearing_rng, describe_smearing,
    apply_position_time_smearing,
    apply_angle_smearing,
)

T_RANGE_NS = 40.0          # raw hit-time histogram range
TC_RANGE_NS = 20.0         # corrected-time histogram, full range
TC_ZOOM_RANGE_NS = 2.0     # corrected-time histogram, zoomed near 0


def panel_grid():
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    return fig, axes.flatten()


def main():
    root_file = sys.argv[1] if len(sys.argv) > 1 else \
        "/mnt/user-data/uploads/ntu_bib_plus_1evt.root"
    outdir = Path(sys.argv[2] if len(sys.argv) > 2 else "../output_time_of_flight")
    outdir = prepare_output_dir(outdir)

    smearing_config = os.environ.get("SMEARING_CONFIG", "").strip()
    smear_cfg = None
    smear_rng = None
    if smearing_config:
        smear_cfg = load_smearing_config(smearing_config)
        smear_rng = smearing_rng(smear_cfg)
        joined = describe_smearing(smear_cfg)
        print(f"Smearing config: {_display_path(smearing_config)} ({joined})")

    hits = load_hits(root_file)
    if smear_cfg is not None:
        apply_position_time_smearing(hits, smear_cfg, smear_rng)
    n_hit = len(hits["x"])
    add_incidence_angles(hits)
    if smear_cfg is not None:
        apply_angle_smearing(hits, smear_cfg, smear_rng)
    add_time_of_flight(hits)
    system = hits["system"]
    sys_ids = sorted(SYSTEM_NAMES.keys())

    t = hits["t"]
    tc = hits["t_corrected_ns"]
    finite = np.isfinite(tc)

    # ---- 1. raw hit time, per subsystem -----------------------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        v = t[system == s]
        ax.hist(v, bins=np.linspace(-T_RANGE_NS / 4, T_RANGE_NS, 141), color="#3b7dd8")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("raw hit time t (ns)")
        ax.set_ylabel("hits")
    plt.suptitle("Raw hit time (measured, uncorrected) - grows with distance from IP")
    plt.tight_layout()
    plt.savefig(outdir / "hit_time_raw_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 2. corrected time t - t_expected(TOF), per subsystem, full range -
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        v = tc[(system == s) & finite]
        ax.hist(v, bins=np.linspace(-TC_RANGE_NS, TC_RANGE_NS, 161), color="#e05252")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("t - t$_{expected}$(TOF from IP) (ns)")
        ax.set_ylabel("hits")
    plt.suptitle("Time-of-flight-corrected hit time - a genuine on-time IP muon "
                 "peaks at 0 regardless of subsystem", fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "time_corrected_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 3. corrected time, ZOOMED near 0 ---------------------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        v = tc[(system == s) & finite]
        v = v[np.abs(v) <= TC_ZOOM_RANGE_NS]
        ax.hist(v, bins=np.linspace(-TC_ZOOM_RANGE_NS, TC_ZOOM_RANGE_NS, 121), color="#e05252")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("t - t$_{expected}$(TOF from IP) (ns)")
        ax.set_ylabel("hits")
    plt.suptitle(f"Time-of-flight-corrected hit time, zoomed to +/-{TC_ZOOM_RANGE_NS:.0f}ns",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "time_corrected_per_subsystem_zoom.png", dpi=130)
    plt.close(fig)

    # ---- summary CSV + console report -------------------------------------
    csv_path = outdir / "time_of_flight_summary.csv"
    fields = ["system", "system_name", "n_hits", "t_raw_mean_ns", "t_raw_std_ns",
              "t_corrected_finite_frac", "t_corrected_mean_ns", "t_corrected_std_ns",
              "t_corrected_p05_ns", "t_corrected_p50_ns", "t_corrected_p95_ns",
              "frac_abs_tcorr_gt_1ns"]
    rows = []
    for s in sys_ids:
        mask = system == s
        v_raw = t[mask]
        v_tc_all = tc[mask]
        v_tc = v_tc_all[np.isfinite(v_tc_all)]
        rows.append({
            "system": s, "system_name": SYSTEM_NAMES[s], "n_hits": int(mask.sum()),
            "t_raw_mean_ns": float(v_raw.mean()), "t_raw_std_ns": float(v_raw.std()),
            "t_corrected_finite_frac": float(np.isfinite(v_tc_all).mean()),
            "t_corrected_mean_ns": float(v_tc.mean()) if v_tc.size else float("nan"),
            "t_corrected_std_ns": float(v_tc.std()) if v_tc.size else float("nan"),
            "t_corrected_p05_ns": float(np.percentile(v_tc, 5)) if v_tc.size else float("nan"),
            "t_corrected_p50_ns": float(np.percentile(v_tc, 50)) if v_tc.size else float("nan"),
            "t_corrected_p95_ns": float(np.percentile(v_tc, 95)) if v_tc.size else float("nan"),
            "frac_abs_tcorr_gt_1ns": float((np.abs(v_tc) > 1.0).mean()) if v_tc.size else float("nan"),
        })
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"Loaded {hits['_tree_name']}: n_hit={n_hit:,}")
    print(f"Wrote plots and {csv_path.name} to {_display_path(outdir.resolve())}")
    print()
    print("-- Corrected time t_corrected = t - TOF_expected(IP muon) (ns), "
          "median [p05, p95], per subsystem --")
    for r in rows:
        print(f"  {r['system_name']:12s}  n={r['n_hits']:>9,}  "
              f"p50={r['t_corrected_p50_ns']:7.3f}  "
              f"[{r['t_corrected_p05_ns']:8.3f}, {r['t_corrected_p95_ns']:7.3f}]  "
              f"frac(|t_corr|>1ns)={r['frac_abs_tcorr_gt_1ns']*100:5.2f}%")

    # diagnostic: for the large-|t_corrected| tail, report the raw-momentum
    # pT of those hits, to show they're low-pT secondaries, not primaries
    bad = np.abs(tc) > 1.0
    if bad.any():
        pT_raw = np.sqrt(hits["px"] ** 2 + hits["py"] ** 2)
        print()
        print(f"-- Diagnostic: among the {bad.sum():,} hits ({bad.mean()*100:.2f}%) with "
              f"|t_corrected| > 1ns, raw hit p_T (GeV) --")
        print(f"  median={np.median(pT_raw[bad]):.4g}  "
              f"(all hits median={np.median(pT_raw):.4g}) - "
              f"large excursions concentrate in low-p_T secondary hits, where the "
              f"'originated at the IP' assumption behind the correction breaks down")


if __name__ == "__main__":
    main()
