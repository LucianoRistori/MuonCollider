"""
Step 3 (part 2): overlay the muon-gun signal sample on the incidence-
angle/curvature and time-of-flight diagnostic plots from
incidence_angle_plots.py / time_of_flight_plots.py that aren't already
covered by signal_overlay_plots.py:
  - z_axis_intercept_per_subsystem      (full range, +/-3000mm)
  - z_axis_intercept_per_subsystem_zoom (+/-100mm)
  - inv_radius_per_subsystem            (full range, relabeled in pT GeV/c)
  - inv_radius_per_subsystem_zoom       (zoomed to |pT| ~ 1 GeV/c window)
  - time_corrected_per_subsystem        (full +/-20ns range)
  - time_corrected_per_subsystem_zoom   (+/-2ns zoom)

Normalization (hits per collision): the merged BIB file (plus+minus) is
taken to represent the full background expected in a single MC collision,
and each signal event represents one collision's worth of the single
muon-gun muon superimposed on it. So both histograms here are put on the
common, addable footing "hits / collision / bin":
  - BIB: raw hit counts per bin (unweighted) - the merged file already IS
    one collision's expected background.
  - signal: hit counts per bin, weighted by 1/n_sig_events - the average
    hits per bin contributed by one signal muon (one collision).
This is different from signal_overlay_plots.py's density plots, which
report BIB density in hits/mm^2 (for the whole sample) side by side with
signal density per event, per mm^2 - both are about spatial density, not
about "how would BIB and signal add up in one collision" the way these
angle/curvature/time histograms now are.

Usage:
    python3 signal_overlay_angle_plots.py <bib.root> <signal.root> [output_dir]
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bib_common import (
    load_hits, add_incidence_angles, add_time_of_flight,
    prepare_output_dir, _display_path, SYSTEM_NAMES, load_cuts,
    load_smearing_config, smearing_rng, apply_position_time_smearing,
    apply_angle_smearing,
)
from incidence_angle_plots import (
    pt_ticks_for_axis, panel_grid,
    Z0_RANGE_MM, Z0_ZOOM_RANGE_MM, PT_ZOOM_RANGE_GEV, B_FIELD_T, GEV_PER_INV_M,
)
from time_of_flight_plots import TC_RANGE_NS, TC_ZOOM_RANGE_NS

BIB_COLOR = "#3b7dd8"
SIG_COLOR = "#2ca858"

CUT_LINE_COLOR = "#a83232"


def draw_symmetric_cut_lines(ax, cuts, name, xmax=None):
    """For a simple |value - center| <= halfwidth cut (t_corrected_ns,
    z_axis_intercept_mm): draw vertical lines at the two edges, if the
    cut is enabled and (optionally) within the visible range."""
    spec = cuts[name]
    if not spec["enabled"]:
        return
    lo, hi = spec["center"] - spec["halfwidth"], spec["center"] + spec["halfwidth"]
    for edge in (lo, hi):
        if xmax is None or abs(edge) <= xmax:
            ax.axvline(edge, color=CUT_LINE_COLOR, linewidth=1.1,
                       linestyle="--", alpha=0.8)


def draw_momentum_cut_lines(ax, cuts, xmax=None):
    """momentum_gev cut accepts |inv_radius_per_mm| <= threshold (a band
    around 1/R=0); draw its two edges on the curvature axis (in 1/m,
    matching the plotted units bv = inv_radius_per_mm*1000)."""
    spec = cuts["momentum_gev"]
    if not spec["enabled"]:
        return
    inv_radius_threshold_per_m = (GEV_PER_INV_M / spec["halfwidth"])
    for edge in (-inv_radius_threshold_per_m, inv_radius_threshold_per_m):
        if xmax is None or abs(edge) <= xmax:
            ax.axvline(edge, color=CUT_LINE_COLOR, linewidth=1.1,
                       linestyle="--", alpha=0.8)



def overlay_hist(ax, bib_v, sig_v, bins, n_sig_events):
    """BIB drawn as raw hit counts per bin (the merged plus+minus file is
    taken to be one collision's full background); signal drawn as hit
    counts per bin divided by the number of signal events, i.e. the
    average contribution of the one signal muon per collision. Both are
    then in the same "hits / collision / bin" units - see module
    docstring."""
    ax.hist(bib_v, bins=bins, histtype="step", color=BIB_COLOR,
            linewidth=1.4, label="BIB (per collision)")
    sig_weights = np.full(len(sig_v), 1.0 / n_sig_events)
    ax.hist(sig_v, bins=bins, weights=sig_weights, histtype="step", color=SIG_COLOR,
            linewidth=1.4, label="signal (per collision, avg. over events)")


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <bib.root> <signal.root> "
              f"[output_dir] [cuts_config]")
        sys.exit(1)
    bib_file = sys.argv[1]
    signal_file = sys.argv[2]
    outdir = Path(sys.argv[3] if len(sys.argv) > 3 else "../output_signal_angle_overlay")
    cuts_config = sys.argv[4] if len(sys.argv) > 4 else \
        str(Path(__file__).resolve().parent / "cuts_config.txt")
    outdir = prepare_output_dir(outdir)

    cuts = load_cuts(cuts_config)
    print(f"Loading cuts (for reference lines): {_display_path(cuts_config)}")
    for name, c in cuts.items():
        state = "enabled" if c["enabled"] else "disabled"
        print(f"  {name:20s}  center={c['center']:+.4g}  halfwidth={c['halfwidth']:.4g}  [{state}]")

    smearing_config = os.environ.get("SMEARING_CONFIG", "").strip()
    smear_cfg = None
    smear_rng = None
    if smearing_config:
        smear_cfg = load_smearing_config(smearing_config)
        smear_rng = smearing_rng(smear_cfg)
        applied = [f"{name} sigma={smear_cfg[name]['sigma']}"
                   for name in ("position", "time", "angle_long", "angle_trans")
                   if smear_cfg[name]["enabled"] and smear_cfg[name]["sigma"] > 0]
        joined = ", ".join(applied) if applied else "all disabled/zero"
        print(f"Smearing config: {_display_path(smearing_config)} ({joined})")

    print(f"Loading BIB file: {_display_path(bib_file)}")
    bib_hits = load_hits(bib_file)
    if smear_cfg is not None:
        apply_position_time_smearing(bib_hits, smear_cfg, smear_rng)
    add_incidence_angles(bib_hits)
    if smear_cfg is not None:
        apply_angle_smearing(bib_hits, smear_cfg, smear_rng)
    add_time_of_flight(bib_hits)
    n_bib_hit = len(bib_hits["x"])
    print(f"  {bib_hits['_tree_name']}: {bib_hits['_n_events']} event(s), n_hit={n_bib_hit:,}")

    print(f"Loading signal file: {_display_path(signal_file)}")
    sig_hits = load_hits(signal_file)
    if smear_cfg is not None:
        apply_position_time_smearing(sig_hits, smear_cfg, smear_rng)
    add_incidence_angles(sig_hits)
    if smear_cfg is not None:
        apply_angle_smearing(sig_hits, smear_cfg, smear_rng)
    add_time_of_flight(sig_hits)
    n_sig_hit = len(sig_hits["x"])
    n_sig_events = sig_hits["_n_events"]
    print(f"  {sig_hits['_tree_name']}: {n_sig_events:,} event(s), n_hit={n_sig_hit:,}")

    bib_sys = bib_hits["system"]
    sig_sys = sig_hits["system"]
    sys_ids = sorted(SYSTEM_NAMES.keys())

    # ---- 1. z_axis_intercept_mm, full range (+/-Z0_RANGE_MM) -------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bv = bib_hits["z_axis_intercept_mm"][bib_sys == s]
        bv = bv[np.isfinite(bv)]
        bv = bv[np.abs(bv) <= Z0_RANGE_MM]
        sv = sig_hits["z_axis_intercept_mm"][sig_sys == s]
        sv = sv[np.isfinite(sv)]
        sv = sv[np.abs(sv) <= Z0_RANGE_MM]
        overlay_hist(ax, bv, sv, bins=np.linspace(-Z0_RANGE_MM, Z0_RANGE_MM, 121),
                     n_sig_events=n_sig_events)
        draw_symmetric_cut_lines(ax, cuts, "z_axis_intercept_mm", xmax=Z0_RANGE_MM)
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("z-axis intercept of meridian-plane track (mm)")
        ax.set_ylabel("hits / collision / bin")
        ax.legend(fontsize=7)
    plt.suptitle("Z-axis intercept: BIB vs. signal, both in hits/collision/bin "
                 "(BIB = merged plus+minus file = one collision's background; "
                 "signal = per-event average, see docstring)", fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "z_axis_intercept_per_subsystem_with_signal.png", dpi=130)
    plt.close(fig)

    # ---- 2. z_axis_intercept_mm, ZOOMED to +/-Z0_ZOOM_RANGE_MM ------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bv = bib_hits["z_axis_intercept_mm"][bib_sys == s]
        bv = bv[np.isfinite(bv)]
        bv = bv[np.abs(bv) <= Z0_ZOOM_RANGE_MM]
        sv = sig_hits["z_axis_intercept_mm"][sig_sys == s]
        sv = sv[np.isfinite(sv)]
        sv = sv[np.abs(sv) <= Z0_ZOOM_RANGE_MM]
        overlay_hist(ax, bv, sv, bins=np.linspace(-Z0_ZOOM_RANGE_MM, Z0_ZOOM_RANGE_MM, 121),
                     n_sig_events=n_sig_events)
        draw_symmetric_cut_lines(ax, cuts, "z_axis_intercept_mm", xmax=Z0_ZOOM_RANGE_MM)
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("z-axis intercept of meridian-plane track (mm)")
        ax.set_ylabel("hits / collision / bin")
        ax.legend(fontsize=7)
    plt.suptitle(f"Z-axis intercept, zoomed to +/-{Z0_ZOOM_RANGE_MM:.0f}mm: "
                 "BIB vs. signal, both in hits/collision/bin", fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "z_axis_intercept_per_subsystem_zoom_with_signal.png", dpi=130)
    plt.close(fig)

    # ---- 3. inv_radius_per_mm, full range, relabeled in pT GeV/c ---------
    INV_R_MAX = 80.0
    pt_ticks, pt_labels = pt_ticks_for_axis(INV_R_MAX, step=20.0)
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bv = bib_hits["inv_radius_per_mm"][bib_sys == s] * 1000.0
        sv = sig_hits["inv_radius_per_mm"][sig_sys == s] * 1000.0
        overlay_hist(ax, bv, sv, bins=np.linspace(-INV_R_MAX, INV_R_MAX, 161),
                     n_sig_events=n_sig_events)
        draw_momentum_cut_lines(ax, cuts, xmax=INV_R_MAX)
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xticks(pt_ticks)
        ax.set_xticklabels(pt_labels)
        ax.set_xlabel("p$_T$ (GeV/c)")
        ax.set_ylabel("hits / collision / bin")
        ax.legend(fontsize=7)
    plt.suptitle(
        f"Transverse momentum p$_T$ = 0.3 B R (B={B_FIELD_T:.0f}T): BIB vs. signal, "
        "both in hits/collision/bin (same nonlinear pT relabeling as the "
        "BIB-only plot)", fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(outdir / "inv_radius_per_subsystem_with_signal.png", dpi=130)
    plt.close(fig)

    # ---- 4. inv_radius_per_mm, ZOOMED to |pT| >= PT_ZOOM_RANGE_GEV -------
    INV_R_ZOOM_MAX = GEV_PER_INV_M / PT_ZOOM_RANGE_GEV
    pt_zoom_ticks, pt_zoom_labels = pt_ticks_for_axis(
        INV_R_ZOOM_MAX, step=INV_R_ZOOM_MAX / 5.0)
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bx = bib_hits["inv_radius_per_mm"][bib_sys == s] * 1000.0
        bv = bx[np.abs(bx) <= INV_R_ZOOM_MAX]
        sx = sig_hits["inv_radius_per_mm"][sig_sys == s] * 1000.0
        sv = sx[np.abs(sx) <= INV_R_ZOOM_MAX]
        overlay_hist(ax, bv, sv, bins=np.linspace(-INV_R_ZOOM_MAX, INV_R_ZOOM_MAX, 121),
                     n_sig_events=n_sig_events)
        draw_momentum_cut_lines(ax, cuts, xmax=INV_R_ZOOM_MAX)
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xticks(pt_zoom_ticks)
        ax.set_xticklabels(pt_zoom_labels)
        ax.set_xlabel("p$_T$ (GeV/c)")
        ax.set_ylabel("hits / collision / bin")
        ax.legend(fontsize=7)
    plt.suptitle(
        f"Transverse momentum, zoomed to |p$_T$| >= {PT_ZOOM_RANGE_GEV:.0f} GeV/c: "
        "BIB vs. signal, both in hits/collision/bin", fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(outdir / "inv_radius_per_subsystem_zoom_with_signal.png", dpi=130)
    plt.close(fig)

    # ---- 5. t_corrected_ns, full +/-TC_RANGE_NS range ---------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bv = bib_hits["t_corrected_ns"][bib_sys == s]
        bv = bv[np.isfinite(bv)]
        sv = sig_hits["t_corrected_ns"][sig_sys == s]
        sv = sv[np.isfinite(sv)]
        overlay_hist(ax, bv, sv, bins=np.linspace(-TC_RANGE_NS, TC_RANGE_NS, 161),
                     n_sig_events=n_sig_events)
        draw_symmetric_cut_lines(ax, cuts, "t_corrected_ns", xmax=TC_RANGE_NS)
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("t - t$_{expected}$(TOF from IP) (ns)")
        ax.set_ylabel("hits / collision / bin")
        ax.legend(fontsize=7)
    plt.suptitle("Time-of-flight-corrected hit time: BIB vs. signal, both in "
                 "hits/collision/bin", fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "time_corrected_per_subsystem_with_signal.png", dpi=130)
    plt.close(fig)

    # ---- 6. t_corrected_ns, ZOOMED to +/-TC_ZOOM_RANGE_NS ------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bv = bib_hits["t_corrected_ns"][bib_sys == s]
        bv = bv[np.isfinite(bv)]
        bv = bv[np.abs(bv) <= TC_ZOOM_RANGE_NS]
        sv = sig_hits["t_corrected_ns"][sig_sys == s]
        sv = sv[np.isfinite(sv)]
        sv = sv[np.abs(sv) <= TC_ZOOM_RANGE_NS]
        overlay_hist(ax, bv, sv, bins=np.linspace(-TC_ZOOM_RANGE_NS, TC_ZOOM_RANGE_NS, 121),
                     n_sig_events=n_sig_events)
        draw_symmetric_cut_lines(ax, cuts, "t_corrected_ns", xmax=TC_ZOOM_RANGE_NS)
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("t - t$_{expected}$(TOF from IP) (ns)")
        ax.set_ylabel("hits / collision / bin")
        ax.legend(fontsize=7)
    plt.suptitle(f"Time-of-flight-corrected hit time, zoomed to +/-{TC_ZOOM_RANGE_NS:.0f}ns: "
                 "BIB vs. signal, both in hits/collision/bin", fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "time_corrected_per_subsystem_zoom_with_signal.png", dpi=130)
    plt.close(fig)

    print(f"\nWrote 6 plots to {_display_path(outdir.resolve())}")


if __name__ == "__main__":
    main()
