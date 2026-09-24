"""
Step 4 (part 3): "N-1" cut-validation plots. For each of the three cut
variables (t_corrected_ns, z_axis_intercept_mm, and the curvature-based
momentum), redo its BIB-vs-signal distribution from
signal_overlay_angle_plots.py, but with all the OTHER enabled cuts
applied first and this one variable's own cut left off - the standard
"N-1" technique: it shows whether the cut threshold actually sits where
signal and BIB separate, without the plot being shaped by the very cut
being judged. With our 3 cuts, "N-1" means 2 of the 3 applied at a time:
  - z_axis_intercept plots    -> t_corrected_ns AND momentum cuts applied
  - inv_radius/pT plots       -> t_corrected_ns AND z_axis_intercept cuts applied
  - t_corrected_ns plots      -> z_axis_intercept AND momentum cuts applied
A disabled cut in cuts_config.txt contributes no restriction either way
(bib_common.apply_cuts() already returns an all-True mask for it), so
this naturally reduces to "apply whichever of the other cuts are
currently enabled" if fewer than 3 are turned on.

Produces the same 6 plots (full range + zoom, each) as
signal_overlay_angle_plots.py, in the same hits/collision/bin BIB-vs-
signal convention, plus vertical dashed lines marking the (enabled) cut
threshold for the variable being plotted, so the current cut can be
judged by eye against the N-1 distribution:
  - z_axis_intercept_per_subsystem_n1.png       (full range, +/-3000mm)
  - z_axis_intercept_per_subsystem_zoom_n1.png  (+/- z_axis_intercept_mm's
                                                  zoom_halfwidth, mm)
  - inv_radius_per_subsystem_n1.png             (full range, pT-relabeled)
  - inv_radius_per_subsystem_zoom_n1.png        (+/- momentum_gev's
                                                  zoom_halfwidth, GeV/c)
  - time_corrected_per_subsystem_n1.png         (full +/-20ns range)
  - time_corrected_per_subsystem_zoom_n1.png    (+/- t_corrected_ns's
                                                  zoom_halfwidth, ns)

Usage:
    python3 n1_cut_plots.py <bib.root> <signal.root> [cuts_config] [output_dir]
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
    load_cuts, apply_cuts, prepare_output_dir, _display_path, SYSTEM_NAMES,
    load_smearing_config, smearing_rng, describe_smearing,
    apply_position_time_smearing,
    apply_angle_smearing,
)
from incidence_angle_plots import (
    pt_ticks_for_axis, panel_grid,
    Z0_RANGE_MM, B_FIELD_T, GEV_PER_INV_M,
)
from time_of_flight_plots import TC_RANGE_NS
from signal_overlay_angle_plots import overlay_hist

# Zoom half-ranges (z_axis_intercept_mm/momentum_gev/t_corrected_ns) come
# from each cut's own `zoom_halfwidth` in cuts_config.txt (see
# bib_common.load_cuts/ZOOM_HALFWIDTH_DEFAULTS) rather than being fixed
# here, so the zoom window can be widened/narrowed alongside a cut's
# halfwidth without a code change (e.g. so the cut-threshold line stays
# inside the visible plot when a halfwidth is loosened).

CUT_LINE_COLOR = "#a83232"


def n1_mask(per_cut_masks, skip_name):
    """AND of every OTHER cut's own mask (disabled cuts are already
    all-True in per_cut_masks, so they drop out on their own)."""
    n = len(next(iter(per_cut_masks.values())))
    mask = np.ones(n, dtype=bool)
    for name, m in per_cut_masks.items():
        if name != skip_name:
            mask &= m
    return mask


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


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <bib.root> <signal.root> "
              f"[cuts_config] [output_dir]")
        sys.exit(1)
    bib_file = sys.argv[1]
    signal_file = sys.argv[2]
    cuts_config = sys.argv[3] if len(sys.argv) > 3 else \
        str(Path(__file__).resolve().parent / "cuts_config.txt")
    outdir = Path(sys.argv[4] if len(sys.argv) > 4 else "../output_n1_cuts")
    outdir = prepare_output_dir(outdir)

    smearing_config = os.environ.get("SMEARING_CONFIG", "").strip()
    smear_cfg = None
    smear_rng = None
    if smearing_config:
        smear_cfg = load_smearing_config(smearing_config)
        smear_rng = smearing_rng(smear_cfg)
        joined = describe_smearing(smear_cfg)
        print(f"Smearing config: {_display_path(smearing_config)} "
              f"({joined})")

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

    print(f"\nLoading cuts: {_display_path(cuts_config)}")
    cuts = load_cuts(cuts_config)
    for name in ("t_corrected_ns", "z_axis_intercept_mm", "momentum_gev"):
        c = cuts[name]
        state = "enabled" if c["enabled"] else "DISABLED"
        print(f"  {name:20s}  center={c['center']:+.4g}  halfwidth={c['halfwidth']:.4g}  "
              f"zoom_halfwidth={c['zoom_halfwidth']:.4g}  [{state}]")

    Z0_ZOOM_RANGE_MM = cuts["z_axis_intercept_mm"]["zoom_halfwidth"]
    PT_ZOOM_RANGE_GEV = cuts["momentum_gev"]["zoom_halfwidth"]
    TC_ZOOM_RANGE_NS = cuts["t_corrected_ns"]["zoom_halfwidth"]

    _, bib_per_cut = apply_cuts(bib_hits, cuts)
    _, sig_per_cut = apply_cuts(sig_hits, cuts)

    # N-1 masks: for each variable, the AND of the other two cuts' masks.
    bib_n1_time = n1_mask(bib_per_cut, "t_corrected_ns")
    sig_n1_time = n1_mask(sig_per_cut, "t_corrected_ns")
    bib_n1_z = n1_mask(bib_per_cut, "z_axis_intercept_mm")
    sig_n1_z = n1_mask(sig_per_cut, "z_axis_intercept_mm")
    bib_n1_p = n1_mask(bib_per_cut, "momentum_gev")
    sig_n1_p = n1_mask(sig_per_cut, "momentum_gev")

    for label, name, n in (
        ("z_axis_intercept plots", "z_axis_intercept_mm", None),
        ("inv_radius/pT plots", "momentum_gev", None),
        ("t_corrected_ns plots", "t_corrected_ns", None),
    ):
        others = [c for c in ("t_corrected_ns", "z_axis_intercept_mm", "momentum_gev") if c != name]
        print(f"  {label}: applying {others[0]} + {others[1]} (not {name})")

    bib_sys = bib_hits["system"]
    sig_sys = sig_hits["system"]
    sys_ids = sorted(SYSTEM_NAMES.keys())

    # ---- 1. z_axis_intercept_mm, full range, N-1 (time + momentum applied) --
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bsel = (bib_sys == s) & bib_n1_z
        ssel = (sig_sys == s) & sig_n1_z
        bv = bib_hits["z_axis_intercept_mm"][bsel]
        bv = bv[np.isfinite(bv)]
        bv = bv[np.abs(bv) <= Z0_RANGE_MM]
        sv = sig_hits["z_axis_intercept_mm"][ssel]
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
    plt.suptitle("Z-axis intercept, N-1 (time + momentum cuts applied, not this one): "
                 "BIB vs. signal, hits/collision/bin; dashed line = this cut's threshold",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "z_axis_intercept_per_subsystem_n1.png", dpi=130)
    plt.close(fig)

    # ---- 2. z_axis_intercept_mm, ZOOMED, N-1 --------------------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bsel = (bib_sys == s) & bib_n1_z
        ssel = (sig_sys == s) & sig_n1_z
        bv = bib_hits["z_axis_intercept_mm"][bsel]
        bv = bv[np.isfinite(bv)]
        bv = bv[np.abs(bv) <= Z0_ZOOM_RANGE_MM]
        sv = sig_hits["z_axis_intercept_mm"][ssel]
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
    plt.suptitle(f"Z-axis intercept, zoomed to +/-{Z0_ZOOM_RANGE_MM:.0f}mm, N-1 "
                 "(time + momentum cuts applied): BIB vs. signal, hits/collision/bin",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "z_axis_intercept_per_subsystem_zoom_n1.png", dpi=130)
    plt.close(fig)

    # ---- 3. inv_radius_per_mm, full range, N-1 (time + z applied) ----------
    INV_R_MAX = 80.0
    pt_ticks, pt_labels = pt_ticks_for_axis(INV_R_MAX, step=20.0)
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bsel = (bib_sys == s) & bib_n1_p
        ssel = (sig_sys == s) & sig_n1_p
        bv = bib_hits["inv_radius_per_mm"][bsel] * 1000.0
        sv = sig_hits["inv_radius_per_mm"][ssel] * 1000.0
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
        f"Transverse momentum p$_T$ = 0.3 B R (B={B_FIELD_T:.0f}T), N-1 (time + z-intercept "
        "cuts applied, not this one): BIB vs. signal, hits/collision/bin; "
        "dashed lines = this cut's threshold", fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(outdir / "inv_radius_per_subsystem_n1.png", dpi=130)
    plt.close(fig)

    # ---- 4. inv_radius_per_mm, ZOOMED, N-1 ----------------------------------
    INV_R_ZOOM_MAX = GEV_PER_INV_M / PT_ZOOM_RANGE_GEV
    pt_zoom_ticks, pt_zoom_labels = pt_ticks_for_axis(
        INV_R_ZOOM_MAX, step=INV_R_ZOOM_MAX / 5.0)
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bsel = (bib_sys == s) & bib_n1_p
        ssel = (sig_sys == s) & sig_n1_p
        bx = bib_hits["inv_radius_per_mm"][bsel] * 1000.0
        bv = bx[np.abs(bx) <= INV_R_ZOOM_MAX]
        sx = sig_hits["inv_radius_per_mm"][ssel] * 1000.0
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
        f"Transverse momentum, zoomed to |p$_T$| >= {PT_ZOOM_RANGE_GEV:.0f} GeV/c, N-1 "
        "(time + z-intercept cuts applied): BIB vs. signal, hits/collision/bin", fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(outdir / "inv_radius_per_subsystem_zoom_n1.png", dpi=130)
    plt.close(fig)

    # ---- 5. t_corrected_ns, full range, N-1 (z + momentum applied) ---------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bsel = (bib_sys == s) & bib_n1_time
        ssel = (sig_sys == s) & sig_n1_time
        bv = bib_hits["t_corrected_ns"][bsel]
        bv = bv[np.isfinite(bv)]
        sv = sig_hits["t_corrected_ns"][ssel]
        sv = sv[np.isfinite(sv)]
        overlay_hist(ax, bv, sv, bins=np.linspace(-TC_RANGE_NS, TC_RANGE_NS, 161),
                     n_sig_events=n_sig_events)
        draw_symmetric_cut_lines(ax, cuts, "t_corrected_ns", xmax=TC_RANGE_NS)
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("t - t$_{expected}$(TOF from IP) (ns)")
        ax.set_ylabel("hits / collision / bin")
        ax.legend(fontsize=7)
    plt.suptitle("Time-of-flight-corrected hit time, N-1 (z-intercept + momentum cuts "
                 "applied, not this one): BIB vs. signal, hits/collision/bin; "
                 "dashed line = this cut's threshold", fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "time_corrected_per_subsystem_n1.png", dpi=130)
    plt.close(fig)

    # ---- 6. t_corrected_ns, ZOOMED, N-1 -------------------------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bsel = (bib_sys == s) & bib_n1_time
        ssel = (sig_sys == s) & sig_n1_time
        bv = bib_hits["t_corrected_ns"][bsel]
        bv = bv[np.isfinite(bv)]
        bv = bv[np.abs(bv) <= TC_ZOOM_RANGE_NS]
        sv = sig_hits["t_corrected_ns"][ssel]
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
    plt.suptitle(f"Time-of-flight-corrected hit time, zoomed to +/-{TC_ZOOM_RANGE_NS:.0f}ns, "
                 "N-1 (z-intercept + momentum cuts applied): BIB vs. signal, "
                 "hits/collision/bin", fontsize=10)
    plt.tight_layout()
    plt.savefig(outdir / "time_corrected_per_subsystem_zoom_n1.png", dpi=130)
    plt.close(fig)

    print(f"\nWrote 6 N-1 plots to {_display_path(outdir.resolve())}")


if __name__ == "__main__":
    main()
