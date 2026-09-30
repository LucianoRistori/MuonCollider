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
Each panel is one subsystem, and shows that subsystem's own cut on the
plotted variable (__cuts_config.txt sets the cuts per subsystem) as
dashed red lines, with its value in the panel title.

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

Display: the signal is typically several orders of magnitude below BIB,
so in each panel it is drawn multiplied by a power of ten that brings
the two distributions to the same height on the log scale (see
signal_scale(); each one's range reaches at most MAX_RANGE_DECADES below
its peak), shown in the legend box ("Signal x 10^n"). The left
axis is in BIB units, the right axis in the signal's own, unscaled
units - both hits / collision / bin.

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
from matplotlib.ticker import LogFormatterSciNotation

from bib_common import (
    add_grid,
    load_hits, add_incidence_angles, add_time_of_flight,
    prepare_output_dir, _display_path, SYSTEM_NAMES, load_cuts,
    load_smearing_config, smearing_rng, describe_smearing,
    under_run_all, short_path, loaded_line, print_table, default_cuts_config,
    apply_position_time_smearing,
    apply_angle_smearing,
    signal_muon_hits_only,
)
from incidence_angle_plots import (
    pt_ticks_for_axis, panel_grid,
    Z0_RANGE_MM, Z0_ZOOM_RANGE_MM, PT_ZOOM_RANGE_GEV, B_FIELD_T, GEV_PER_INV_M,
)
from time_of_flight_plots import TC_RANGE_NS, TC_ZOOM_RANGE_NS

BIB_COLOR = "#3b7dd8"
SIG_COLOR = "#2ca858"

CUT_LINE_COLOR = "#a83232"


def draw_symmetric_cut_lines(ax, cuts, name, system, xmax=None):
    """For a |value| <= limit cut (t_corrected_ns, z_axis_intercept_mm):
    draw vertical lines at +/- the limit of this panel's own subsystem
    (cuts are set per subsystem), unless the cut is off there or
    (optionally) the lines fall outside the visible range."""
    limit = cuts[name]["per_system"].get(system)
    if limit is None:
        return
    for edge in (-limit, limit):
        if xmax is None or abs(edge) <= xmax:
            ax.axvline(edge, color=CUT_LINE_COLOR, linewidth=1.1,
                       linestyle="--", alpha=0.8)


def draw_momentum_cut_lines(ax, cuts, system, xmax=None):
    """momentum_gev cut accepts |inv_radius_per_mm| <= threshold (a band
    around 1/R=0); draw its two edges on the curvature axis (in 1/m,
    matching the plotted units bv = inv_radius_per_mm*1000), for this
    panel's own subsystem."""
    pt_min = cuts["momentum_gev"]["per_system"].get(system)
    if not pt_min:          # off (None), or pT >= 0, which keeps every hit
        return
    inv_radius_threshold_per_m = GEV_PER_INV_M / pt_min
    for edge in (-inv_radius_threshold_per_m, inv_radius_threshold_per_m):
        if xmax is None or abs(edge) <= xmax:
            ax.axvline(edge, color=CUT_LINE_COLOR, linewidth=1.1,
                       linestyle="--", alpha=0.8)


def panel_title(cuts, name, system):
    """A panel's title: its subsystem, and that subsystem's own value of
    the cut on the plotted variable, e.g. 'OT endcap   (cut |z0| <= 40 mm)'."""
    x = cuts[name]["per_system"].get(system)
    if name == "t_corrected_ns":
        cut = "no time cut" if x is None else f"cut |t - t$_{{TOF}}$| \u2264 {x:g} ns"
    elif name == "z_axis_intercept_mm":
        cut = "no z$_0$ cut" if x is None else f"cut |z$_0$| \u2264 {x:g} mm"
    else:
        cut = "no p$_T$ cut" if x is None else f"cut p$_T$ \u2265 {x:g} GeV/c"
    return f"{SYSTEM_NAMES[system]}   ({cut})"


Y_MARGIN = 0.05            # of the y span, left free above and below the histograms
# The y axis spans at least this many decades, so that it always shows at
# least two labelled ticks (minor ticks are labelled below 2 decades):
MIN_SPAN_DECADES = 0.5
# Each histogram's range goes down at most this many decades below its
# highest bin: lower bins are left out of the range, and so of the plot
# (None: the full range, down to the smallest non-empty bin).
MAX_RANGE_DECADES = 3.0
# for the figure titles:
SCALE_NOTE = " (signal $\\times$ 10$^n$ as marked; right axis unscaled)"


def signal_scale(bib_counts, sig_counts, margin=Y_MARGIN):
    """
    For one panel: the power of ten n such that the signal, multiplied
    by 10**n, sits at the same height as BIB on the log scale, and the
    y-axis limits (in BIB units). Each histogram's range runs from its
    smallest to its largest non-empty bin, but at most MAX_RANGE_DECADES
    below the largest; 10**n brings the middle of the
    signal range (on the log scale) to the middle of the BIB range,
    rounded to the nearest power of ten. The axis covers both ranges -
    the larger of the two spans, plus up to half a decade from the
    rounding - with `margin` of the span left free above and below (and
    at least MIN_SPAN_DECADES). n = 0 if either histogram is empty.
    Returns (n, (ymin, ymax)), or (0, None) if both are empty.
    """
    def log_range(counts):          # (lowest, highest) log10 of the range, or None
        v = np.log10(counts[counts > 0])
        if not v.size:
            return None
        hi = v.max()
        lo = v.min() if MAX_RANGE_DECADES is None else max(v.min(), hi - MAX_RANGE_DECADES)
        return lo, hi
    b, s = log_range(bib_counts), log_range(sig_counts)
    n = 0
    if b and s:
        n = int(np.round((b[0] + b[1]) / 2 - (s[0] + s[1]) / 2))
    shown = [x for x in (b, s and (s[0] + n, s[1] + n)) if x]
    if not shown:
        return 0, None
    lo, hi = min(x[0] for x in shown), max(x[1] for x in shown)
    span = max(hi - lo, MIN_SPAN_DECADES)
    mid = (lo + hi) / 2
    half = (0.5 + margin) * span
    return n, (10.0 ** (mid - half), 10.0 ** (mid + half))


def times_power_of_ten(n):
    """'1', '10', or '10$^{n}$' (mathtext)."""
    return "1" if n == 0 else "10" if n == 1 else f"10$^{{{n}}}$"


def overlay_hist(ax, bib_v, sig_v, bins, n_sig_events):
    """
    One panel: BIB and signal histograms, both in hits / collision / bin -
    BIB as raw hit counts per bin (the merged plus+minus file is taken to
    be one collision's full background), signal as hit counts per bin
    divided by the number of signal events, i.e. the average
    contribution of the one signal muon per collision (see module
    docstring). Log y scale. The signal is drawn multiplied by 10**n from
    signal_scale(), with "Signal x 10^n" at the top of the legend box; the
    left axis is in BIB units, the right one in the signal's own, unscaled
    units. Returns n.
    """
    edges = np.asarray(bins, dtype=float)
    bib_counts = np.histogram(bib_v, bins=edges)[0].astype(float)
    sig_counts = np.histogram(sig_v, bins=edges)[0] / n_sig_events
    n, ylim = signal_scale(bib_counts, sig_counts)
    k = 10.0 ** n
    # (hist() of the bin edges, weighted by the contents: a step outline
    # of already-filled histograms)
    ax.hist(edges[:-1], bins=edges, weights=bib_counts, histtype="step",
            color=BIB_COLOR, linewidth=1.4, label="BIB (left scale)")
    ax.hist(edges[:-1], bins=edges, weights=sig_counts * k, histtype="step",
            color=SIG_COLOR, linewidth=1.4, label="signal (right scale)")
    ax.set_yscale("log")
    right = ax.twinx()
    right.set_yscale("log")
    if ylim is not None:
        ax.set_ylim(*ylim)
        right.set_ylim(ylim[0] / k, ylim[1] / k)
    ax.set_ylabel("BIB hits / collision / bin", color=BIB_COLOR)
    ax.tick_params(axis="y", which="both", labelcolor=BIB_COLOR)
    right.set_ylabel("signal hits / collision / bin", color=SIG_COLOR)
    right.tick_params(axis="y", which="both", labelcolor=SIG_COLOR)
    for a in (ax, right):   # label minor ticks (2, 3, 4, 6 x 10^n) below 2 decades
        a.yaxis.set_minor_formatter(
            LogFormatterSciNotation(labelOnlyBase=False, minor_thresholds=(2, 0.5)))
    if not sig_counts.any():
        label = "no signal hits"
    else:
        label = f"Signal \u00d7 {times_power_of_ten(n)}"
        if not bib_counts.any():
            label += "  (no BIB hits)"
    # the factor heads the legend box, which goes where it covers the least
    # of the histograms (the curves can reach any corner)
    legend = ax.legend(fontsize=7, loc="best", title=label, title_fontsize=10)
    legend.get_title().set_color(SIG_COLOR)
    return n


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <bib.root> <signal.root> "
              f"[output_dir] [cuts_config]")
        sys.exit(1)
    bib_file = sys.argv[1]
    signal_file = sys.argv[2]
    outdir = Path(sys.argv[3] if len(sys.argv) > 3 else "../output_signal_angle_overlay")
    cuts_config = sys.argv[4] if len(sys.argv) > 4 else \
        default_cuts_config()
    outdir = prepare_output_dir(outdir)

    cuts = load_cuts(cuts_config)
    if not under_run_all():
        print(f"Cuts (for the reference lines): {short_path(cuts_config)}")

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
    add_incidence_angles(bib_hits)
    if smear_cfg is not None:
        apply_angle_smearing(bib_hits, smear_cfg, smear_rng)
    add_time_of_flight(bib_hits)
    n_bib_hit = len(bib_hits["x"])
    print(loaded_line(bib_file, bib_hits, "BIB"))
    sig_hits = load_hits(signal_file, muon_hits_only=signal_muon_hits_only(cuts_config))
    if smear_cfg is not None:
        apply_position_time_smearing(sig_hits, smear_cfg, smear_rng)
    add_incidence_angles(sig_hits)
    if smear_cfg is not None:
        apply_angle_smearing(sig_hits, smear_cfg, smear_rng)
    add_time_of_flight(sig_hits)
    n_sig_hit = len(sig_hits["x"])
    n_sig_events = sig_hits["_n_events"]
    print(loaded_line(signal_file, sig_hits, "signal"))

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
        draw_symmetric_cut_lines(ax, cuts, "z_axis_intercept_mm", s, xmax=Z0_RANGE_MM)
        ax.set_title(panel_title(cuts, "z_axis_intercept_mm", s))
        ax.set_xlabel("z-axis intercept of meridian-plane track (mm)")
    plt.suptitle("Z-axis intercept: BIB vs. signal, both in hits/collision/bin" + SCALE_NOTE,
                 fontsize=10)
    plt.tight_layout()
    add_grid()
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
        draw_symmetric_cut_lines(ax, cuts, "z_axis_intercept_mm", s, xmax=Z0_ZOOM_RANGE_MM)
        ax.set_title(panel_title(cuts, "z_axis_intercept_mm", s))
        ax.set_xlabel("z-axis intercept of meridian-plane track (mm)")
    plt.suptitle(f"Z-axis intercept, zoomed to +/-{Z0_ZOOM_RANGE_MM:.0f}mm: "
                 "BIB vs. signal, both in hits/collision/bin" + SCALE_NOTE, fontsize=10)
    plt.tight_layout()
    add_grid()
    plt.savefig(outdir / "z_axis_intercept_per_subsystem_zoom_with_signal.png", dpi=130)
    plt.close(fig)

    # ---- 3. inv_radius_per_mm, full range, relabeled in pT GeV/c ---------
    INV_R_MAX = 80.0
    pt_ticks, pt_labels = pt_ticks_for_axis(INV_R_MAX, step=40.0)   # (labels are long)
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bv = bib_hits["inv_radius_per_mm"][bib_sys == s] * 1000.0
        sv = sig_hits["inv_radius_per_mm"][sig_sys == s] * 1000.0
        overlay_hist(ax, bv, sv, bins=np.linspace(-INV_R_MAX, INV_R_MAX, 161),
                     n_sig_events=n_sig_events)
        draw_momentum_cut_lines(ax, cuts, s, xmax=INV_R_MAX)
        ax.set_title(panel_title(cuts, "momentum_gev", s))
        ax.set_xticks(pt_ticks)
        ax.set_xticklabels(pt_labels)
        ax.set_xlabel("p$_T$ (GeV/c)")
    plt.suptitle(
        f"Transverse momentum p$_T$ = 0.3 B R (B={B_FIELD_T:.0f}T): BIB vs. signal, "
        "both in hits/collision/bin" + SCALE_NOTE + "; same nonlinear pT relabeling as the "
        "BIB-only plot", fontsize=10,
    )
    plt.tight_layout()
    add_grid()
    plt.savefig(outdir / "inv_radius_per_subsystem_with_signal.png", dpi=130)
    plt.close(fig)

    # ---- 4. inv_radius_per_mm, ZOOMED to |pT| >= PT_ZOOM_RANGE_GEV -------
    INV_R_ZOOM_MAX = GEV_PER_INV_M / PT_ZOOM_RANGE_GEV
    pt_zoom_ticks, pt_zoom_labels = pt_ticks_for_axis(
        INV_R_ZOOM_MAX, step=INV_R_ZOOM_MAX / 4.0)
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        bx = bib_hits["inv_radius_per_mm"][bib_sys == s] * 1000.0
        bv = bx[np.abs(bx) <= INV_R_ZOOM_MAX]
        sx = sig_hits["inv_radius_per_mm"][sig_sys == s] * 1000.0
        sv = sx[np.abs(sx) <= INV_R_ZOOM_MAX]
        overlay_hist(ax, bv, sv, bins=np.linspace(-INV_R_ZOOM_MAX, INV_R_ZOOM_MAX, 121),
                     n_sig_events=n_sig_events)
        draw_momentum_cut_lines(ax, cuts, s, xmax=INV_R_ZOOM_MAX)
        ax.set_title(panel_title(cuts, "momentum_gev", s))
        ax.set_xticks(pt_zoom_ticks)
        ax.set_xticklabels(pt_zoom_labels)
        ax.set_xlabel("p$_T$ (GeV/c)")
    plt.suptitle(
        f"Transverse momentum, zoomed to |p$_T$| >= {PT_ZOOM_RANGE_GEV:.0f} GeV/c: "
        "BIB vs. signal, both in hits/collision/bin" + SCALE_NOTE, fontsize=10,
    )
    plt.tight_layout()
    add_grid()
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
        draw_symmetric_cut_lines(ax, cuts, "t_corrected_ns", s, xmax=TC_RANGE_NS)
        ax.set_title(panel_title(cuts, "t_corrected_ns", s))
        ax.set_xlabel("t - t$_{expected}$(TOF from IP) (ns)")
    plt.suptitle("Time-of-flight-corrected hit time: BIB vs. signal, both in "
                 "hits/collision/bin" + SCALE_NOTE, fontsize=10)
    plt.tight_layout()
    add_grid()
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
        draw_symmetric_cut_lines(ax, cuts, "t_corrected_ns", s, xmax=TC_ZOOM_RANGE_NS)
        ax.set_title(panel_title(cuts, "t_corrected_ns", s))
        ax.set_xlabel("t - t$_{expected}$(TOF from IP) (ns)")
    plt.suptitle(f"Time-of-flight-corrected hit time, zoomed to +/-{TC_ZOOM_RANGE_NS:.0f}ns: "
                 "BIB vs. signal, both in hits/collision/bin" + SCALE_NOTE, fontsize=10)
    plt.tight_layout()
    add_grid()
    plt.savefig(outdir / "time_corrected_per_subsystem_zoom_with_signal.png", dpi=130)
    plt.close(fig)

    print(f"Wrote 6 plots to {short_path(outdir)}/")


if __name__ == "__main__":
    main()
