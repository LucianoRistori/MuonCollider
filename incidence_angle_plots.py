"""
Step 2 (part 1): compute and visualize the track incidence angle at each
hit, decomposed into a longitudinal and transverse component relative to
the local sensor-plane normal (see bib_common.add_incidence_angles() for
the exact definition). This is a diagnostic/exploratory step to confirm
the two angle definitions behave sensibly, before using them to define
selection cuts and re-examine hit densities.

Usage:
    python3 incidence_angle_plots.py <input.root> [output_dir]
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
    load_hits, add_incidence_angles, prepare_output_dir, _display_path,
    SYSTEM_NAMES,
    load_smearing_config, smearing_rng, apply_position_time_smearing,
    apply_angle_smearing,
)

ANGLE_KEYS = ("theta_long_deg", "theta_trans_deg", "theta_full_deg")
Z0_RANGE_MM = 3000.0  # plot range for z_axis_intercept_mm; see note at its histogram
Z0_ZOOM_RANGE_MM = 100.0  # +/-10cm zoomed range for z_axis_intercept_mm
PT_ZOOM_RANGE_GEV = 1.0  # +/-1 GeV/c zoomed range for transverse momentum

# For the inv_radius_per_subsystem.png x-axis: pT[GeV/c] = 0.3*B[T]*R[m],
# so with x = 1/R in 1/m, pT = GEV_PER_INV_M / |x|. Assumed field, used only
# for this axis label (the underlying curvature data/binning is unaffected).
B_FIELD_T = 5.0
GEV_PER_INV_M = 0.3 * B_FIELD_T


def pt_ticks_for_axis(x_max, step):
    """Tick positions (in 1/R, 1/m), evenly spaced by `step` as a normal
    linear curvature axis would be, each relabeled with its corresponding
    pT = GEV_PER_INV_M/|x| in GeV/c. Ticks stay evenly spaced in 1/R (so
    they don't crowd together) - only the printed numbers are nonlinear,
    since pT itself is a reciprocal, not linear, function of 1/R. The
    center tick (1/R=0, a perfectly straight/infinite-pT track) is labeled
    with an infinity symbol rather than a number."""
    ticks = np.arange(-x_max, x_max + 1e-9, step)
    labels = []
    for x in ticks:
        if abs(x) < 1e-9:
            labels.append(r"$\pm\infty$")
        else:
            pt = GEV_PER_INV_M / abs(x)
            labels.append(f"{'-' if x < 0 else ''}{pt:.3g}")
    return ticks, labels


def panel_grid():
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    return fig, axes.flatten()


def main():
    root_file = sys.argv[1] if len(sys.argv) > 1 else \
        "/mnt/user-data/uploads/ntu_bib_plus_1evt.root"
    outdir = Path(sys.argv[2] if len(sys.argv) > 2 else "../output_incidence")
    outdir = prepare_output_dir(outdir)

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

    hits = load_hits(root_file)
    if smear_cfg is not None:
        apply_position_time_smearing(hits, smear_cfg, smear_rng)
    n_hit = len(hits["x"])
    add_incidence_angles(hits)
    if smear_cfg is not None:
        apply_angle_smearing(hits, smear_cfg, smear_rng)
    system = hits["system"]
    sys_ids = sorted(SYSTEM_NAMES.keys())

    # ---- 1. theta_long_deg, one panel per subsystem ------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        v = hits["theta_long_deg"][system == s]
        ax.hist(v, bins=np.linspace(-180, 180, 121), color="#3b7dd8")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("theta_long (deg)")
        ax.set_ylabel("hits")
    plt.suptitle("Longitudinal incidence angle (meridian-plane, relative to sensor normal)")
    plt.tight_layout()
    plt.savefig(outdir / "theta_long_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 1b. z_axis_intercept_mm, companion to theta_long, one panel per
    # subsystem: for each hit, its meridian-plane (rho, z) trajectory -
    # a straight line through the hit position with local slope
    # (p_rho, p_z) - extrapolated back to where it crosses the Z axis
    # (rho=0). Undefined for tracks running along z within the meridian
    # plane (p_rho ~ 0); those are dropped (reported as excluded below)
    # rather than plotted at +/-infinity. A small tail of large-but-finite
    # values (near-grazing tracks) also falls outside +/-Z0_RANGE_MM and
    # is excluded from the histogram (not clipped into the edge bins) -
    # the excluded fraction is printed per subsystem for transparency.
    fig, axes = panel_grid()
    z0_excluded_frac = {}
    for ax, s in zip(axes, sys_ids):
        v = hits["z_axis_intercept_mm"][system == s]
        finite = v[np.isfinite(v)]
        in_range = finite[np.abs(finite) <= Z0_RANGE_MM]
        z0_excluded_frac[s] = 1.0 - (in_range.size / v.size if v.size else 0.0)
        ax.hist(in_range, bins=np.linspace(-Z0_RANGE_MM, Z0_RANGE_MM, 121),
                color="#8a5bc7")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("z-axis intercept of meridian-plane track (mm)")
        ax.set_ylabel("hits")
    plt.suptitle("Z where the hit's local (meridian-plane) track direction "
                 "crosses the Z axis")
    plt.tight_layout()
    plt.savefig(outdir / "z_axis_intercept_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 1c. z_axis_intercept_mm, ZOOMED to +/-Z0_ZOOM_RANGE_MM ----------
    # Same quantity as above, just a much narrower range with finer bins,
    # to resolve the structure right around z=0 that the full +/-3000mm
    # view compresses into a couple of central bins.
    fig, axes = panel_grid()
    z0_zoom_excluded_frac = {}
    for ax, s in zip(axes, sys_ids):
        v = hits["z_axis_intercept_mm"][system == s]
        finite = v[np.isfinite(v)]
        in_range = finite[np.abs(finite) <= Z0_ZOOM_RANGE_MM]
        z0_zoom_excluded_frac[s] = 1.0 - (in_range.size / v.size if v.size else 0.0)
        ax.hist(in_range, bins=np.linspace(-Z0_ZOOM_RANGE_MM, Z0_ZOOM_RANGE_MM, 121),
                color="#8a5bc7")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("z-axis intercept of meridian-plane track (mm)")
        ax.set_ylabel("hits")
    plt.suptitle(f"Z-axis intercept, zoomed to +/-{Z0_ZOOM_RANGE_MM:.0f}mm "
                 "(same quantity as above, narrower range/finer bins)")
    plt.tight_layout()
    plt.savefig(outdir / "z_axis_intercept_per_subsystem_zoom.png", dpi=130)
    plt.close(fig)

    # ---- 2. theta_trans_deg, one panel per subsystem ------------------
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        v = hits["theta_trans_deg"][system == s]
        ax.hist(v, bins=np.linspace(-180, 180, 121), color="#e05252")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("theta_trans (deg)")
        ax.set_ylabel("hits")
    plt.suptitle("Transverse incidence angle (X-Y plane, relative to sensor normal)")
    plt.tight_layout()
    plt.savefig(outdir / "theta_trans_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 2b. inv_radius_per_mm, companion to theta_trans, one panel per
    # subsystem: curvature 1/R (plotted in 1/m) of the transverse-plane
    # circle through the origin and the hit, tangent there to the track's
    # transverse momentum direction - i.e. the curvature a particle
    # produced at the origin would need to reach this hit going in this
    # transverse direction. Bounded (unlike z_axis_intercept_mm), so no
    # exclusion is needed. Binning stays in 1/R (1/m); the x-axis is
    # relabeled in transverse momentum pT[GeV/c] = 0.3*B*R[m] assuming
    # B=5T (BIB_FIELD_T), since pT vs. 1/R is a reciprocal (nonlinear)
    # relationship, not a rescaling - see pt_ticks_for_axis().
    INV_R_MAX = 80.0
    pt_ticks, pt_labels = pt_ticks_for_axis(INV_R_MAX, step=20.0)
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        v = hits["inv_radius_per_mm"][system == s] * 1000.0  # 1/mm -> 1/m
        ax.hist(v, bins=np.linspace(-INV_R_MAX, INV_R_MAX, 161), color="#c78a2f")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xticks(pt_ticks)
        ax.set_xticklabels(pt_labels)
        ax.set_xlabel("p$_T$ (GeV/c)")
        ax.set_ylabel("hits")
    plt.suptitle(
        f"Transverse curvature, axis relabeled as transverse momentum "
        f"p$_T$ = 0.3 B R (B={B_FIELD_T:.0f}T)\n"
        "tick positions are evenly spaced in 1/R (as the underlying binning is); "
        "only the printed p$_T$ values are nonlinear - sign follows curvature "
        "sign; center (1/R=0) is a straight, infinite-p$_T$ track",
        fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(outdir / "inv_radius_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 2c. transverse curvature, ZOOMED so the pT axis reads to
    # +/-PT_ZOOM_RANGE_GEV at the edges -----------------------------------
    # Same convention as the full-range plot above (still binned and
    # centered in 1/R, decreasing |pT| going outward from the center) -
    # NOT a re-binning directly in pT, which would put pT=0 at the center
    # and make pT increase outward (backwards relative to the full-range
    # plot; an earlier version of this script made that mistake). Here we
    # just crop the curvature axis to the narrow window
    # |1/R| <= GEV_PER_INV_M / PT_ZOOM_RANGE_GEV (where the relabeled pT
    # tick would read +/-PT_ZOOM_RANGE_GEV) and rebin at much finer
    # resolution within it, revealing the internal structure of what was
    # a single central spike in the full-range plot. Hits outside this
    # window (|1/R| larger, i.e. |pT| < PT_ZOOM_RANGE_GEV - the more
    # sharply curved, lower-momentum tail) are excluded from the
    # histogram; the excluded fraction is printed per subsystem.
    INV_R_ZOOM_MAX = GEV_PER_INV_M / PT_ZOOM_RANGE_GEV
    pt_zoom_ticks, pt_zoom_labels = pt_ticks_for_axis(
        INV_R_ZOOM_MAX, step=INV_R_ZOOM_MAX / 5.0)
    fig, axes = panel_grid()
    pt_zoom_excluded_frac = {}
    for ax, s in zip(axes, sys_ids):
        x = hits["inv_radius_per_mm"][system == s] * 1000.0  # 1/mm -> 1/m
        in_range = x[np.abs(x) <= INV_R_ZOOM_MAX]
        pt_zoom_excluded_frac[s] = 1.0 - (in_range.size / x.size if x.size else 0.0)
        ax.hist(in_range, bins=np.linspace(-INV_R_ZOOM_MAX, INV_R_ZOOM_MAX, 121),
                color="#c78a2f")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xticks(pt_zoom_ticks)
        ax.set_xticklabels(pt_zoom_labels)
        ax.set_xlabel("p$_T$ (GeV/c)")
        ax.set_ylabel("hits")
    plt.suptitle(
        f"Transverse momentum p$_T$ = 0.3 B R (B={B_FIELD_T:.0f}T), zoomed to "
        f"the curvature window where |p$_T$| >= {PT_ZOOM_RANGE_GEV:.0f} GeV/c\n"
        "same center-high/edge-low p$_T$ sense as the full-range plot - "
        "just cropped and finely rebinned near the center",
        fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(outdir / "inv_radius_per_subsystem_zoom.png", dpi=130)
    plt.close(fig)

    # ---- 3. theta_full_deg (plain 3D angle vs. normal), for reference -
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        v = hits["theta_full_deg"][system == s]
        ax.hist(v, bins=np.linspace(0, 180, 91), color="#4a9c4a")
        ax.set_yscale("log")
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("theta_full (deg)")
        ax.set_ylabel("hits")
    plt.suptitle("Full 3D incidence angle vs. sensor normal (reference / cross-check)")
    plt.tight_layout()
    plt.savefig(outdir / "theta_full_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 4. 2D: theta_long vs theta_trans, one panel per subsystem ----
    fig, axes = panel_grid()
    for ax, s in zip(axes, sys_ids):
        mask = system == s
        h = ax.hist2d(
            hits["theta_long_deg"][mask], hits["theta_trans_deg"][mask],
            bins=(np.linspace(-180, 180, 91), np.linspace(-180, 180, 91)),
            norm=matplotlib.colors.LogNorm(), cmap="viridis",
        )
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_xlabel("theta_long (deg)")
        ax.set_ylabel("theta_trans (deg)")
        plt.colorbar(h[3], ax=ax, label="hits")
    plt.suptitle("Longitudinal vs. transverse incidence angle")
    plt.tight_layout()
    plt.savefig(outdir / "theta_long_vs_trans_2d.png", dpi=130)
    plt.close(fig)

    # ---- per-subsystem summary CSV + module-tilt diagnostic -----------
    summary_path = outdir / "incidence_angle_summary.csv"
    fields = ["system", "system_name", "n_hits"]
    for k in ANGLE_KEYS + ("n_tilt_deg",):
        fields += [f"{k}_mean", f"{k}_std", f"{k}_p01", f"{k}_p50", f"{k}_p99"]
    fields += ["z_axis_intercept_mm_finite_frac", "z_axis_intercept_mm_p01",
               "z_axis_intercept_mm_p50", "z_axis_intercept_mm_p99"]
    fields += ["inv_radius_per_m_mean", "inv_radius_per_m_std",
               "inv_radius_per_m_p01", "inv_radius_per_m_p50", "inv_radius_per_m_p99"]
    rows = []
    for s in sys_ids:
        mask = system == s
        row = {"system": s, "system_name": SYSTEM_NAMES[s], "n_hits": int(mask.sum())}
        for k in ANGLE_KEYS + ("n_tilt_deg",):
            v = hits[k][mask]
            row[f"{k}_mean"] = float(v.mean())
            row[f"{k}_std"] = float(v.std())
            row[f"{k}_p01"] = float(np.percentile(v, 1))
            row[f"{k}_p50"] = float(np.percentile(v, 50))
            row[f"{k}_p99"] = float(np.percentile(v, 99))
        z0 = hits["z_axis_intercept_mm"][mask]
        z0_finite = z0[np.isfinite(z0)]
        row["z_axis_intercept_mm_finite_frac"] = float(z0_finite.size / z0.size) if z0.size else float("nan")
        for p in (1, 50, 99):
            row[f"z_axis_intercept_mm_p{p:02d}"] = (
                float(np.percentile(z0_finite, p)) if z0_finite.size else float("nan")
            )
        invr = hits["inv_radius_per_mm"][mask] * 1000.0  # 1/mm -> 1/m
        row["inv_radius_per_m_mean"] = float(invr.mean())
        row["inv_radius_per_m_std"] = float(invr.std())
        for p in (1, 50, 99):
            row[f"inv_radius_per_m_p{p:02d}"] = float(np.percentile(invr, p))
        rows.append(row)
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"Loaded {hits['_tree_name']}: n_hit={n_hit:,}")
    print(f"Wrote plots and {summary_path.name} to {_display_path(outdir.resolve())}")
    print()
    print("-- Module out-of-meridian-plane tilt (n_tilt_deg): should be ~0 for a "
          "perfectly axially-symmetric layer; nonzero means theta_long/theta_trans "
          "don't exactly add up (in quadrature, in tan) to theta_full for that hit --")
    for r in rows:
        print(f"  {r['system_name']:12s}  n_tilt mean={r['n_tilt_deg_mean']:6.3f} deg  "
              f"p99={r['n_tilt_deg_p99']:6.3f} deg")
    print()
    print("-- theta_long / theta_trans / theta_full (deg): mean +/- std --")
    for r in rows:
        print(f"  {r['system_name']:12s}  "
              f"long={r['theta_long_deg_mean']:7.2f}+/-{r['theta_long_deg_std']:6.2f}  "
              f"trans={r['theta_trans_deg_mean']:7.2f}+/-{r['theta_trans_deg_std']:6.2f}  "
              f"full={r['theta_full_deg_mean']:7.2f}+/-{r['theta_full_deg_std']:6.2f}")
    print()
    print(f"-- z-axis intercept of meridian-plane track (mm), median [p01, p99] "
          f"(only shown/plotted within +/-{Z0_RANGE_MM:.0f}mm; excluded = undefined "
          f"p_rho~0 tracks + out-of-range tail) --")
    for r in rows:
        s = r["system"]
        print(f"  {r['system_name']:12s}  "
              f"p50={r['z_axis_intercept_mm_p50']:9.1f}  "
              f"[{r['z_axis_intercept_mm_p01']:9.1f}, {r['z_axis_intercept_mm_p99']:9.1f}]  "
              f"excluded={z0_excluded_frac[s]*100:5.2f}%")
    print()
    print("-- transverse curvature 1/R, circle through origin and hit (1/m), "
          "median [p01, p99] --")
    for r in rows:
        print(f"  {r['system_name']:12s}  "
              f"p50={r['inv_radius_per_m_p50']:8.3f}  "
              f"[{r['inv_radius_per_m_p01']:8.3f}, {r['inv_radius_per_m_p99']:8.3f}]")
    print()
    print(f"-- z-axis intercept, ZOOMED to +/-{Z0_ZOOM_RANGE_MM:.0f}mm: fraction of "
          f"hits excluded (outside this narrower window, or undefined) --")
    for r in rows:
        s = r["system"]
        print(f"  {r['system_name']:12s}  excluded={z0_zoom_excluded_frac[s]*100:6.2f}%")
    print()
    print(f"-- transverse momentum p_T (GeV/c, B={B_FIELD_T:.0f}T), ZOOMED to the "
          f"curvature window where |p_T| >= {PT_ZOOM_RANGE_GEV:.0f} GeV/c: fraction "
          f"of hits excluded (|p_T| below this window, i.e. the sharply-curved, "
          f"low-momentum tail) --")
    for r in rows:
        s = r["system"]
        print(f"  {r['system_name']:12s}  excluded={pt_zoom_excluded_frac[s]*100:6.2f}%")


if __name__ == "__main__":
    main()
