"""
fake_rate.py - expected number of fake tracks, E[#fakes], in the 6-layer
IT + OT barrel tower of the fake-rate framework (Sections 14-16 of the
paper; calibration in fake_rate_framework/real_tower_calibration.py),
from the BIB hit density of each IT/OT barrel layer after the cuts.

E[#fakes] = P_true x prod_i n_i,  n_i = (BIB hits after cuts / area)_i x W_i^2
with W_i^2 the plane areas of the projective tower (1 m^2 outermost plane)
and P_true the calibrated probability for 3 track models.

P_true depends steeply on the hit position resolution: it was calibrated at
sigma_ref = 0.1 mm in both coordinates and scales as
    P_true(sigma_u, sigma_v) = P_ref x (sigma_u/sigma_ref)^nb x (sigma_v/sigma_ref)^nd
(paper Sections 5.1, 10.4, 13.2), with nb, nd the degrees of freedom of the
bending-view and depth-view fits: 4 + 4 for the helix models, 5 + 4 for the
origin line (B = 0). In the barrel, u is the r-phi (bending) coordinate and
v the z (depth) coordinate. The resolutions come from __smearing_config.txt
([position] sigma_u_mm, sigma_v_mm); with position smearing off, sigma_ref
is used.

Used by resolution_scan.py (./run_scan), and by ./run_all at the end of
its cuts step:
    python3 fake_rate.py <step4_cuts folder> [<__smearing_config.txt>]
(without the second argument, $SMEARING_CONFIG is used if set)
reads density_per_layer_before_after_cuts.csv there, prints E[#fakes] and
the headroom (the common factor on all six densities that gives one
expected fake, mu_E1 = E^(-1/6)), writes fake_rate.csv and plots E[#fakes]
vs. that factor (fake_rate_vs_density_multiplier.png) next to it.
"""
import csv
import os
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
# Calibrated fake probability P_true for this tower, per track model,
# at the reference position resolution SIGMA_REF_MM in both coordinates
P_TRUE = {"exact_helix": 7.3054e-26, "conservative_quad": 6.2186e-26, "line_B0": 3.6895e-28}
SIGMA_REF_MM = 0.1
# (bending-view ndof, depth-view ndof) per model: P_true ~ sigma_u^nb sigma_v^nd
NDOF_VIEWS = {"exact_helix": (4, 4), "conservative_quad": (4, 4), "line_B0": (5, 4)}
MODEL_LABEL = {"exact_helix": "exact helix", "conservative_quad": "conservative quad",
               "line_B0": "line (B=0)"}


def p_true(sigma_u=SIGMA_REF_MM, sigma_v=SIGMA_REF_MM):
    """Fake probability per model at position resolutions sigma_u (r-phi)
    and sigma_v (z), in mm, scaled from the calibration at SIGMA_REF_MM."""
    return {m: p * (sigma_u / SIGMA_REF_MM) ** NDOF_VIEWS[m][0]
               * (sigma_v / SIGMA_REF_MM) ** NDOF_VIEWS[m][1] for m, p in P_TRUE.items()}


def position_sigmas(smear_cfg):
    """(sigma_u, sigma_v) in mm for the fake probability, from a loaded
    smearing config; SIGMA_REF_MM for both if position smearing is off."""
    pos = (smear_cfg or {}).get("position") or {}
    if pos.get("enabled") and pos.get("sigma_u", 0) > 0 and pos.get("sigma_v", 0) > 0:
        return float(pos["sigma_u"]), float(pos["sigma_v"])
    return SIGMA_REF_MM, SIGMA_REF_MM


def efakes(n_after, area=AREA_MM2, sigma_u=SIGMA_REF_MM, sigma_v=SIGMA_REF_MM):
    """n_after: {(system, layer): BIB hits after the cuts} for the 6 tower
    layers; sigma_u, sigma_v: hit position resolutions (mm).
    Returns (prod_n, {model: E[#fakes]})."""
    prod_n = 1.0
    for key in AREA_MM2:
        prod_n *= n_after[key] / area[key] * W2_SYNTH[key]
    return prod_n, {m: p * prod_n for m, p in p_true(sigma_u, sigma_v).items()}


N_LAYERS = len(AREA_MM2)


def mu_one_fake(e):
    """Common factor on all six densities that gives E[#fakes] = 1 (E scales as mu^6)."""
    return e ** (-1.0 / N_LAYERS) if e > 0 else float("inf")


def summary_line(e):
    return ("Expected fake tracks E[#fakes], 6-layer IT+OT barrel tower (exact helix): "
            f"{e['exact_helix']:.2e};  1 fake at {mu_one_fake(e['exact_helix']):.0f}x these densities")


def plot_headroom(e, out_png):
    """E[#fakes] vs. a common factor mu on all six layer densities (E ~ mu^6),
    for the exact-helix model (the other two models stay in fake_rate.csv),
    marking this run (mu = 1) and mu at E = 1."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    e0 = e["exact_helix"]
    mu1 = mu_one_fake(e0)
    mu = np.geomspace(0.1, max(3 * mu1, 10), 300)
    fig, ax = plt.subplots(figsize=(8, 5.2))
    c = "#1b9e77"
    ax.plot(mu, e0 * mu ** N_LAYERS, color=c, lw=2.2)
    ax.plot([mu1], [1.0], "o", color=c, ms=8, zorder=5)
    ax.annotate(f"1 fake at {mu1:.0f}x", (mu1, 1.0), xytext=(-12, 10), textcoords="offset points",
                ha="right", va="bottom", fontsize=11, weight="bold", color=c)
    ax.plot([1.0], [e0], "s", color="black", ms=6, zorder=5)
    ax.annotate(f"this run: {e0:.1e}", (1.0, e0), xytext=(8, -4), textcoords="offset points",
                ha="left", va="top", fontsize=9)
    ax.axhline(1.0, color="gray", lw=0.9)
    ax.axvline(1.0, color="black", lw=1.0, alpha=0.6)
    ax.text(1.06, 0.03, "this run's\nBIB densities", transform=ax.get_xaxis_transform(),
            fontsize=8, va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"common factor $\mu$ on the BIB hit density of all six layers")
    ax.set_ylabel(r"$E[\#\mathrm{fakes}]$ (6-layer IT+OT barrel tower)")
    ax.set_title(r"Fake-track headroom: $E[\#\mathrm{fakes}] \propto \mu^6$", fontsize=11)
    try:
        import bib_common
        bib_common.add_grid(fig)
    except Exception:
        ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main(argv):
    if len(argv) not in (1, 2):
        print(__doc__.split("Used by")[1], file=sys.stderr)
        return 2
    folder = Path(argv[0])
    smear_path = argv[1] if len(argv) == 2 else os.environ.get("SMEARING_CONFIG")
    smear_cfg = None
    if smear_path:
        import bib_common
        smear_cfg = bib_common.load_smearing_config(smear_path)
    else:
        print(f"NOTE: no smearing config given - fake probability at the reference "
              f"position resolution {SIGMA_REF_MM} mm", file=sys.stderr)
    su, sv = position_sigmas(smear_cfg)
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
    prod_n, e = efakes(n_after, area, su, sv)
    with open(folder / "fake_rate.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sigma_u_mm", "sigma_v_mm", "prod_n"] + [f"efakes_{m}" for m in P_TRUE]
                   + [f"mu_one_fake_{m}" for m in P_TRUE]
                   + [f"bib_hits_after_{'IT' if s == 3 else 'OT'}_barrel_L{l}" for s, l in AREA_MM2])
        w.writerow([su, sv, prod_n] + [e[m] for m in P_TRUE] + [mu_one_fake(e[m]) for m in P_TRUE]
                   + [n_after[k] for k in AREA_MM2])
    plot_headroom(e, folder / "fake_rate_vs_density_multiplier.png")
    print("BIB hits after the cuts, IT/OT barrel layers:  "
          + "  ".join(f"{'IT' if s == 3 else 'OT'} L{l} {n_after[(s, l)]:,}" for s, l in AREA_MM2))
    print(summary_line(e))
    print(f"  (density after cuts of each layer x tower plane area, times the fake probability "
          f"calibrated at {SIGMA_REF_MM} mm and scaled to the position resolution "
          f"{su:g}/{sv:g} mm (r-phi/z); Sections 14-16)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
