"""
fake_rate.py - expected number of fake tracks, E[#fakes], in the 6-layer
IT + OT barrel tower of the fake-rate framework, from the BIB hit density of
each IT/OT barrel layer after the cuts.

A track counts as found if at least k = min_hits_found of its hits pass the
cuts (__cuts_config.txt, [track]; default 5). The matching fake definition is
"at least k of the 6 layers": by linearity of expectation

    E[#fakes] = sum over the sets S of >= k layers of  P_S x prod_{i in S} n_i,
    n_i = (BIB hits after cuts / area)_i x W_i^2,

with W_i^2 the plane areas of the projective tower (1 m^2 outermost plane) and
P_S the probability that one random combination of hits on the layers of S
passes the exact-helix fit (one-sided 3-sigma chi2 cut for its ndof, p_T > 10
GeV/c in 5 T). For k = 5: the 6-of-6 term plus six 5-of-6 terms. The sum
counts a fake that also passes on a subset more than once; for the densities
of interest the overlap is negligible (window occupancy << 1).

P_S is computed from the geometry (fake_rate_framework/geometric_K.py: volume
of the accepted track manifold, validated by direct Monte Carlo) at
sigma_ref = 0.1 mm in both coordinates, and scales with the position
resolution as
    P_S(sigma_u, sigma_v) = P_S,ref x (sigma_u/sigma_ref)^(m-2) x (sigma_v/sigma_ref)^(m-2)
for a set of m layers (m - 2 degrees of freedom per view). In the barrel, u
is the r-phi (bending) coordinate and v the z (depth) coordinate. The
resolutions come from __smearing_config.txt ([position] sigma_u_mm,
sigma_v_mm); with position smearing off, sigma_ref is used.

History: up to October 2026 P_true came from a Monte Carlo calibration
(fake_rate_framework/real_tower_calibration.py) whose depth-view fit had r and
z swapped (residuals in r instead of z); it was about 450 times too low, and
only the 6-of-6 term was counted.

Used by resolution_scan.py (./run_scan), and by ./run_all at the end of
its cuts step:
    python3 fake_rate.py <step4_cuts folder> [<__smearing_config.txt> [<__cuts_config.txt>]]
(without them, $SMEARING_CONFIG / $CUTS_CONFIG are used if set; with no cuts
config, k = 5) reads density_per_layer_before_after_cuts.csv there, prints
E[#fakes] and the headroom (the common factor mu on all six densities that
gives one expected fake), writes fake_rate.csv and plots E[#fakes] vs. that
factor (fake_rate_vs_density_multiplier.png) next to it.
"""
import csv
import itertools
import math
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
LAYER_KEYS = list(AREA_MM2)               # tower layers, inside out (index 0..5)
N_LAYERS = len(LAYER_KEYS)
# Plane areas W_i^2 (mm^2) of the projective tower (1 m^2 outermost plane)
Y_TRUE = {(IT_BARREL, 0): 164.0, (IT_BARREL, 1): 354.0, (IT_BARREL, 2): 554.0,
          (OT_BARREL, 0): 819.0, (OT_BARREL, 1): 1153.0, (OT_BARREL, 2): 1486.0}
W2_SYNTH = {k: (1000.0 * y / 1486.0) ** 2 for k, y in Y_TRUE.items()}

SIGMA_REF_MM = 0.1
# P_S at SIGMA_REF_MM in both coordinates, exact helix, for every set S of
# >= 4 tower layers (indices into LAYER_KEYS); from fake_rate_framework/geometric_K.py
P_REF_SETS = {
    (0, 1, 2, 3, 4, 5): 3.33483e-23,   # ndof 8, K 8.06119e-29
    (0, 1, 2, 3, 4): 2.97627e-17,   # ndof 6, K 2.89694e-21
    (0, 1, 2, 3, 5): 4.49772e-17,   # ndof 6, K 4.37784e-21
    (0, 1, 2, 4, 5): 2.39827e-17,   # ndof 6, K 2.33435e-21
    (0, 1, 3, 4, 5): 1.09356e-17,   # ndof 6, K 1.06442e-21
    (0, 2, 3, 4, 5): 4.50751e-18,   # ndof 6, K 4.38737e-22
    (1, 2, 3, 4, 5): 2.36728e-18,   # ndof 6, K 2.30418e-22
    (0, 1, 2, 3): 1.02458e-11,   # ndof 4, K 3.23355e-14
    (0, 1, 2, 4): 1.86523e-11,   # ndof 4, K 5.88660e-14
    (0, 1, 2, 5): 2.75282e-11,   # ndof 4, K 8.68781e-14
    (0, 1, 3, 4): 8.55362e-12,   # ndof 4, K 2.69949e-14
    (0, 1, 3, 5): 1.31998e-11,   # ndof 4, K 4.16581e-14
    (0, 1, 4, 5): 6.43619e-12,   # ndof 4, K 2.03124e-14
    (0, 2, 3, 4): 3.48263e-12,   # ndof 4, K 1.09911e-14
    (0, 2, 3, 5): 5.58671e-12,   # ndof 4, K 1.76315e-14
    (0, 2, 4, 5): 2.80722e-12,   # ndof 4, K 8.85949e-15
    (0, 3, 4, 5): 1.23921e-12,   # ndof 4, K 3.91089e-15
    (1, 2, 3, 4): 1.88552e-12,   # ndof 4, K 5.95064e-15
    (1, 2, 3, 5): 2.96479e-12,   # ndof 4, K 9.35676e-15
    (1, 2, 4, 5): 1.53319e-12,   # ndof 4, K 4.83868e-15
    (1, 3, 4, 5): 6.55615e-13,   # ndof 4, K 2.06910e-15
    (2, 3, 4, 5): 4.67907e-13,   # ndof 4, K 1.47670e-15
}
K_DEFAULT = 5
MODEL_LABEL = {"exact_helix": "exact helix, >= k of 6 layers", "exact_helix_6of6": "exact helix, 6 of 6"}
# kept for scripts that iterate over the models: the CSV columns efakes_<model>
P_TRUE = {m: None for m in MODEL_LABEL}


def p_sets(sigma_u=SIGMA_REF_MM, sigma_v=SIGMA_REF_MM):
    """P_S per layer set at position resolutions sigma_u (r-phi), sigma_v (z), in mm."""
    return {S: p * (sigma_u / SIGMA_REF_MM) ** (len(S) - 2) * (sigma_v / SIGMA_REF_MM) ** (len(S) - 2)
            for S, p in P_REF_SETS.items()}


def position_sigmas(smear_cfg):
    """(sigma_u, sigma_v) in mm for the fake probability, from a loaded
    smearing config; SIGMA_REF_MM for both if position smearing is off."""
    pos = (smear_cfg or {}).get("position") or {}
    if pos.get("enabled") and pos.get("sigma_u", 0) > 0 and pos.get("sigma_v", 0) > 0:
        return float(pos["sigma_u"]), float(pos["sigma_v"])
    return SIGMA_REF_MM, SIGMA_REF_MM


def tower_n(n_after, area=AREA_MM2):
    """Hits per tower plane: density after cuts x plane area, per layer (inside out)."""
    return [n_after[k] / area[k] * W2_SYNTH[k] for k in LAYER_KEYS]


def efake_terms(n_after, area=AREA_MM2, sigma_u=SIGMA_REF_MM, sigma_v=SIGMA_REF_MM, k_min=K_DEFAULT):
    """{m: sum over the sets of exactly m layers of P_S prod n_i}, for m = 6 .. k_min."""
    n = tower_n(n_after, area)
    terms = {}
    for S, p in p_sets(sigma_u, sigma_v).items():
        if len(S) >= k_min:
            terms[len(S)] = terms.get(len(S), 0.0) + p * math.prod(n[i] for i in S)
    return terms


def efakes(n_after, area=AREA_MM2, sigma_u=SIGMA_REF_MM, sigma_v=SIGMA_REF_MM, k_min=K_DEFAULT):
    """n_after: {(system, layer): BIB hits after the cuts} for the 6 tower
    layers; sigma_u, sigma_v: hit position resolutions (mm); k_min: hits needed.
    Returns (prod_n over the 6 layers, {model: E[#fakes]}): "exact_helix" is
    at least k_min of 6, "exact_helix_6of6" the 6-of-6 term alone."""
    terms = efake_terms(n_after, area, sigma_u, sigma_v, k_min)
    prod_n = math.prod(tower_n(n_after, area))
    return prod_n, {"exact_helix": sum(terms.values()), "exact_helix_6of6": terms.get(N_LAYERS, 0.0)}


def mu_one_fake(e):
    """Common factor mu on all six densities that gives one expected fake.
    e: {m: E_m} (terms scale as mu^m), or a number (pure 6-of-6 term)."""
    terms = e if isinstance(e, dict) else {N_LAYERS: e}
    terms = {m: v for m, v in terms.items() if v > 0}
    if not terms:
        return float("inf")
    f = lambda lm: math.log(sum(v * math.exp(m * lm) for m, v in terms.items()))
    lo, hi = -200.0, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return math.exp(0.5 * (lo + hi))


def k_label(k):
    return "6 of 6 layers" if k >= N_LAYERS else f"at least {k} of 6 layers"


def summary_line(e, terms, k):
    return (f"Expected fake tracks E[#fakes], 6-layer IT+OT barrel tower (exact helix, {k_label(k)}): "
            f"{e['exact_helix']:.2e};  1 fake at {mu_one_fake(terms):.1f}x these densities"
            + (f"  (6 of 6 alone: {e['exact_helix_6of6']:.1e})" if k < N_LAYERS else ""))


def plot_headroom(terms, k, out_png):
    """E[#fakes] vs. a common factor mu on all six layer densities (each
    m-layer term grows as mu^m), marking this run (mu = 1) and mu at E = 1."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    e0 = sum(terms.values())
    mu1 = mu_one_fake(terms)
    mu = np.geomspace(min(0.1, mu1 / 30), max(3 * mu1, 10), 300)
    tot = sum(v * mu ** m for m, v in terms.items())
    fig, ax = plt.subplots(figsize=(8, 5.2))
    c = "#1b9e77"
    ax.plot(mu, tot, color=c, lw=2.2, label=f"exact helix, {k_label(k)}")
    if k < N_LAYERS and terms.get(N_LAYERS, 0) > 0:
        ax.plot(mu, terms[N_LAYERS] * mu ** N_LAYERS, color=c, lw=1.2, ls="--", label="6 of 6 alone")
        ax.legend(fontsize=9, loc="lower right")
    ax.plot([mu1], [1.0], "o", color=c, ms=8, zorder=5)
    ax.annotate(f"1 fake at {mu1:.1f}x", (mu1, 1.0), xytext=(-12, 10), textcoords="offset points",
                ha="right", va="bottom", fontsize=11, weight="bold", color=c)
    ax.plot([1.0], [e0], "s", color="black", ms=6, zorder=5)
    ax.annotate(f"this run: {e0:.1e}", (1.0, e0), xytext=(8, -4), textcoords="offset points",
                ha="left", va="top", fontsize=9)
    ax.axhline(1.0, color="gray", lw=0.9)
    ax.axvline(1.0, color="black", lw=1.0, alpha=0.6)
    ax.text(1.06, 0.03, "this run's\nBIB densities", transform=ax.get_xaxis_transform(),
            fontsize=8, va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(min(tot.min(), e0) / 3, tot.max() * 3)
    ax.set_xlabel(r"common factor $\mu$ on the BIB hit density of all six layers")
    ax.set_ylabel(r"$E[\#\mathrm{fakes}]$ (6-layer IT+OT barrel tower)")
    pw = r"\mu^6" if k >= N_LAYERS else " + ".join(r"a_{%d}\mu^{%d}" % (m, m) for m in range(N_LAYERS, k - 1, -1))
    ax.set_title(r"Fake-track headroom: $E[\#\mathrm{fakes}] \propto %s$" % pw, fontsize=11)
    try:
        import bib_common
        bib_common.add_grid(fig)
    except Exception:
        ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def min_hits_from(cuts_path):
    if not cuts_path:
        return K_DEFAULT
    import bib_common
    return int(bib_common.load_track_params(cuts_path)["min_hits_found"])


def clamp_k(k):
    if k > N_LAYERS:
        return N_LAYERS
    if k < 4:
        print(f"WARNING: min_hits_found = {k}: fake probabilities exist only for >= 4 of the 6 "
              f"tower layers; using 4", file=sys.stderr)
        return 4
    return k


def main(argv):
    if len(argv) not in (1, 2, 3):
        print(__doc__.split("Used by")[1], file=sys.stderr)
        return 2
    folder = Path(argv[0])
    smear_path = argv[1] if len(argv) >= 2 else os.environ.get("SMEARING_CONFIG")
    cuts_path = argv[2] if len(argv) == 3 else os.environ.get("CUTS_CONFIG")
    smear_cfg = None
    if smear_path:
        import bib_common
        smear_cfg = bib_common.load_smearing_config(smear_path)
    else:
        print(f"NOTE: no smearing config given - fake probability at the reference "
              f"position resolution {SIGMA_REF_MM} mm", file=sys.stderr)
    if not cuts_path:
        print(f"NOTE: no cuts config given - fakes counted for at least {K_DEFAULT} of 6 layers",
              file=sys.stderr)
    k = clamp_k(min_hits_from(cuts_path))
    su, sv = position_sigmas(smear_cfg)
    n_after, area = {}, {}
    with open(folder / "density_per_layer_before_after_cuts.csv") as f:
        for r in csv.DictReader(f):
            key = (int(r["system"]), int(r["layer"]))
            if key in AREA_MM2:
                n_after[key] = int(r["bib_n_hits_after"])
                area[key] = float(r["area_mm2"])
    missing = [x for x in AREA_MM2 if x not in n_after]
    if missing:
        print(f"ERROR: no IT/OT barrel layer(s) {missing} in the density table", file=sys.stderr)
        return 1
    for x in AREA_MM2:      # the tower assumes these exact layer areas
        if abs(area[x] / AREA_MM2[x] - 1) > 1e-3:
            print(f"WARNING: layer {x} area {area[x]:.0f} mm^2 differs from the tower's "
                  f"{AREA_MM2[x]:.0f} mm^2 (different geometry?)", file=sys.stderr)
    terms = efake_terms(n_after, area, su, sv, k)
    prod_n, e = efakes(n_after, area, su, sv, k)
    mu_k, mu_6 = mu_one_fake(terms), mu_one_fake(e["exact_helix_6of6"])
    with open(folder / "fake_rate.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sigma_u_mm", "sigma_v_mm", "k_min", "prod_n", "efakes_exact_helix",
                    "efakes_exact_helix_6of6"] + [f"efakes_{m}_layers" for m in sorted(terms, reverse=True)]
                   + ["mu_one_fake_exact_helix", "mu_one_fake_exact_helix_6of6"]
                   + [f"bib_hits_after_{'IT' if s == 3 else 'OT'}_barrel_L{l}" for s, l in AREA_MM2])
        w.writerow([su, sv, k, prod_n, e["exact_helix"], e["exact_helix_6of6"]]
                   + [terms[m] for m in sorted(terms, reverse=True)] + [mu_k, mu_6]
                   + [n_after[x] for x in AREA_MM2])
    plot_headroom(terms, k, folder / "fake_rate_vs_density_multiplier.png")
    print("BIB hits after the cuts, IT/OT barrel layers:  "
          + "  ".join(f"{'IT' if s == 3 else 'OT'} L{l} {n_after[(s, l)]:,}" for s, l in AREA_MM2))
    print(summary_line(e, terms, k))
    print(f"  (density after cuts of each layer x tower plane area, times the fake probability of "
          f"each layer set from the geometry, at the position resolution {su:g}/{sv:g} mm (r-phi/z))")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
