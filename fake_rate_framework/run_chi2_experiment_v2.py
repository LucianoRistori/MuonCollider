"""
run_chi2_experiment_v2.py

Variant of the first experiment, but with much smaller planes (1x1 mm)
and coarser resolution (sigma=100 micron), so the plane is only ~10 sigma
across. This makes the "interesting" chi2 region (near ndof) directly
reachable by plain brute-force Monte Carlo, unlike the 1000x1000 mm /
10 micron case where it was ~1e9 sigma-equivalents away.
"""

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tracksim import Detector, one_random_combination_chi2

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
N_PLANES = 8
SIZE_X = 1.0         # mm
SIZE_Y = 1.0         # mm
SIGMA_X = 0.100      # mm (100 micron)
SIGMA_Y = 0.100      # mm (100 micron)
Z_SPACING = 100.0    # mm between consecutive planes (hardwired choice)

N_TRIALS = 200_000
SEED = 12345

NDOF = 2 * N_PLANES - 4  # = 12 for N=8

# Cuts near the "real track" scale (chi2 ~ ndof), for reference only --
# with sigma=10 micron and 1000mm planes we expect ~0 random combos to
# pass these; the interesting cuts are found empirically below.
CHI2_CUTS_SMALL = [10, 15, 20, 25, 30, 40, 50, NDOF * 2, NDOF * 3]

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def main():
    rng = np.random.default_rng(SEED)
    detector = Detector.make_uniform(
        n_planes=N_PLANES,
        size_x=SIZE_X,
        size_y=SIZE_Y,
        sigma_x=SIGMA_X,
        sigma_y=SIGMA_Y,
        z_spacing=Z_SPACING,
    )

    chi2_values = np.empty(N_TRIALS)
    for i in range(N_TRIALS):
        chi2_values[i] = one_random_combination_chi2(detector, rng)

    # ---- summary stats ----
    print(f"N_planes = {N_PLANES}, ndof = {NDOF}")
    print(f"N_trials = {N_TRIALS}")
    print(f"Plane size = {SIZE_X} x {SIZE_Y} mm, sigma = {SIGMA_X*1000:.1f} micron")
    print()
    print(f"Empirical mean(chi2)   = {chi2_values.mean():.4g}")
    print(f"Empirical median(chi2) = {np.median(chi2_values):.4g}")
    print(f"Empirical min/max      = {chi2_values.min():.4g} / {chi2_values.max():.4g}")
    print(f"(theory for a REAL (Gaussian-scattered) track: mean=ndof={NDOF}, var=2*ndof={2*NDOF})")
    print("(Note: random UNIFORM-noise combinations are NOT expected to follow the")
    print(" chi2(ndof) distribution exactly -- that theoretical curve describes Gaussian")
    print(" measurement smearing around a true line, not coincidental alignment of")
    print(" uniformly-scattered noise. It's shown only as a reference scale.)")
    print()
    print("Fraction of random combos passing various chi2 cuts")
    print(f"{'chi2 cut':>10} | {'P(chi2 < cut) empirical':>24} | {'theoretical CDF (real trk)':>26}")
    for cut in CHI2_CUTS_SMALL:
        p_emp = np.mean(chi2_values < cut)
        p_theory = stats.chi2.cdf(cut, df=NDOF)
        print(f"{cut:>10.1f} | {p_emp:>24.5f} | {p_theory:>26.5f}")

    print()
    print("Empirical percentiles of the random-combination chi2 distribution:")
    for pct in [0.01, 0.1, 1, 5, 25, 50, 75, 95, 99, 99.9]:
        val = np.percentile(chi2_values, pct)
        print(f"  {pct:>6.2f}th percentile: chi2 = {val:.4g}")

    # save raw values + summary for downstream use (e.g. combining with hit density)
    np.save("chi2_values_N8_v2.npy", chi2_values)

    # ---- plot ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # left: histogram (linear scale, now directly comparable to chi2(ndof))
    ax = axes[0]
    hi = np.percentile(chi2_values, 99.5)
    bins = np.linspace(0, max(hi, NDOF * 3), 120)
    ax.hist(chi2_values, bins=bins, density=True, alpha=0.6,
            label=f"empirical (N={N_TRIALS})", color="#4C72B0")
    xx = np.linspace(bins[0], bins[-1], 500)
    ax.plot(xx, stats.chi2.pdf(xx, df=NDOF), "r-", lw=2,
            label=f"$\\chi^2$ pdf, ndof={NDOF} (Gaussian-track reference)")
    ax.set_xlabel(r"$\chi^2$")
    ax.set_ylabel("probability density")
    ax.set_title(f"Random-combination $\\chi^2$ distribution\n(N={N_PLANES} planes, "
                 f"{SIZE_X:.1f}x{SIZE_Y:.1f} mm, $\\sigma$={SIGMA_X*1000:.0f} $\\mu$m)")
    ax.legend()

    # right: empirical survival function P(chi2 >= cut) vs cut, on log-log axes
    ax2 = axes[1]
    sorted_chi2 = np.sort(chi2_values)
    n = len(sorted_chi2)
    survival = 1.0 - np.arange(n) / n  # P(chi2 >= sorted_chi2[i])
    ax2.loglog(sorted_chi2, survival, color="#4C72B0", lw=1.5, label="empirical")
    ax2.loglog(xx[xx > 0], stats.chi2.sf(xx[xx > 0], df=NDOF), "r--", lw=2,
               label=f"$\\chi^2$ sf, ndof={NDOF} (Gaussian-track reference)")
    ax2.set_xlabel(r"$\chi^2$ cut")
    ax2.set_ylabel(r"$P(\chi^2 \geq \mathrm{cut})$  (fraction of random combos passing)")
    ax2.set_title("Empirical survival function (log-log)")
    ax2.legend()
    ax2.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig("chi2_distribution_N8_v2.png", dpi=150)
    print("\nSaved plot to chi2_distribution_N8_v2.png")
    print("Saved raw chi2 values to chi2_values_N8_v2.npy")


if __name__ == "__main__":
    main()
