"""
run_chi2_experiment.py

First deliverable: for N=8 planes (1000x1000 mm, sigma_x=sigma_y=10 micron),
repeatedly pick ONE random uniformly-distributed hit per plane (pure
combinatorial noise, no real track), fit a straight 3D line, and record
the resulting chi2. This characterizes the chi2 distribution of purely
random ("fake") hit combinations, which is the basis for later estimating
how many fake tracks to expect at a given chi2 cut, once a hit density
is specified.
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
SIZE_X = 1000.0     # mm
SIZE_Y = 1000.0     # mm
SIGMA_X = 0.010      # mm (10 micron)
SIGMA_Y = 0.010      # mm (10 micron)
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
    print(f"(theory for a REAL track: mean=ndof={NDOF}, var=2*ndof={2*NDOF})")
    print()
    print("--> Random-combination chi2 is ~8-9 orders of magnitude above the")
    print("    real-track scale (ndof=12), because sigma (10 micron) is tiny")
    print("    compared to the plane size (1000 mm) hits are scattered over.")
    print("    A chi2 cut anywhere near ndof will reject ~100% of fakes here.")
    print()
    print("Reference: fraction of random combos passing 'real-track-scale' cuts")
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
    np.save("chi2_values_N8.npy", chi2_values)

    # ---- plot ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    log_chi2 = np.log10(chi2_values)

    # left: histogram of log10(chi2) -- the natural scale here
    ax = axes[0]
    ax.hist(log_chi2, bins=100, color="#4C72B0", alpha=0.75)
    ax.axvline(np.log10(NDOF), color="r", ls="--", lw=2,
               label=f"real-track scale ($\\chi^2$=ndof={NDOF})")
    ax.set_xlabel(r"$\log_{10}(\chi^2)$")
    ax.set_ylabel("trials")
    ax.set_title(f"Random-combination $\\chi^2$ distribution\n(N={N_PLANES} planes, "
                 f"{SIZE_X:.0f}x{SIZE_Y:.0f} mm, $\\sigma$={SIGMA_X*1000:.0f} $\\mu$m)")
    ax.legend()

    # right: empirical survival function P(chi2 >= cut) vs cut, on log-log axes,
    # covering the full observed range -- this is what you read a fake rate off of
    ax2 = axes[1]
    sorted_chi2 = np.sort(chi2_values)
    n = len(sorted_chi2)
    survival = 1.0 - np.arange(n) / n  # P(chi2 >= sorted_chi2[i])
    ax2.loglog(sorted_chi2, survival, color="#4C72B0", lw=1.5)
    ax2.axvline(NDOF, color="r", ls="--", lw=2, label=f"real-track scale (ndof={NDOF})")
    ax2.set_xlabel(r"$\chi^2$ cut")
    ax2.set_ylabel(r"$P(\chi^2 \geq \mathrm{cut})$  (fraction of random combos passing)")
    ax2.set_title("Empirical survival function (log-log)")
    ax2.legend()
    ax2.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig("chi2_distribution_N8.png", dpi=150)
    print("\nSaved plot to chi2_distribution_N8.png")
    print("Saved raw chi2 values to chi2_values_N8.npy")


if __name__ == "__main__":
    main()
