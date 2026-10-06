"""
fit_tail_powerlaw_v2.py

High-statistics version: uses the vectorized batch_random_chi2 to generate
very large numbers of trials (hundreds of millions) in chunks, then fits
the power law P(chi2<cut) ~ K*cut^(ndof/2) using a lower, more asymptotic
cut range than was reachable with the 200k/2M-trial runs. This checks
whether the fitted exponent converges to the theoretical ndof/2 = 6 once
we're deep enough in the tail, and gives a more trustworthy extrapolation
to very small cut values.
"""

import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tracksim import Detector, batch_random_chi2

N_PLANES = 8
SIZE_X = 1.0
SIZE_Y = 1.0
SIGMA_X = 0.100
SIGMA_Y = 0.100
Z_SPACING = 100.0

TOTAL_TRIALS = 200_000_000
CHUNK = 5_000_000
SEED = 2024

NDOF = 2 * N_PLANES - 4
POWER = NDOF / 2.0  # = 6


def main():
    rng = np.random.default_rng(SEED)
    detector = Detector.make_uniform(
        n_planes=N_PLANES, size_x=SIZE_X, size_y=SIZE_Y,
        sigma_x=SIGMA_X, sigma_y=SIGMA_Y, z_spacing=Z_SPACING,
    )

    n_chunks = TOTAL_TRIALS // CHUNK
    all_chi2 = np.empty(TOTAL_TRIALS, dtype=np.float64)

    t0 = time.time()
    for i in range(n_chunks):
        all_chi2[i * CHUNK:(i + 1) * CHUNK] = batch_random_chi2(detector, CHUNK, rng)
        if (i + 1) % 8 == 0 or i == n_chunks - 1:
            elapsed = time.time() - t0
            print(f"  chunk {i+1}/{n_chunks}  ({(i+1)*CHUNK:,} trials)  "
                  f"elapsed={elapsed:.0f}s")

    print(f"Total time: {time.time()-t0:.0f}s for {TOTAL_TRIALS:,} trials")

    sorted_chi2 = np.sort(all_chi2)
    n = len(sorted_chi2)
    print(f"min chi2 observed: {sorted_chi2[0]:.4f}")
    print(f"chi2 at 1e-6 level (index {int(1e-6*n)}): {sorted_chi2[int(1e-6*n)] if int(1e-6*n)>0 else sorted_chi2[0]:.4f}")

    # Fit range: pick a window deep in the tail but with enough events for a
    # stable fit (require at least ~200 events at the low end of the window).
    fit_lo_idx = 200
    fit_hi_idx = int(1e-4 * n)  # up to ~0.01 percentile
    idx = np.unique(np.geomspace(fit_lo_idx, fit_hi_idx, 50).astype(int))

    x_fit = sorted_chi2[idx]
    y_fit = idx / n

    log_x = np.log(x_fit)
    log_y = np.log(y_fit)
    A = np.vstack([log_x, np.ones_like(log_x)]).T
    power_fit, log_K_fit = np.linalg.lstsq(A, log_y, rcond=None)[0]
    K_fit = np.exp(log_K_fit)

    log_K_fixed = np.mean(log_y - POWER * log_x)
    K_fixed = np.exp(log_K_fixed)

    print()
    print(f"Fit range: chi2 in [{x_fit.min():.3f}, {x_fit.max():.3f}], "
          f"P in [{y_fit.min():.2e}, {y_fit.max():.2e}]  ({len(idx)} points)")
    print(f"Free power-law fit:  P(chi2<cut) = {K_fit:.4e} * cut^{power_fit:.3f}")
    print(f"Fixed-exponent fit:  P(chi2<cut) = {K_fixed:.4e} * cut^{POWER:.1f}")
    print()
    print("Extrapolated P(chi2 < cut), fixed-exponent power law:")
    for cut in [12, 10, 8, 6, 4, 2, 1, 0.5, 0.1]:
        p_extrap = K_fixed * cut ** POWER
        print(f"  cut = {cut:>6.2f}  ->  P = {p_extrap:.3e}")

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    show_idx = np.unique(np.geomspace(1, n - 1, 600).astype(int))
    ax.loglog(sorted_chi2[show_idx], show_idx / n, ".", ms=2.5, color="#4C72B0",
              alpha=0.4, label=f"empirical (N={TOTAL_TRIALS:,} trials)")
    ax.loglog(x_fit, y_fit, "o", ms=5, color="#DD8452", label="fit window")

    xx = np.geomspace(0.05, sorted_chi2.max(), 300)
    ax.loglog(xx, K_fixed * xx ** POWER, "r-", lw=2,
              label=f"fixed exponent: $K \\cdot \\chi^{{2\\,{POWER:.0f}}}$")
    ax.loglog(xx, K_fit * xx ** power_fit, "g--", lw=1.5,
              label=f"free fit (exponent={power_fit:.2f})")
    ax.axvline(NDOF, color="k", ls=":", lw=1, alpha=0.6, label=f"ndof={NDOF}")

    ax.set_xlabel(r"$\chi^2$ cut")
    ax.set_ylabel(r"$P(\chi^2 < \mathrm{cut})$")
    ax.set_title(f"Low-$\\chi^2$ tail power law (high statistics)\n"
                 f"N={N_PLANES} planes, {SIZE_X:.1f}x{SIZE_Y:.1f} mm, "
                 f"$\\sigma$={SIGMA_X*1000:.0f} $\\mu$m")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig("chi2_tail_powerlaw_fit_v2.png", dpi=150)
    print("\nSaved plot to chi2_tail_powerlaw_fit_v2.png")


if __name__ == "__main__":
    main()
