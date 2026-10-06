"""
fit_tail_fixed_exponent.py

Test whether fixing the power-law exponent to the theoretical value
ndof/2 = 6 (and fitting only the normalization K) is justified, by
computing the "local" implied K = P_empirical(cut) / cut^6 across the
available tail and checking whether it converges to a stable plateau as
cut gets smaller (deeper into the asymptotic regime) rather than still
trending -- a trend would mean either more statistics or a correction
term is needed before trusting a fixed-exponent extrapolation.
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

TOTAL_TRIALS = 100_000_000
CHUNK = 2_000_000
SEED = 777

NDOF = 2 * N_PLANES - 4
POWER = NDOF / 2.0  # = 6, fixed by theory


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
    print(f"Generated {TOTAL_TRIALS:,} trials in {time.time()-t0:.0f}s")

    sorted_chi2 = np.sort(all_chi2)
    n = len(sorted_chi2)
    del all_chi2

    # Look at implied K = P(cut)/cut^POWER at many depths in the tail,
    # from very few events (deep, noisy) to many events (shallow, precise
    # but possibly biased away from the asymptotic power law).
    event_counts = np.unique(np.geomspace(20, int(0.02 * n), 60).astype(int))
    cuts = sorted_chi2[event_counts]
    p_emp = event_counts / n
    K_local = p_emp / cuts ** POWER
    # Poisson-ish relative uncertainty on p_emp (and hence K_local): 1/sqrt(count)
    rel_err = 1.0 / np.sqrt(event_counts)

    print(f"{'n_events':>10} {'cut(chi2)':>12} {'P_emp':>12} {'K_local=P/cut^6':>18} {'rel.err':>10}")
    for ne, c, p, k, re in zip(event_counts, cuts, p_emp, K_local, rel_err):
        print(f"{ne:>10d} {c:>12.4f} {p:>12.3e} {k:>18.4e} {re*100:>9.1f}%")

    # Weighted mean of K_local over the deepest (smallest-cut) half of the
    # sampled range, weighted by inverse variance (~ event_counts), as a
    # robust plateau estimate.
    deep_half = cuts < np.median(cuts)
    weights = event_counts[deep_half].astype(float)
    K_plateau = np.sum(K_local[deep_half] * weights) / np.sum(weights)
    print(f"\nWeighted-mean K over deepest half of tail (cut < {np.median(cuts):.2f}): "
          f"K = {K_plateau:.4e}")

    # Compare to shallow half, to see how much K drifts across the range
    shallow_half = ~deep_half
    weights_s = event_counts[shallow_half].astype(float)
    K_shallow = np.sum(K_local[shallow_half] * weights_s) / np.sum(weights_s)
    print(f"Weighted-mean K over shallow half of tail (cut >= {np.median(cuts):.2f}): "
          f"K = {K_shallow:.4e}")
    print(f"Ratio (shallow/deep): {K_shallow/K_plateau:.3f}  "
          f"(close to 1.0 => exponent=6 is a good fit across this whole range;"
          f" far from 1.0 => still trending, need deeper tail or a correction term)")

    # ---- plot: K_local(cut) vs cut ----
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.errorbar(cuts, K_local, yerr=K_local * rel_err, fmt="o", ms=4,
                color="#4C72B0", ecolor="#4C72B0", alpha=0.7, capsize=2)
    ax.axhline(K_plateau, color="r", ls="--", lw=1.5,
               label=f"weighted mean (deep half) K={K_plateau:.3e}")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\chi^2$ cut")
    ax.set_ylabel(r"implied $K = P(\chi^2<\mathrm{cut}) / \mathrm{cut}^6$")
    ax.set_title("Checking the fixed-exponent (ndof/2=6) assumption:\n"
                 "is the implied normalization K flat vs. cut?")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig("K_local_vs_cut.png", dpi=150)
    print("\nSaved plot to K_local_vs_cut.png")


if __name__ == "__main__":
    main()
