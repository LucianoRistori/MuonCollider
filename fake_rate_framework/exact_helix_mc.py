# CORRECTION (10 October 2026): the depth-view fit in this framework had r and z
# swapped (it fitted r = p + q z, with residuals in r instead of z), so the fake
# probabilities computed or used here are about 450 times too low for the helix.
# Kept as the record of the earlier paper draft; the corrected values come from
# fake_rate_framework/geometric_K.py and are used by fake_rate.py.
"""
exact_helix_mc.py

The TRUE 4-parameter beam-origin-constrained helix fit (b, kappa, p, q),
replacing Section 10.4/11's "conservative" independent-(b,c) quadratic
stand-in with the actual physics: bending view fit via the exact circle
X(Y;b,kappa) (exact_helix_fit_core.py, vectorized multi-start LM), depth
view unchanged (r=p+qZ, linear, 2 params). This is the fake-rate
calibration for the model Luciano will actually use on real data.

Same geometry, same whole-geometry shrink (Section 11.3's fix), same
sigma as barrel_helix_tower.py, so results are directly comparable to
Section 11.5's table and the Section-11.5-corrected origin-constrained
straight-line baseline.
"""

import time
import numpy as np
from scipy import stats
from exact_helix_fit_core import fit_exact_bend

N_PLANES = 6
Y_SPACING = 200.0
SIGMA = 0.100

Y_true = np.array([Y_SPACING * i for i in range(1, N_PLANES + 1)])
W_true = np.array([1000.0 * i / N_PLANES for i in range(1, N_PLANES + 1)])

NDOF = 2 * N_PLANES - 4   # = 8 (same as Section 10.4: 2 params/view x 2 views)
POWER = NDOF / 2.0        # = 4

TARGET_RATIO_NEAREST = 10.0
LAMBDA = W_true[0] / (TARGET_RATIO_NEAREST * SIGMA)
W_shrunk = W_true / LAMBDA
Y_shrunk = Y_true / LAMBDA
inv_sigma2 = 1.0 / SIGMA ** 2

SEEDS5 = np.array([1e-6, -0.3, -0.12, 0.12, 0.3])  # scaled per-call to kappa_clip


def batch_fit_exact(n_trials, rng, W, Y, n_iter=22):
    X = rng.uniform(-0.5, 0.5, size=(n_trials, len(Y))) * W[None, :]
    Z = rng.uniform(-0.5, 0.5, size=(n_trials, len(Y))) * W[None, :]

    Y_max = np.max(np.abs(Y))
    kmax = 3.0 / Y_max
    kappa_seeds = kmax * np.array([1e-6, -0.6, -0.24, 0.24, 0.6]) / 0.3  # rescale template to this Y's kmax
    # (template SEEDS5 was tuned for Y_max=7.2mm, kmax~0.417; rescale proportionally)
    kappa_seeds = SEEDS5 * (kmax / (3.0 / 7.2))

    b_fit, k_fit, chi2_bend = fit_exact_bend(X, Y, inv_sigma2, n_iter=n_iter, kappa_seeds=kappa_seeds)

    # depth view: unchanged linear r = p + q*Z
    r = np.sqrt(X ** 2 + Y[None, :] ** 2)
    n_pts = len(Y)
    Sz = np.sum(Z, axis=1); Szz = np.sum(Z * Z, axis=1)
    Sr = np.sum(r, axis=1); Srz = np.sum(r * Z, axis=1)
    det = n_pts * Szz - Sz ** 2
    q = (n_pts * Srz - Sz * Sr) / det
    p = (Sr - q * Sz) / n_pts
    resid_depth = r - (p[:, None] + q[:, None] * Z)
    chi2_depth = np.sum(resid_depth ** 2, axis=1) * inv_sigma2

    return chi2_bend + chi2_depth, chi2_bend, chi2_depth, b_fit, k_fit


if __name__ == "__main__":
    print(f"ndof={NDOF}, power={POWER}, LAMBDA={LAMBDA:.4f}")

    # --- Step 0: scaling law check ---
    print("\n--- Step 0: scaling law check (exact 4-param model) ---")
    for lam_test in [LAMBDA, LAMBDA * 7.0]:
        Wt = W_true / lam_test
        Yt = Y_true / lam_test
        rng0 = np.random.default_rng(1)
        chi2_t, *_ = batch_fit_exact(3000, rng0, Wt, Yt)
        chi2_sorted = np.sort(chi2_t)
        cut1pct = chi2_sorted[30]  # ~1st percentile of 3000
        s = cut1pct * lam_test ** 2
        print(f"  lambda={lam_test:9.2f}  chi2[1pct]*lambda^2 = {s:.6f}")

    # --- Step 1: deep-tail power-law plateau (LEFT tail, same convention as barrel_helix_tower.py) ---
    print("\n--- Step 1: chi2_total power-law plateau (power=4), exact 4-param model ---")
    N_TRIALS = 4_000_000
    BATCH = 200_000
    SEED = 20261002
    rng = np.random.default_rng(SEED)
    all_chi2 = np.empty(N_TRIALS, dtype=np.float64)
    all_kappa = np.empty(N_TRIALS, dtype=np.float64)
    t0 = time.time()
    done = 0
    while done < N_TRIALS:
        nb = min(BATCH, N_TRIALS - done)
        c2, cb, cd, bf, kf = batch_fit_exact(nb, rng, W_shrunk, Y_shrunk)
        all_chi2[done:done+nb] = c2
        all_kappa[done:done+nb] = kf
        done += nb
        if done % 1_000_000 < BATCH:
            print(f"  ... {done:,}/{N_TRIALS:,} trials, {time.time()-t0:.0f}s elapsed")
    print(f"Generated {N_TRIALS:,} trials in {time.time()-t0:.0f}s")

    order = np.argsort(all_chi2)
    sorted_chi2 = all_chi2[order]
    sorted_kappa = all_kappa[order]
    n = len(sorted_chi2)

    event_counts = np.unique(np.geomspace(20, int(0.02 * n), 50).astype(int))
    cuts = sorted_chi2[event_counts]
    p_emp = event_counts / n
    K_local = p_emp / cuts ** POWER
    rel_err = 1.0 / np.sqrt(event_counts)

    print(f"\n{'n_events':>10} {'cut(chi2)':>12} {'P_emp':>12} {'K_local':>14} {'rel.err':>8}")
    for ne, cc, pp, kk, re in zip(event_counts, cuts, p_emp, K_local, rel_err):
        print(f"{ne:>10d} {cc:>12.6f} {pp:>12.3e} {kk:>14.4e} {re*100:>7.1f}%")

    deep = cuts < np.percentile(cuts, 25)
    w_ = event_counts[deep].astype(float)
    K_EXACT = float(np.sum(K_local[deep] * w_) / np.sum(w_))
    print(f"\nPlateau K_shrunk (deep quartile) = {K_EXACT:.4e}")

    np.save('/tmp/exact_helix_sorted_chi2.npy', sorted_chi2)
    np.save('/tmp/exact_helix_sorted_kappa.npy', sorted_kappa)

    # --- Step 2: kappa_max acceptance, characterized directly (no proxy needed) ---
    print("\n--- Step 2: kappa_max cut suppression f(kappa_max), exact kappa (no proxy) ---")
    kappa_abs = np.abs(all_kappa)
    pct_vals = [0.5, 1, 2, 5, 10, 20]
    kappa_pcts = np.percentile(kappa_abs, pct_vals)
    print("kappa_abs (shrunk) percentiles:", {p: f"{v:.4e}" for p, v in zip(pct_vals, kappa_pcts)})

    chi2_test_cuts_n = [200, 500, 1200, 3000]
    print(f"\n{'kappa_max':>12} {'chi2_cut':>10} {'P(joint)':>12} {'P(chi2 only)':>13} {'f=ratio':>9}")
    f_results = {}
    for kmax_s in kappa_pcts:
        fs = []
        for ne in chi2_test_cuts_n:
            cut = sorted_chi2[ne]
            mask_joint = (sorted_chi2 <= cut) & (np.abs(sorted_kappa) <= kmax_s)
            p_joint = np.sum(mask_joint) / n
            p_chi2only = ne / n
            f = p_joint / p_chi2only if p_chi2only > 0 else np.nan
            fs.append(f)
            print(f"{kmax_s:12.4e} {cut:10.3f} {p_joint:12.3e} {p_chi2only:13.3e} {f:9.4f}")
        f_results[kmax_s] = np.mean(fs)
        print(f"  --> mean f over chi2 cuts tested: {np.mean(fs):.4f}  (spread {np.std(fs):.4f})")

    print(f"\n{'kappa_max':>12} {'f(kmax)':>10} {'f/kmax':>12}")
    ratios = []
    for kmax_s, f in f_results.items():
        ratio = f / kmax_s
        ratios.append(ratio)
        print(f"{kmax_s:12.4e} {f:10.4f} {ratio:12.4e}")
    C_KAPPA_EXACT = ratios[0]
    print(f"\nC_kappa (linear-law, smallest kmax) = {C_KAPPA_EXACT:.4e}")

    # --- Step 3: combine to true-detector numbers for p_T>10 GeV/c, B=5T ---
    print("\n--- Step 3: combine to true-detector numbers ---")
    B_TESLA = 5.0
    PT_MIN_GEV = 10.0
    R_MIN_MM = (PT_MIN_GEV / (0.2998 * B_TESLA)) * 1000.0
    KAPPA_MAX_TRUE = 1.0 / R_MIN_MM
    KAPPA_MAX_SHRUNK = KAPPA_MAX_TRUE * LAMBDA
    F_KAPPA = C_KAPPA_EXACT * KAPPA_MAX_SHRUNK

    P_3SIGMA = stats.norm.sf(3)
    CUT_3SIGMA = stats.chi2.isf(P_3SIGMA, NDOF)
    print(f"kappa_max (shrunk) = {KAPPA_MAX_SHRUNK:.4e},  F_kappa = {F_KAPPA:.4e}")
    print(f"ndof={NDOF}, matched cut = {CUT_3SIGMA:.4f}")

    p_true = F_KAPPA * K_EXACT * (CUT_3SIGMA / LAMBDA ** 2) ** POWER
    n1 = (1.0 / p_true) ** (1.0 / N_PLANES)
    n10 = (10.0 / p_true) ** (1.0 / N_PLANES)
    print(f"\nP_true(exact 4-param helix, p_T>10GeV, chi2<{CUT_3SIGMA:.2f}) = {p_true:.4e}")
    print(f"n(E=1)  = {n1:.0f} hits/layer")
    print(f"n(E=10) = {n10:.0f} hits/layer")

    print(f"\nK_EXACT = {K_EXACT:.6e}")
    print(f"C_KAPPA_EXACT = {C_KAPPA_EXACT:.6e}")
