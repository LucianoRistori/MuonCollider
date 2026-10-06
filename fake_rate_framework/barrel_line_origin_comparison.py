"""
barrel_line_origin_comparison.py

Correct control for Section 11's B-field comparison.

Luciano's question: the Section 11.5 "straight-line baseline" (K_LINE =
1.122e-11, ndof=8) was pulled from the OLD generic 3D line model
(x(z)=a+b*z, y(z)=c+d*z, 4 free parameters, NO origin constraint). That
is not the right B=0 control for the beam-origin-constrained curved-track
model of Section 10/11 -- the origin constraint (X0=Y0=0) by itself
already removes two degrees of freedom and tightens the fit, independent
of B.

The correct control is the B=0 LIMIT of the Section 10 model itself: same
beam-origin-constrained geometry, curvature term c dropped identically --
  bending view:  X(Y) = b*Y            (no intercept, NO curvature; 1 free param)
  depth view:    r(Z) = p + q*Z        (2 free params)
ndof = (N-1) + (N-2) = 2N-3 = 9  (one MORE dof than the curved model's
2N-4=8, since the curved model spends one extra parameter -- c -- fitting
the same data, so naive chi2-only rejection should be WEAKER for the
curved model, not stronger; the kappa_max cut is what has to do the work).

Same geometry, same shrink (LAMBDA=166.6667), same sigma, same N=6.
No kappa_max cut here (no curvature parameter to cut on).
Tail convention (matches barrel_helix_tower.py exactly): fakes are the
LEFT tail -- P(chi2 < cut), noise hits accidentally giving a GOOD fit --
so cut = sorted_chi2[n_events] with chi2 sorted ASCENDING, not the right
tail.
"""

import time
import numpy as np
from scipy import stats

N_PLANES = 6
Y_SPACING = 200.0
SIGMA = 0.100

Y_true = np.array([Y_SPACING * i for i in range(1, N_PLANES + 1)])
W_true = np.array([1000.0 * i / N_PLANES for i in range(1, N_PLANES + 1)])

NDOF = 2 * N_PLANES - 3   # = 9
POWER = NDOF / 2.0        # = 4.5

TARGET_RATIO_NEAREST = 10.0
LAMBDA = W_true[0] / (TARGET_RATIO_NEAREST * SIGMA)
W_shrunk = W_true / LAMBDA
Y_shrunk = Y_true / LAMBDA

print(f"ndof = {NDOF}, power = {POWER}")
print(f"LAMBDA = {LAMBDA:.4f}")


def batch_fit_line(n_trials, rng, W, Y):
    n = len(Y)
    X = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    Z = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    w = 1.0 / SIGMA ** 2

    b = np.sum(X * Y[None, :], axis=1) / np.sum(Y * Y)
    resid_bend = X - b[:, None] * Y[None, :]
    chi2_bend = np.sum(resid_bend ** 2, axis=1) * w

    r = np.sqrt(X ** 2 + Y[None, :] ** 2)
    n_pts = n
    Sz = np.sum(Z, axis=1); Szz = np.sum(Z * Z, axis=1)
    Sr = np.sum(r, axis=1); Srz = np.sum(r * Z, axis=1)
    det = n_pts * Szz - Sz ** 2
    q = (n_pts * Srz - Sz * Sr) / det
    p = (Sr - q * Sz) / n_pts
    resid_depth = r - (p[:, None] + q[:, None] * Z)
    chi2_depth = np.sum(resid_depth ** 2, axis=1) * w

    return chi2_bend + chi2_depth


# --- Step 0: scaling law sanity check ---
print("\n--- Step 0: scaling law check ---")
for lam_test in [LAMBDA, LAMBDA * 7.0]:
    Wt = W_true / lam_test
    Yt = Y_true / lam_test
    c2 = batch_fit_line(500_000, np.random.default_rng(1), Wt, Yt)
    # take a LEFT-tail (small chi2) percentile this time
    s = np.sort(c2)[5000] * lam_test ** 2   # 1st percentile, left tail
    print(f"  lambda={lam_test:9.2f}  chi2_shrunk[1pct]*lambda^2 = {s:.6f}")

# --- Step 1: high-statistics power-law plateau fit (LEFT tail) ---
print("\n--- Step 1: chi2_total power-law plateau (fixed power=4.5), LEFT tail ---")
t0 = time.time()
N_TRIALS = 60_000_000
BATCH = 3_000_000
SEED = 20261001
rng = np.random.default_rng(SEED)
chi2_all = np.empty(N_TRIALS, dtype=np.float64)
done = 0
while done < N_TRIALS:
    nb = min(BATCH, N_TRIALS - done)
    chi2_all[done:done + nb] = batch_fit_line(nb, rng, W_shrunk, Y_shrunk)
    done += nb
print(f"Generated {N_TRIALS:,} trials in {time.time()-t0:.0f}s")

chi2_sorted = np.sort(chi2_all)
n = len(chi2_sorted)

event_counts = np.unique(np.geomspace(20, int(0.02 * n), 60).astype(int))
cuts = chi2_sorted[event_counts]          # LEFT tail: ne-th SMALLEST chi2
p_emp = event_counts / n
K_local = p_emp / cuts ** POWER
rel_err = 1.0 / np.sqrt(event_counts)

print(f"\n{'n_events':>10} {'cut(chi2)':>12} {'P_emp':>12} {'K_local':>14} {'rel.err':>8}")
for ne, cc, pp, kk, re in zip(event_counts, cuts, p_emp, K_local, rel_err):
    print(f"{ne:>10d} {cc:>12.6f} {pp:>12.3e} {kk:>14.4e} {re*100:>7.1f}%")

deep = cuts < np.percentile(cuts, 25)
w_ = event_counts[deep].astype(float)
K_LINE_ORIGIN = float(np.sum(K_local[deep] * w_) / np.sum(w_))
print(f"\nPlateau K_shrunk (deep quartile, cut<{np.percentile(cuts,25):.4f}) = {K_LINE_ORIGIN:.4e}")

# --- Step 2: combine to true-detector numbers, matched 3-sigma cut ---
print("\n--- Step 2: combine to true-detector numbers ---")
P_3SIGMA = stats.norm.sf(3)
CUT = stats.chi2.isf(P_3SIGMA, NDOF)
print(f"ndof={NDOF}, matched one-sided-3-sigma chi2 cut = {CUT:.4f}")


def P_true_line_origin(cut, lam=LAMBDA):
    return K_LINE_ORIGIN * (cut / lam ** 2) ** POWER


p_true = P_true_line_origin(CUT)
n1 = (1.0 / p_true) ** (1.0 / N_PLANES)
n10 = (10.0 / p_true) ** (1.0 / N_PLANES)
print(f"P_true(origin-constrained straight line, chi2<{CUT:.2f}) = {p_true:.4e}")
print(f"n(E[#fakes]=1)  = {n1:.0f} hits/layer")
print(f"n(E[#fakes]=10) = {n10:.0f} hits/layer")

print(f"\nK_LINE_ORIGIN (save this) = {K_LINE_ORIGIN:.6e}")
