# CORRECTION (10 October 2026): the depth-view fit in this framework had r and z
# swapped (it fitted r = p + q z, with residuals in r instead of z), so the fake
# probabilities computed or used here are about 450 times too low for the helix.
# Kept as the record of the earlier paper draft; the corrected values come from
# fake_rate_framework/geometric_K.py and are used by fake_rate.py.
"""
barrel_helix_tower.py

Worked numerical example for paper_draft.md Section 10/11: the same
6-plane projective-tower barrel geometry used throughout Section 7,
now with a uniform B=5 T axial field and the physics cut "muon,
p_T > 10 GeV/c" (the large-transverse-momentum / small-curvature regime
of Section 10), for tracks constrained to originate on the beam axis
(X0=Y0=0, only Z0 free).

Track model (Section 10.4): the fit splits into two independent linear
least-squares problems --
  bending view:  X(Y) = b*Y + c*Y^2         (no intercept; 2 free params)
  depth view:    r(Z) = p + q*Z,  r=sqrt(X^2+Y^2)   (2 free params)
with ndof = (N-2) + (N-2) = 2N-4, identical in form to the original
straight-line model (Section 2.2).

Geometry: N=6 planes at Y_i=200,400,...,1200 mm (same spacing as
Section 7.4/9.1's projective tower), sizes W_i=1000*i/6 mm (X- and
Z-extent, projective cone from the origin), sigma=100 micron uniform.
Noise hits: (X_i,Z_i) uniform over plane i's W_i x W_i active area.

Curvature acceptance cut: kappa_proxy = 2*c (dropping the (1+b^2)^1.5
slope-correction factor of Section 10.3 for simplicity, same
simplification used in the earlier forward-case curvature cut,
milestone 9/helix_tower.py); accept if |kappa_proxy| < kappa_max,
kappa_max set by B=5 T, p_T>10 GeV/c: R_min = p_T/(0.2998*B) [m],
kappa_max = 1/R_min.

Depth-view weighting: r_i = sqrt(X_i^2+Y_i^2) is a derived quantity;
Section 10.4 notes the exact propagated resolution sigma_r depends on
the (random) hit itself. Here we use the simplest fixed-weight choice,
weighting the depth-view residual by 1/sigma_X^2 (same as the bending
view) -- an approximation flagged explicitly in the writeup, kept
because it preserves the exact data-linearity the scaling-law argument
needs.
"""

import time
import numpy as np

# ---------------------------------------------------------------------------
# Geometry (same projective tower as Section 7.4/7.5/9.1)
# ---------------------------------------------------------------------------
N_PLANES = 6
Y_SPACING = 200.0  # mm
SIGMA = 0.100       # mm (100 micron), same for X and Z, every plane

Y_true = np.array([Y_SPACING * i for i in range(1, N_PLANES + 1)])          # 200..1200
W_true = np.array([1000.0 * i / N_PLANES for i in range(1, N_PLANES + 1)])  # 166.7..1000

NDOF = 2 * N_PLANES - 4   # = 8, same formula as the straight-line model
POWER = NDOF / 2.0        # = 4

# ---------------------------------------------------------------------------
# Physics: B=5T, muon p_T > 10 GeV/c  ->  kappa_max
# ---------------------------------------------------------------------------
B_TESLA = 5.0
PT_MIN_GEV = 10.0
R_MIN_M = PT_MIN_GEV / (0.2998 * B_TESLA)     # standard R[m]=pT[GeV]/(0.2998*B[T])
R_MIN_MM = R_MIN_M * 1000.0
KAPPA_MAX_TRUE = 1.0 / R_MIN_MM               # 1/mm, at the TRUE (unshrunk) Y-scale

print(f"B = {B_TESLA} T, p_T,min = {PT_MIN_GEV} GeV/c")
print(f"R_min = {R_MIN_M:.4f} m = {R_MIN_MM:.1f} mm")
print(f"kappa_max (true scale) = {KAPPA_MAX_TRUE:.6e} /mm")

# Sagitta at the last plane for a radial (b=0) track at kappa_max, for scale:
sagitta_last = 0.5 * KAPPA_MAX_TRUE * Y_true[-1] ** 2
print(f"Reference sagitta at Y={Y_true[-1]:.0f} mm, kappa_max, b=0: {sagitta_last:.2f} mm "
      f"(vs plane half-width {W_true[-1]/2:.1f} mm)")

# ---------------------------------------------------------------------------
# Shrunk geometry -- IMPORTANT, and different from every earlier milestone:
# r=sqrt(X^2+Y^2) mixes the "noise width" coordinate X with the fixed plane
# position Y in a way that is NOT homogeneous if only W is shrunk while Y is
# held at its true value (verified below: doing that gives a spurious extra
# factor in the scaling check). The reason is that r carries an intrinsic
# ABSOLUTE length scale (distance to the beam axis) that the bending-view
# model X=bY+cY^2 never had (that relation is satisfied at any absolute
# Y-scale by suitable b,c). The fix: shrink X, Y, AND Z together by the same
# common factor lambda (the whole geometry, not just the transverse width),
# holding sigma fixed at its true physical value. Under this transformation
# r_shrunk = r_true/lambda exactly, restoring the needed homogeneity for
# BOTH views simultaneously (verified in Step 0 below).
# ---------------------------------------------------------------------------
TARGET_RATIO_NEAREST = 10.0
LAMBDA = W_true[0] / (TARGET_RATIO_NEAREST * SIGMA)
W_shrunk = W_true / LAMBDA
Y_shrunk_geom = Y_true / LAMBDA
print(f"\nLAMBDA = {LAMBDA:.4f}")
print(f"W_shrunk = {np.round(W_shrunk, 4)} mm, Y_shrunk = {np.round(Y_shrunk_geom, 4)} mm "
      f"(sigma={SIGMA} mm fixed -- whole geometry shrunk together)")


def batch_fit(n_trials, rng, W, Y):
    """
    Vectorized: one random hit (X_i,Z_i) per plane, fit bending view
    (no-intercept quadratic X=b*Y+c*Y^2) and depth view (r=p+q*Z,
    r=sqrt(X^2+Y^2)), both weighted by 1/SIGMA^2. Returns chi2_total,
    chi2_bend, chi2_depth, b, c, p, q (each shape (n_trials,) except chi2
    which is summed already where noted).
    """
    n = len(Y)
    X = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    Z = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    w = 1.0 / SIGMA ** 2  # scalar, same for every plane/coord

    # --- bending view: X = b*Y + c*Y^2, no intercept ---
    Y1 = Y
    Y2 = Y ** 2
    S11 = np.sum(Y1 * Y1) * w
    S12 = np.sum(Y1 * Y2) * w
    S22 = np.sum(Y2 * Y2) * w
    RX1 = np.sum(X * Y1[None, :], axis=1) * w
    RX2 = np.sum(X * Y2[None, :], axis=1) * w
    det = S11 * S22 - S12 * S12
    b = (RX1 * S22 - RX2 * S12) / det
    c = (RX2 * S11 - RX1 * S12) / det
    model_bend = b[:, None] * Y1[None, :] + c[:, None] * Y2[None, :]
    chi2_bend = w * np.sum((X - model_bend) ** 2, axis=1)

    # --- depth view: r = p + q*Z, ordinary intercept+slope fit ---
    r = np.sqrt(X ** 2 + Y[None, :] ** 2)
    Sw = n * w
    Swz = w * np.sum(Z, axis=1)
    Swzz = w * np.sum(Z * Z, axis=1)
    Swr = w * np.sum(r, axis=1)
    Swzr = w * np.sum(Z * r, axis=1)
    denom = Sw * Swzz - Swz ** 2
    q = (Sw * Swzr - Swz * Swr) / denom
    p = (Swr - q * Swz) / Sw
    model_depth = p[:, None] + q[:, None] * Z
    chi2_depth = w * np.sum((r - model_depth) ** 2, axis=1)

    chi2_total = chi2_bend + chi2_depth
    return chi2_total, chi2_bend, chi2_depth, b, c, p, q


# ---------------------------------------------------------------------------
# Step 0: sanity check -- exact scaling law for this new (two-view) model
# ---------------------------------------------------------------------------
print("\n--- Step 0: exact scaling law check (two shrink factors) ---")
print("(chi2_shrunk * lambda^2 should be lambda-independent -- NOT chi2_shrunk/lambda^2,")
print(" since chi2_shrunk is already chi2_true/lambda^2 by construction)")
rng0 = np.random.default_rng(42)
chi2_a, *_ = batch_fit(2_000_000, rng0, W_shrunk, Y_shrunk_geom)
rng0b = np.random.default_rng(43)
LAMBDA2 = LAMBDA * 7.0
chi2_b, *_ = batch_fit(2_000_000, rng0b, W_true / LAMBDA2, Y_true / LAMBDA2)
s_a = chi2_a * LAMBDA ** 2
s_b = chi2_b * LAMBDA2 ** 2
for pct in [1, 10, 50, 90, 99]:
    va, vb = np.percentile(s_a, pct), np.percentile(s_b, pct)
    print(f"  pct={pct:>3d}  s(lambda={LAMBDA:.1f})={va:.5f}   "
          f"s(lambda={LAMBDA2:.1f})={vb:.5f}   ratio={va/vb:.4f}")

# ---------------------------------------------------------------------------
# Step 0b: validity of the first-order-in-kappa approximation at kappa_max,
# for a real track (not noise) at the farthest plane. Exact circle in the
# transverse (X,Y) plane vs the first-order X(Y)=bY+cY^2 formula.
# ---------------------------------------------------------------------------
print("\n--- Step 0b: first-order approximation validity at kappa_max, Y_max ---")
Y_max = Y_true[-1]
for b_test in [0.0, 0.3, 0.6]:
    kap = KAPPA_MAX_TRUE
    # exact: solve s(Y) by Newton from Y(s)=Y0+ (1/kap)[sin(phi0+kap s)-sin(phi0)]
    phi0 = np.arctan(b_test)  # dX/dY=b, dY/ds=cos(phi0) convention (Sec.10 notes)
    from scipy.optimize import brentq
    def Yfun(s):
        return (np.sin(phi0 + kap * s) - np.sin(phi0)) / kap
    s_guess = Y_max  # first-order guess
    s_exact = brentq(lambda s: Yfun(s) - Y_max, 0, 3 * Y_max)
    X_exact = (np.cos(phi0) - np.cos(phi0 + kap * s_exact)) / kap
    c_pred = 0.5 * kap * (1 + b_test ** 2) ** 1.5
    X_approx = b_test * Y_max + c_pred * Y_max ** 2
    print(f"  b={b_test:.2f}: X_exact={X_exact:8.3f} mm  X_1st-order={X_approx:8.3f} mm  "
          f"diff={X_exact-X_approx:+7.3f} mm ({SIGMA*1000:.0f} um resolution)")

# ---------------------------------------------------------------------------
# Step 1: chi2_total tail plateau (fixed exponent POWER=4), on shrunk geometry
# ---------------------------------------------------------------------------
print("\n--- Step 1: chi2_total power-law plateau (fixed power=4) ---")
TOTAL_TRIALS = 60_000_000
CHUNK = 3_000_000
SEED = 20261001

rng = np.random.default_rng(SEED)
n_chunks = TOTAL_TRIALS // CHUNK
all_chi2 = np.empty(TOTAL_TRIALS, dtype=np.float64)
all_kappa_proxy = np.empty(TOTAL_TRIALS, dtype=np.float64)

t0 = time.time()
for i in range(n_chunks):
    chi2_t, chi2_bend, chi2_depth, b, c, p, q = batch_fit(CHUNK, rng, W_shrunk, Y_shrunk_geom)
    all_chi2[i * CHUNK:(i + 1) * CHUNK] = chi2_t
    all_kappa_proxy[i * CHUNK:(i + 1) * CHUNK] = 2.0 * c
print(f"Generated {TOTAL_TRIALS:,} trials in {time.time()-t0:.0f}s")

order = np.argsort(all_chi2)
sorted_chi2 = all_chi2[order]
sorted_kappa = all_kappa_proxy[order]
n = len(sorted_chi2)

event_counts = np.unique(np.geomspace(20, int(0.02 * n), 60).astype(int))
cuts = sorted_chi2[event_counts]
p_emp = event_counts / n
K_local = p_emp / cuts ** POWER
rel_err = 1.0 / np.sqrt(event_counts)

print(f"\n{'n_events':>10} {'cut(chi2)':>12} {'P_emp':>12} {'K_local':>14} {'rel.err':>8}")
for ne, cc, pp, kk, re in zip(event_counts, cuts, p_emp, K_local, rel_err):
    print(f"{ne:>10d} {cc:>12.4f} {pp:>12.3e} {kk:>14.4e} {re*100:>7.1f}%")

deep = cuts < np.percentile(cuts, 25)
w_ = event_counts[deep].astype(float)
K_shrunk = np.sum(K_local[deep] * w_) / np.sum(w_)
print(f"\nPlateau K_shrunk (deep quartile, cut<{np.percentile(cuts,25):.2f}) = {K_shrunk:.4e}")

# ---------------------------------------------------------------------------
# Step 2: kappa_max cut suppression factor f(kappa_max) -- expect LINEAR
# scaling (1D interval cut on c, unlike the 2D-disk R_min cut of the old
# forward-case model), since only ONE view carries curvature here.
# ---------------------------------------------------------------------------
print("\n--- Step 2: kappa_max cut suppression f(kappa_max) ---")
print(f"kappa_proxy (shrunk) percentiles: "
      f"1%={np.percentile(np.abs(all_kappa_proxy),1):.3e}  "
      f"10%={np.percentile(np.abs(all_kappa_proxy),10):.3e}  "
      f"50%={np.percentile(np.abs(all_kappa_proxy),50):.3e}")

kappa_test_values = np.percentile(np.abs(all_kappa_proxy), [0.5, 1, 2, 5, 10, 20])
chi2_cuts_test = [np.percentile(all_chi2, p) for p in [1, 5, 20, 50]]

print(f"\n{'kappa_max':>12} {'chi2_cut':>10} {'P(joint)':>12} {'P(chi2 only)':>12} {'f=ratio':>10}")
f_by_kappa = {}
for kmax in kappa_test_values:
    ratios = []
    for ccut in chi2_cuts_test:
        mask_k = np.abs(all_kappa_proxy) < kmax
        p_chi2_only = np.mean(all_chi2 < ccut)
        p_joint = np.mean((all_chi2 < ccut) & mask_k)
        f = p_joint / p_chi2_only
        ratios.append(f)
        print(f"{kmax:>12.3e} {ccut:>10.3f} {p_joint:>12.3e} {p_chi2_only:>12.3e} {f:>10.4f}")
    f_by_kappa[kmax] = np.mean(ratios)
    print(f"  --> mean f over chi2 cuts tested: {np.mean(ratios):.4f}  (spread {np.std(ratios):.4f})")

print(f"\n{'kappa_max':>12} {'f(kmax)':>10} {'f/kmax':>14} {'f/kmax^2':>14}")
kmax_arr = np.array(list(f_by_kappa.keys()))
f_arr = np.array(list(f_by_kappa.values()))
for kmax, f in zip(kmax_arr, f_arr):
    print(f"{kmax:>12.3e} {f:>10.4f} {f/kmax:>14.4e} {f/kmax**2:>14.4e}")

# Use the smallest-kmax point as the calibration point (deepest into the
# small-kmax regime, consistent with house style)
C_KAPPA_LINEAR = f_arr[0] / kmax_arr[0]
print(f"\nC_kappa (linear-law calibration, f(kmax)=C_kappa*kmax, smallest kmax) = {C_KAPPA_LINEAR:.4e}")

# ---------------------------------------------------------------------------
# Step 3: combine, translate to the TRUE detector + physical kappa_max
# ---------------------------------------------------------------------------
print("\n--- Step 3: combine to true-detector numbers ---")
KAPPA_SHRUNK_MAX = KAPPA_MAX_TRUE * LAMBDA  # kappa has units 1/length -> scales AS lambda
                                             # (opposite of the old forward-case model,
                                             # because here the WHOLE geometry -- incl. Y --
                                             # is shrunk, not just the transverse width)
F_KAPPA = C_KAPPA_LINEAR * KAPPA_SHRUNK_MAX
print(f"kappa_max (shrunk scale) = {KAPPA_SHRUNK_MAX:.4e} /mm")
print(f"F_kappa (p_T>10 GeV/c, B=5T acceptance fraction) = {F_KAPPA:.4e}")

from scipy import stats
CUT_3SIGMA = stats.chi2.isf(stats.norm.sf(3), NDOF)
print(f"ndof={NDOF}, matched one-sided-3-sigma chi2 cut = {CUT_3SIGMA:.4f}")


def P_true_helix_barrel(cut):
    return F_KAPPA * K_shrunk * (cut / LAMBDA ** 2) ** POWER


p_true = P_true_helix_barrel(CUT_3SIGMA)
print(f"P_true(chi2<{CUT_3SIGMA:.2f} AND p_T>10GeV) = {p_true:.4e}")

n1 = (1.0 / p_true) ** (1.0 / N_PLANES)
n10 = (10.0 / p_true) ** (1.0 / N_PLANES)
print(f"n(E[#fakes]=1)  = {n1:.0f} hits/layer")
print(f"n(E[#fakes]=10) = {n10:.0f} hits/layer")

# --- compare to the existing straight-line (B=0) numbers, same geometry/ndof ---
K_LINE = 1.122e-11  # from Section 7.4 / milestone 6-7, same geometry, same LAMBDA
p_true_line = K_LINE * (CUT_3SIGMA / LAMBDA ** 2) ** POWER
n1_line = (1.0 / p_true_line) ** (1.0 / N_PLANES)
n10_line = (10.0 / p_true_line) ** (1.0 / N_PLANES)
print(f"\n[comparison] straight-line (B=0), same geometry/ndof/cut:")
print(f"P_true(line) = {p_true_line:.4e}   n(E=1)={n1_line:.0f}   n(E=10)={n10_line:.0f}")
print(f"ratio P_true(helix,5T,pT>10)/P_true(line) = {p_true/p_true_line:.3f}")
print(f"ratio n(E=1): line/helix = {n1_line/n1:.3f}")

np.save("barrel_helix_chi2.npy", all_chi2)
np.save("barrel_helix_kappa.npy", all_kappa_proxy)
print("\nSaved raw chi2/kappa arrays.")
