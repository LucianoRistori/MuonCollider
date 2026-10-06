"""
helix_tower.py

Same 6-layer projective tower geometry as projective_tower.py, but now
fitting a first-order-in-curvature ("parabolic helix") model instead of a
straight line:
    x(z) = a + b*z + c*z^2
    y(z) = d + e*z + f*z^2
fit independently (6 free parameters total -- an unconstrained/conservative
stand-in for the true 5-parameter helix, see discussion in chat), so
ndof = 2*N - 6 (= 6 for N=6), power = ndof/2 = 3.

Adds a curvature (radius-of-curvature) acceptance cut: only accept a
combination if its effective radius of curvature R = 1/(2*sqrt(c^2+f^2))
exceeds R_MIN = 10000 mm (10 m) -- i.e. reject anything that curves too
sharply to be a plausible high-momentum real track.

Step 1: verify the chi2 tail power-law (fixed exponent=3) plateaus, same
        method as milestone 4/6/7.
Step 2: characterize how a curvature cut suppresses the chi2<cut
        probability -- check the kappa_max^2 scaling (expected from a
        "locally flat near kappa=0" argument, analogous to the chi2
        cut^(ndof/2) argument) and rough independence from the chi2 cut
        (same style check as the slope-cut factorization, milestone 5).
Step 3: combine to get numbers for the true (10m-plane) geometry.
"""

import time
import numpy as np
from tracksim import Plane, Detector, batch_random_poly_fit

N_PLANES = 6
SIGMA = 0.100  # mm, same for every layer (unchanged from projective_tower.py)
Z_SPACING = 200.0  # mm
DEGREE = 2
NDOF = 2 * N_PLANES - 2 * (DEGREE + 1)
POWER = NDOF / 2.0

# True geometry (identical to projective_tower.py)
z_true = np.array([Z_SPACING * i for i in range(1, N_PLANES + 1)])
W_true = np.array([1000.0 * i / N_PLANES for i in range(1, N_PLANES + 1)])

TARGET_RATIO_NEAREST = 10.0
LAMBDA = W_true[0] / (TARGET_RATIO_NEAREST * SIGMA)
W_shrunk = W_true / LAMBDA

print(f"N_PLANES={N_PLANES}  DEGREE={DEGREE}  NDOF={NDOF}  POWER={POWER}")
print(f"LAMBDA = {LAMBDA:.4f}")

detector_shrunk = Detector(planes=[
    Plane(z=z_true[i], size_x=W_shrunk[i], size_y=W_shrunk[i], sigma_x=SIGMA, sigma_y=SIGMA)
    for i in range(N_PLANES)
])

# ---------------------------------------------------------------------------
# Step 1: chi2 tail plateau check (fixed exponent = 3), ignoring curvature
# ---------------------------------------------------------------------------
TOTAL_TRIALS = 40_000_000
CHUNK = 2_000_000
SEED = 909090

rng = np.random.default_rng(SEED)
n_chunks = TOTAL_TRIALS // CHUNK
all_chi2 = np.empty(TOTAL_TRIALS, dtype=np.float64)
all_kappa = np.empty(TOTAL_TRIALS, dtype=np.float64)

t0 = time.time()
for i in range(n_chunks):
    chi2, cx, cy = batch_random_poly_fit(detector_shrunk, CHUNK, rng, degree=DEGREE, return_coeffs=True)
    all_chi2[i*CHUNK:(i+1)*CHUNK] = chi2
    # curvature-vector magnitude from the z^2 coefficients (index 2)
    all_kappa[i*CHUNK:(i+1)*CHUNK] = np.sqrt(cx[:, 2]**2 + cy[:, 2]**2)
print(f"Generated {TOTAL_TRIALS:,} trials in {time.time()-t0:.0f}s")

order = np.argsort(all_chi2)
sorted_chi2 = all_chi2[order]
sorted_kappa_by_chi2 = all_kappa[order]
n = len(sorted_chi2)

event_counts = np.unique(np.geomspace(20, int(0.02 * n), 60).astype(int))
cuts = sorted_chi2[event_counts]
p_emp = event_counts / n
K_local = p_emp / cuts ** POWER
rel_err = 1.0 / np.sqrt(event_counts)

print(f"\n{'n_events':>10} {'cut(chi2)':>12} {'P_emp':>12} {'K_local':>14} {'rel.err':>8}")
for ne, c, p, k, re in zip(event_counts, cuts, p_emp, K_local, rel_err):
    print(f"{ne:>10d} {c:>12.4f} {p:>12.3e} {k:>14.4e} {re*100:>7.1f}%")

deep = cuts < np.percentile(cuts, 25)
w = event_counts[deep].astype(float)
K_shrunk = np.sum(K_local[deep] * w) / np.sum(w)
print(f"\nPlateau K_shrunk (deep quartile, cut<{np.percentile(cuts,25):.2f}) = {K_shrunk:.4e}")

deep2 = (cuts >= np.percentile(cuts, 25)) & (cuts < np.percentile(cuts, 50))
w2 = event_counts[deep2].astype(float)
K_shrunk_2 = np.sum(K_local[deep2] * w2) / np.sum(w2)
print(f"Next quartile K_shrunk = {K_shrunk_2:.4e}  (ratio to plateau: {K_shrunk_2/K_shrunk:.3f})")

# ---------------------------------------------------------------------------
# Step 2: curvature-cut suppression factor f(kappa_max)
# ---------------------------------------------------------------------------
print("\n--- curvature-cut suppression check ---")
# Typical scale of kappa for these random combinations (for reference)
print(f"kappa (shrunk) percentiles: 1%={np.percentile(all_kappa,1):.3e}  "
      f"10%={np.percentile(all_kappa,10):.3e}  50%={np.percentile(all_kappa,50):.3e}")

# Test several moderate kappa_max values (loose enough for good statistics)
kappa_test_values = np.percentile(all_kappa, [0.5, 1, 2, 5, 10, 20])
chi2_cuts_test = [np.percentile(all_chi2, p) for p in [1, 5, 20, 50]]

print(f"\n{'kappa_max':>12} {'chi2_cut':>10} {'P(joint)':>12} {'P(chi2 only)':>12} {'f=ratio':>10}")
f_by_kappa = {}
for kmax in kappa_test_values:
    ratios = []
    for ccut in chi2_cuts_test:
        mask_k = all_kappa < kmax
        p_chi2_only = np.mean(all_chi2 < ccut)
        p_joint = np.mean((all_chi2 < ccut) & mask_k)
        f = p_joint / p_chi2_only
        ratios.append(f)
        print(f"{kmax:>12.3e} {ccut:>10.3f} {p_joint:>12.3e} {p_chi2_only:>12.3e} {f:>10.4f}")
    f_by_kappa[kmax] = np.mean(ratios)
    print(f"  --> mean f over chi2 cuts tested: {np.mean(ratios):.4f}  (spread {np.std(ratios):.4f})")

# Check kappa_max^2 scaling: f(kmax)/kmax^2 should be ~constant
print(f"\n{'kappa_max':>12} {'f(kmax)':>10} {'f/kmax^2':>14}")
kmax_arr = np.array(list(f_by_kappa.keys()))
f_arr = np.array(list(f_by_kappa.values()))
for kmax, f in zip(kmax_arr, f_arr):
    print(f"{kmax:>12.3e} {f:>10.4f} {f/kmax**2:>14.4e}")

np.save("helix_tower_chi2.npy", all_chi2)
np.save("helix_tower_kappa.npy", all_kappa)
print("\nSaved raw chi2/kappa arrays for follow-up analysis.")
