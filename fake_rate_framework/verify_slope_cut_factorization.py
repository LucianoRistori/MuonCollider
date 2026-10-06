"""
verify_slope_cut_factorization.py

Check the claim that a slope/angle acceptance cut factorizes exactly out
of the chi2<cut probability:

  P(chi2<cut AND |slope_x|<smax AND |slope_y|<smax) = f(smax) * P(chi2<cut)

with f(smax) = (Delta/W)*(2-Delta/W) per axis, Delta = smax*(z_last-z_first),
independent of the chi2 cut value.

Since the slope acceptance only "bites" when Delta is comparable to or
smaller than the plane size W, we test at the 1mm plane scale with a
Delta chosen to be a sizeable, easily-measurable fraction of 1mm (not the
physical 5-degree example, which only matters for the full 1000mm plane
-- this is purely a test of the mathematical factorization claim).
"""

import numpy as np
from tracksim import Detector, batch_random_chi2

N_PLANES = 8
W = 1.0          # mm
SIGMA = 0.100    # mm
Z_SPACING = 100.0
N_TRIALS = 4_000_000
SEED = 55555

detector = Detector.make_uniform(n_planes=N_PLANES, size_x=W, size_y=W,
                                  sigma_x=SIGMA, sigma_y=SIGMA, z_spacing=Z_SPACING)

rng = np.random.default_rng(SEED)
z = np.array([p.z for p in detector.planes])
z_first, z_last = z[0], z[-1]
baseline = z_last - z_first

# Generate trials and get chi2 AND the fitted slopes (need a small addition:
# recompute slopes directly here using the same weighted-fit formulas, since
# batch_random_chi2 only returns chi2).
n = N_PLANES
sigma_x = np.array([p.sigma_x for p in detector.planes])
sigma_y = np.array([p.sigma_y for p in detector.planes])
size_x = np.array([p.size_x for p in detector.planes])
size_y = np.array([p.size_y for p in detector.planes])

x = rng.uniform(-0.5, 0.5, size=(N_TRIALS, n)) * size_x[None, :]
y = rng.uniform(-0.5, 0.5, size=(N_TRIALS, n)) * size_y[None, :]
wx = 1.0 / sigma_x ** 2
wy = 1.0 / sigma_y ** 2

def batch_fit(coord, w):
    Sw = np.sum(w)
    Swz = np.sum(w * z)
    Swzz = np.sum(w * z * z)
    Swc = np.sum(w[None, :] * coord, axis=1)
    Swzc = np.sum(w[None, :] * z[None, :] * coord, axis=1)
    denom = Sw * Swzz - Swz ** 2
    b = (Sw * Swzc - Swz * Swc) / denom
    a = (Swc - b * Swz) / Sw
    resid = coord - (a[:, None] + b[:, None] * z[None, :])
    chi2 = np.sum(w[None, :] * resid ** 2, axis=1)
    return a, b, chi2

a_x, b_x, chi2_x = batch_fit(x, wx)
a_y, b_y, chi2_y = batch_fit(y, wy)
chi2 = chi2_x + chi2_y

# Test several Delta values (comparable to W=1mm, so the cut actually bites)
for Delta in [0.05, 0.1, 0.2, 0.3, 0.5]:
    smax = Delta / baseline
    f_theory_1d = (Delta / W) * (2 - Delta / W)
    f_theory_2d = f_theory_1d ** 2

    slope_ok = (np.abs(b_x) < smax) & (np.abs(b_y) < smax)

    for cut in [15, 20, 30, 50]:
        pass_chi2 = chi2 < cut
        p_unrestricted = pass_chi2.mean()
        p_restricted = (pass_chi2 & slope_ok).mean()
        ratio_empirical = p_restricted / p_unrestricted if p_unrestricted > 0 else float('nan')
        print(f"Delta={Delta:>5.2f}mm (smax={smax:.4f} rad)  cut={cut:>3d}  "
              f"empirical f={ratio_empirical:.4f}  theory f={f_theory_2d:.4f}  "
              f"(n_pass_unrestricted={pass_chi2.sum()})")
    print()
