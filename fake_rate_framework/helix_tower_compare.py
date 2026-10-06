"""
helix_tower_compare.py

Final comparison for the 6-layer projective tower: straight-line fit
(ndof=8) vs first-order-in-curvature ("parabolic helix") fit with an
R > 10 m acceptance cut (ndof=6), both at chi2 cut=24, sigma=100um.

Uses the calibration from helix_tower.py:
  chi2 tail:      P_true(chi2<cut) = K_shrunk * (cut/LAMBDA^2)^3
  curvature cut:  f(kappa_max) = C_kappa * kappa_shrunk_max^2   (small-kmax
                   quadratic law, verified numerically in helix_tower.py --
                   f/kmax^2 flat to ~1% at the smallest kmax tested)
  joint:          P_true(chi2<cut AND R>10m) = f_kappa * P_true(chi2<cut)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_PLANES = 6
LAMBDA = 166.6667

# --- straight-line calibration (milestone 7/8) ---
NDOF_LINE = 8
POWER_LINE = 4.0
K_LINE = 1.122e-11

def P_true_line(cut):
    return K_LINE * (cut / LAMBDA**2) ** POWER_LINE

# --- helix (parabolic, unconstrained quadratic) calibration (this milestone) ---
NDOF_HELIX = 6
POWER_HELIX = 3.0
K_HELIX = 9.95e-9           # deep-quartile plateau, chi2-only (curvature-unrestricted)
C_KAPPA = 2.395e10          # f(kmax)/kmax^2, from smallest-kmax numerical check
R_MIN_TRUE = 10000.0        # mm (10 m)
KAPPA_TRUE_MAX = 1.0 / (2 * R_MIN_TRUE)
KAPPA_SHRUNK_MAX = KAPPA_TRUE_MAX / LAMBDA
F_KAPPA = C_KAPPA * KAPPA_SHRUNK_MAX**2

def P_true_helix(cut):
    return F_KAPPA * K_HELIX * (cut / LAMBDA**2) ** POWER_HELIX

CUT = 24.0
p_line = P_true_line(CUT)
p_helix = P_true_helix(CUT)
print(f"chi2 cut = {CUT}")
print(f"P_true(line)              = {p_line:.4e}")
print(f"P_true(helix, R>10m)      = {p_helix:.4e}")
print(f"F_kappa (R>10m suppression) = {F_KAPPA:.4e}")
print(f"ratio helix/line          = {p_helix/p_line:.1f}x more per-combination fake probability")

n_range = np.geomspace(300, 100000, 400)
configs = [
    ("straight line (ndof=8)", p_line, "#4C72B0"),
    ("parabolic helix, R>10m (ndof=6)", p_helix, "#C44E52"),
]

fig, ax = plt.subplots(figsize=(8.5, 6.5))
for label, p, color in configs:
    efakes = (n_range ** N_PLANES) * p
    n1 = (1.0 / p) ** (1.0 / N_PLANES)
    n10 = (10.0 / p) ** (1.0 / N_PLANES)
    ax.loglog(n_range, efakes, color=color, lw=2.5, label=f"{label}  (n@E=1: {n1:.0f})")
    ax.axvline(n1, color=color, ls="--", lw=1.1, alpha=0.7)
    ax.axvline(n10, color=color, ls=":", lw=1.1, alpha=0.7)
    print(f"{label}: n(E=1)={n1:.0f}  n(E=10)={n10:.0f}")

ax.axhline(1, color="gray", ls=":", lw=1, alpha=0.6)
ax.axhline(10, color="gray", ls=":", lw=1, alpha=0.6)
ax.set_xlabel("hits per layer (same for all 6 layers)")
ax.set_ylabel("E[# fake tracks]")
ax.set_title("Expected fake tracks vs hit count per layer\n"
              f"6-layer projective tower, sigma=100μm, chi2 cut={CUT:.0f}\n"
              "straight line vs. parabolic-helix fit with R>10m acceptance")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("helix_tower_compare.png", dpi=150)
print("\nSaved plot to helix_tower_compare.png")
