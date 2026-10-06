"""
fakes_vs_hitcount.py

For the 6-layer projective tower (vertex at origin, farthest layer 1m x
1m, sigma=100um all layers, 200mm z-spacing), plot the expected number
of fake tracks vs the number of hits per plane (same n for every plane),
at a fixed chi2 cut of 24 (Luciano's "3-sigma" cut).

Uses the calibration from projective_tower.py / check_deep_tail_scatter.py:
  P_true(chi2<cut) = K * (cut/lambda^2)^(ndof/2)
with K = 1.122e-11 (replica-averaged plateau, 25 independent replicas)
and lambda = 166.6667 (shrink factor used in that calibration).

E[#fakes](n) = n^N_PLANES * P_true(cut)   [N_PLANES=6 planes, same n each]
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_PLANES = 6
NDOF = 2 * N_PLANES - 4       # = 8
POWER = NDOF / 2.0            # = 4
LAMBDA = 166.6667
K = 1.122e-11                 # replica-averaged plateau (more robust than single-run value)
CUT = 24.0                    # "3-sigma" chi2 cut

def P_true(cut):
    return K * (cut / LAMBDA**2) ** POWER

p_cut = P_true(CUT)
print(f"P_true(chi2<{CUT}) = {p_cut:.4e}")

n_values = np.geomspace(1000, 30000, 300)
efakes = (n_values ** N_PLANES) * p_cut

# find where E[#fakes] crosses 1 and 10
n_cross_1 = (1.0 / p_cut) ** (1.0 / N_PLANES)
n_cross_10 = (10.0 / p_cut) ** (1.0 / N_PLANES)
print(f"n (hits/layer) where E[#fakes]=1:  {n_cross_1:.1f}")
print(f"n (hits/layer) where E[#fakes]=10: {n_cross_10:.1f}")

print(f"\n{'n hits/layer':>14} {'E[#fakes]':>14}")
for n in [1000, 3000, 10000, int(round(n_cross_10)), 20000, 30000]:
    print(f"{n:>14d} {(n**N_PLANES)*p_cut:>14.3e}")

fig, ax = plt.subplots(figsize=(8, 6))
ax.loglog(n_values, efakes, color="#4C72B0", lw=2.5)
ax.axhline(1, color="gray", ls=":", lw=1, alpha=0.7)
ax.axhline(10, color="gray", ls=":", lw=1, alpha=0.7)
ax.axvline(n_cross_1, color="#DD8452", ls="--", lw=1.3,
           label=f"E[#fakes]=1 at n≈{n_cross_1:.0f} hits/layer")
ax.axvline(n_cross_10, color="#C44E52", ls="--", lw=1.3,
           label=f"E[#fakes]=10 at n≈{n_cross_10:.0f} hits/layer")
ax.set_xlabel("hits per layer (same for all 6 layers)")
ax.set_ylabel("E[# fake tracks]")
ax.set_title(f"Expected fake tracks vs hit count per layer\n"
             f"6-layer projective tower, chi2 cut={CUT:.0f} (~3 sigma)")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("fakes_vs_hitcount.png", dpi=150)
print("\nSaved plot to fakes_vs_hitcount.png")
