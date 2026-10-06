"""
fakes_vs_hitcount_sigma_compare.py

Extends fakes_vs_hitcount.py: compares the 6-layer projective tower
expected-fakes-vs-hitcount curve (chi2 cut=24, "3-sigma") for two
resolutions, sigma=100um (the original calibration) and sigma=10um
(10x finer), using the universal scaling law:

    chi2 = (W/sigma)^2 * s   =>   P(chi2<cut) is obtained, for fixed cut,
    by using an effective lambda that is 10x LARGER when sigma is 10x
    SMALLER (same true W, same shrink target of ~10 sigma at the nearest
    layer -- shrinking to a fixed multiple of sigma means the shrunk
    geometry is identical either way, so P_shrunk(.) is literally the
    SAME calibrated function; only lambda = W_true/W_shrunk changes,
    scaling as 1/sigma).

Equivalently: P_true(cut) = K * (cut/lambda^2)^POWER with
    lambda_10um = lambda_100um * 10
so at fixed cut, P_true scales as lambda^-2*POWER = lambda^-ndof,
i.e. P ~ sigma^ndof (sigma^8 here), exactly as derived analytically.

No new MC needed -- this reuses the existing K (replica-averaged
plateau) and LAMBDA calibration from projective_tower.py /
check_deep_tail_scatter.py, just evaluated at two different lambda
values.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_PLANES = 6
NDOF = 2 * N_PLANES - 4       # = 8
POWER = NDOF / 2.0            # = 4
K = 1.122e-11                 # replica-averaged plateau
CUT = 24.0                    # "3-sigma" chi2 cut

LAMBDA_100um = 166.6667
LAMBDA_10um = LAMBDA_100um * 10.0   # sigma 10x smaller -> lambda 10x larger

configs = [
    ("sigma = 100 μm", LAMBDA_100um, "#4C72B0"),
    ("sigma = 10 μm", LAMBDA_10um, "#C44E52"),
]

def P_true(cut, lam):
    return K * (cut / lam**2) ** POWER

print(f"{'config':>14} {'P_true(24)':>14} {'n(E=1)':>12} {'n(E=10)':>12}")
for label, lam, _ in configs:
    p_cut = P_true(CUT, lam)
    n1 = (1.0 / p_cut) ** (1.0 / N_PLANES)
    n10 = (10.0 / p_cut) ** (1.0 / N_PLANES)
    print(f"{label:>14} {p_cut:>14.4e} {n1:>12.0f} {n10:>12.0f}")

print(f"\nRatio of crossing points (10um/100um): "
      f"predicted 10^(ndof/N_PLANES) = 10^{NDOF/N_PLANES:.3f} = {10**(NDOF/N_PLANES):.2f}")

fig, ax = plt.subplots(figsize=(8.5, 6.5))

n_range = np.geomspace(1000, 300000, 400)
for label, lam, color in configs:
    p_cut = P_true(CUT, lam)
    efakes = (n_range ** N_PLANES) * p_cut
    n1 = (1.0 / p_cut) ** (1.0 / N_PLANES)
    n10 = (10.0 / p_cut) ** (1.0 / N_PLANES)
    ax.loglog(n_range, efakes, color=color, lw=2.5, label=label)
    ax.axvline(n1, color=color, ls="--", lw=1.1, alpha=0.7)
    ax.axvline(n10, color=color, ls=":", lw=1.1, alpha=0.7)

ax.axhline(1, color="gray", ls=":", lw=1, alpha=0.6)
ax.axhline(10, color="gray", ls=":", lw=1, alpha=0.6)
ax.set_xlabel("hits per layer (same for all 6 layers)")
ax.set_ylabel("E[# fake tracks]")
ax.set_title("Expected fake tracks vs hit count per layer\n"
              f"6-layer projective tower, chi2 cut={CUT:.0f} (~3 sigma)\n"
              "dashed = E[#fakes]=1 crossing, dotted = E[#fakes]=10 crossing")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("fakes_vs_hitcount_sigma_compare.png", dpi=150)
print("\nSaved plot to fakes_vs_hitcount_sigma_compare.png")
