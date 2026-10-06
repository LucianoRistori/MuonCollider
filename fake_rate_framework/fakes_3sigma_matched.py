"""
fakes_3sigma_matched.py

Redo the projective-tower fake-rate worked examples (milestones 7, 8, 9)
using a chi2 cut matched to a genuine ONE-SIDED 3-sigma tail probability
for each model's ndof, instead of the previously-used fixed chi2=24
(which was only an approximate round-number stand-in, as Luciano's
follow-up question flagged).

One-sided 3-sigma tail probability: p = norm.sf(3) = 1.34990e-3
Matched chi2 cuts (from scipy chi2.isf(p, ndof)):
  ndof=6 (helix, unconstrained-quadratic model): cut = 21.7392
  ndof=7 (true constrained 5-param helix, for reference only):  cut = 23.5802
  ndof=8 (straight line):                        cut = 25.3611

Reuses all previously-calibrated constants (K, lambda, F_kappa) -- no new
MC needed, since the cut only enters the already-derived closed-form
P_true(cut) formulas.
"""

import numpy as np
from scipy import stats

N_PLANES = 6
LAMBDA = 166.6667

P_3SIGMA = stats.norm.sf(3)
CUT_LINE = stats.chi2.isf(P_3SIGMA, 8)     # ndof=8
CUT_HELIX = stats.chi2.isf(P_3SIGMA, 6)    # ndof=6 (unconstrained quadratic model used in milestone 9)

print(f"one-sided 3-sigma tail probability p = {P_3SIGMA:.6e}")
print(f"matched chi2 cut, ndof=8 (line)  = {CUT_LINE:.4f}")
print(f"matched chi2 cut, ndof=6 (helix) = {CUT_HELIX:.4f}")
print(f"(for reference, ndof=7 true constrained helix -> {stats.chi2.isf(P_3SIGMA,7):.4f})")

# --- straight-line calibration (milestone 7/8) ---
POWER_LINE = 4.0
K_LINE = 1.122e-11

def P_true_line(cut, lam=LAMBDA):
    return K_LINE * (cut / lam**2) ** POWER_LINE

# --- helix (parabolic, unconstrained quadratic) + R>10m calibration (milestone 9) ---
POWER_HELIX = 3.0
K_HELIX = 9.95e-9
F_KAPPA = 2.1558e-3   # R>10m suppression, ~cut-independent (checked in helix_tower.py)

def P_true_helix(cut, lam=LAMBDA):
    return F_KAPPA * K_HELIX * (cut / lam**2) ** POWER_HELIX

print("\n--- Milestone 7/8 redo: straight line, matched cut = {:.2f} ---".format(CUT_LINE))
for sigma_label, lam in [("100um", LAMBDA), ("10um", LAMBDA * 10.0)]:
    p = P_true_line(CUT_LINE, lam)
    n1 = (1.0 / p) ** (1.0 / N_PLANES)
    n10 = (10.0 / p) ** (1.0 / N_PLANES)
    print(f"  sigma={sigma_label:>6}: P_true(cut)={p:.4e}   n(E=1)={n1:.0f}   n(E=10)={n10:.0f}")

print("\n--- Milestone 9 redo: helix + R>10m, matched cut = {:.2f} ---".format(CUT_HELIX))
p_helix = P_true_helix(CUT_HELIX)
n1_helix = (1.0 / p_helix) ** (1.0 / N_PLANES)
n10_helix = (10.0 / p_helix) ** (1.0 / N_PLANES)
print(f"  sigma=100um: P_true(cut)={p_helix:.4e}   n(E=1)={n1_helix:.0f}   n(E=10)={n10_helix:.0f}")

p_line_100 = P_true_line(CUT_LINE, LAMBDA)
print(f"\nratio helix/line per-combination probability (matched-3sigma cuts): "
      f"{p_helix/p_line_100:.1f}x")
print(f"ratio of n(E=1) crossing points (line/helix): {(n1_helix and (1.0/p_line_100)**(1.0/N_PLANES)/n1_helix) or 0:.2f}x fewer hits/layer tolerated for helix")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Left panel: line, sigma=100um vs 10um, matched cut
ax = axes[0]
n_range = np.geomspace(1000, 300000, 400)
for sigma_label, lam, color in [("sigma = 100 μm", LAMBDA, "#4C72B0"),
                                  ("sigma = 10 μm", LAMBDA * 10.0, "#C44E52")]:
    p = P_true_line(CUT_LINE, lam)
    efakes = (n_range ** N_PLANES) * p
    n1 = (1.0 / p) ** (1.0 / N_PLANES)
    ax.loglog(n_range, efakes, color=color, lw=2.5, label=f"{sigma_label} (n@E=1: {n1:.0f})")
    ax.axvline(n1, color=color, ls="--", lw=1.1, alpha=0.7)
ax.axhline(1, color="gray", ls=":", lw=1, alpha=0.6)
ax.axhline(10, color="gray", ls=":", lw=1, alpha=0.6)
ax.set_xlabel("hits per layer")
ax.set_ylabel("E[# fake tracks]")
ax.set_title(f"Straight line, matched 3σ cut (chi2={CUT_LINE:.1f}, ndof=8)")
ax.legend()
ax.grid(True, which="both", alpha=0.3)

# Right panel: line vs helix+Rcut, each at its own matched cut, sigma=100um
ax2 = axes[1]
n_range2 = np.geomspace(300, 100000, 400)
p_line = P_true_line(CUT_LINE)
p_hel = P_true_helix(CUT_HELIX)
for label, p, color in [(f"straight line (ndof=8, cut={CUT_LINE:.1f})", p_line, "#4C72B0"),
                          (f"helix+R>10m (ndof=6, cut={CUT_HELIX:.1f})", p_hel, "#C44E52")]:
    efakes = (n_range2 ** N_PLANES) * p
    n1 = (1.0 / p) ** (1.0 / N_PLANES)
    ax2.loglog(n_range2, efakes, color=color, lw=2.5, label=f"{label} (n@E=1: {n1:.0f})")
    ax2.axvline(n1, color=color, ls="--", lw=1.1, alpha=0.7)
ax2.axhline(1, color="gray", ls=":", lw=1, alpha=0.6)
ax2.axhline(10, color="gray", ls=":", lw=1, alpha=0.6)
ax2.set_xlabel("hits per layer")
ax2.set_ylabel("E[# fake tracks]")
ax2.set_title("Line vs. helix+R>10m, each at its OWN matched 3σ cut\n(sigma=100μm)")
ax2.legend()
ax2.grid(True, which="both", alpha=0.3)

fig.suptitle("Fake-track estimates re-derived with a true one-sided 3σ-equivalent chi2 cut per ndof\n"
             "(replaces the earlier fixed chi2=24 approximation)", y=1.02)
fig.tight_layout()
fig.savefig("fakes_3sigma_matched.png", dpi=150, bbox_inches="tight")
print("\nSaved plot to fakes_3sigma_matched.png")
