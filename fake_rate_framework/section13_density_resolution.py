"""
section13_density_resolution.py

Section 13 worked example: E[#fakes] vs hit DENSITY (hits/mm^2, uniform
across the tower -- each plane's hit count = density * its own area)
and vs detector RESOLUTION, for the same B=5T/p_T>10GeV/c projective
tower calibrated in Sections 11-12. Both axes are fully analytic from
the existing calibration -- no new Monte Carlo.

Key facts used:
- E[#fakes] = P_true(cut) * prod_i(n_i), n_i = density * W_i^2  (Section 4).
- P_true(sigma) = P_true(sigma_ref) * (sigma/sigma_ref)^ndof  (Section 5.1),
  with F_kappa UNCHANGED by sigma: for fixed-weight (uniform sigma) least
  squares, sigma enters chi2 only as an overall positive multiplicative
  constant (1/sigma^2), which cannot move the location of the chi2
  minimum -- so the fitted (b,kappa,p,q) for a given noise realization,
  and hence the marginal kappa_fit distribution underlying F_kappa, does
  not depend on sigma at all. Only the chi2 VALUE at the minimum scales
  with sigma, which is exactly the sigma^ndof law above.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

N_PLANES = 6
Y_SPACING = 200.0
SIGMA_REF = 0.100  # mm

W_true = np.array([1000.0 * i / N_PLANES for i in range(1, N_PLANES + 1)])
Y_true = np.array([Y_SPACING * i for i in range(1, N_PLANES + 1)])
AREAS = W_true ** 2
PROD_AREAS = np.prod(AREAS)

LAMBDA_REF = W_true[0] / (10.0 * SIGMA_REF)

P_3SIGMA = stats.norm.sf(3)

# --- calibrated constants from Sections 11-12 (B=5T, p_T>10 GeV/c) ---
models = {
    "exact 4-param helix": dict(ndof=8, K=1.895663e-13, C_kappa=4.9335, is_line=False),
    "conservative quadratic": dict(ndof=8, K=2.2163e-13, C_kappa=5.2139, is_line=False),
    "origin-constrained line (B=0)": dict(ndof=9, K=5.872817e-15, C_kappa=None, is_line=True),
}

B_TESLA = 5.0
PT_MIN_GEV = 10.0
R_MIN_MM = (PT_MIN_GEV / (0.2998 * B_TESLA)) * 1000.0
KAPPA_MAX_TRUE = 1.0 / R_MIN_MM  # 1/mm, true scale


def P_true_at_sigma(name, sigma):
    m = models[name]
    ndof = m["ndof"]
    power = ndof / 2.0
    cut = stats.chi2.isf(P_3SIGMA, ndof)
    lam = LAMBDA_REF * (SIGMA_REF / sigma)
    chi2_part = m["K"] * (cut / lam ** 2) ** power
    if m["is_line"]:
        return chi2_part
    kappa_max_shrunk_ref = KAPPA_MAX_TRUE * LAMBDA_REF  # F_kappa independent of sigma
    F_kappa = m["C_kappa"] * kappa_max_shrunk_ref
    return F_kappa * chi2_part


def E_fakes(name, density, sigma=SIGMA_REF):
    p_true = P_true_at_sigma(name, sigma)
    return p_true * (density ** N_PLANES) * PROD_AREAS


# sanity check against Sections 11/12 at sigma=100um
print("Sanity check: P_true(sigma_ref) should match Sections 11/12 values")
for name in models:
    print(f"  {name:32s} P_true={P_true_at_sigma(name, SIGMA_REF):.4e}")

# --- find critical densities (E=1, E=10) at sigma_ref ---
print("\nCritical densities (hits/mm^2) at sigma=100um:")
crit_density = {}
for name in models:
    p_true = P_true_at_sigma(name, SIGMA_REF)
    rho1 = (1.0 / (p_true * PROD_AREAS)) ** (1.0 / N_PLANES)
    rho10 = (10.0 / (p_true * PROD_AREAS)) ** (1.0 / N_PLANES)
    crit_density[name] = (rho1, rho10)
    n1_nearest = rho1 * AREAS[0]
    n1_farthest = rho1 * AREAS[-1]
    print(f"  {name:32s} rho(E=1)={rho1:.4e} /mm^2  rho(E=10)={rho10:.4e} /mm^2  "
          f"[n_nearest={n1_nearest:.0f}, n_farthest={n1_farthest:.0f} hits at E=1]")

# --- find critical sigma (E=1) at a reference density ---
RHO_REF = crit_density["exact 4-param helix"][0]  # use exact model's own E(1) density as reference
print(f"\nUsing reference density = exact model's own E(1) density = {RHO_REF:.4e} /mm^2")
print("Critical resolution sigma (mm) at this fixed density:")
from scipy.optimize import brentq
crit_sigma = {}
for name in models:
    f = lambda s: E_fakes(name, RHO_REF, sigma=s) - 1.0
    # bracket
    lo, hi = 1e-4, 10.0
    try:
        s1 = brentq(f, lo, hi, xtol=1e-10, rtol=1e-12)
    except ValueError:
        s1 = np.nan
    f10 = lambda s: E_fakes(name, RHO_REF, sigma=s) - 10.0
    try:
        s10 = brentq(f10, lo, hi, xtol=1e-10, rtol=1e-12)
    except ValueError:
        s10 = np.nan
    crit_sigma[name] = (s1, s10)
    print(f"  {name:32s} sigma(E=1)={s1*1000:.2f} um  sigma(E=10)={s10*1000:.2f} um")


# --- Plot 1: E[#fakes] vs hit density, sigma fixed at 100um ---
colors = {"exact 4-param helix": "#1b4f8c", "conservative quadratic": "#c66a1f", "origin-constrained line (B=0)": "#5a5a5a"}
styles = {"exact 4-param helix": "-", "conservative quadratic": "--", "origin-constrained line (B=0)": "-."}

rho_grid = np.geomspace(2e-3, 6e-1, 400)
fig, ax = plt.subplots(figsize=(7.5, 5.5))
for name in models:
    E = np.array([E_fakes(name, r, sigma=SIGMA_REF) for r in rho_grid])
    ax.loglog(rho_grid, E, styles[name], color=colors[name], lw=2, label=name)
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.axhline(10.0, color="gray", lw=0.8, ls=":")
ax.set_xlabel(r"noise hit density $\rho$ (hits / mm$^2$, uniform across the tower)")
ax.set_ylabel(r"$E[\#\mathrm{fakes}]$")
ax.set_title(r"$E[\#\mathrm{fakes}]$ vs. hit density  ($B=5$ T, $p_T>10$ GeV/c, $\sigma=100\ \mu$m, 3$\sigma$ cut)")
ax.legend(loc="upper left", fontsize=9)
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig("section13_density.png", dpi=150)
print("saved section13_density.png")

# --- Plot 2: E[#fakes] vs resolution sigma, density fixed at RHO_REF ---
sigma_grid = np.geomspace(0.02, 0.5, 400)  # mm, i.e. 20um to 500um
fig2, ax2 = plt.subplots(figsize=(7.5, 5.5))
for name in models:
    E = np.array([E_fakes(name, RHO_REF, sigma=s) for s in sigma_grid])
    ax2.loglog(sigma_grid * 1000, E, styles[name], color=colors[name], lw=2, label=name)
ax2.axhline(1.0, color="gray", lw=0.8, ls=":")
ax2.axhline(10.0, color="gray", lw=0.8, ls=":")
ax2.set_xlabel(r"single-hit resolution $\sigma$ ($\mu$m)")
ax2.set_ylabel(r"$E[\#\mathrm{fakes}]$")
ax2.set_title(r"$E[\#\mathrm{fakes}]$ vs. resolution  ($B=5$ T, $p_T>10$ GeV/c, "
              r"$\rho=%.3g$ hits/mm$^2$, 3$\sigma$ cut)" % RHO_REF)
ax2.legend(loc="upper left", fontsize=9)
ax2.grid(True, which="both", alpha=0.25)
fig2.tight_layout()
fig2.savefig("section13_resolution.png", dpi=150)
print("saved section13_resolution.png")
