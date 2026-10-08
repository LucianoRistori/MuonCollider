"""
section14_real_tower.py

Plots for the new Section 14 worked example: the real-detector-informed
projective tower (IT barrel L0-L2 + OT barrel L0-L2 radii, 1 m^2
outermost plane, real non-uniform post-cut BIB hit densities).

Uses the calibration constants produced by real_tower_calibration.py
(copied here as literals -- the MC run is expensive, this script only
makes the two scan plots from the already-calibrated P_true values).

The six post-cut BIB densities are read from a run's per-layer table
(step4_cuts/density_per_layer_before_after_cuts.csv of ./run_all), by
default the paper's copy, docs/paper/data/section14_density_per_layer.csv:
    cd docs/paper/figures
    python3 ../../../fake_rate_framework/section14_real_tower.py [table.csv]
"""
import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_PLANES = 6
SIGMA_REF = 0.100  # mm: the resolution P_true was calibrated at
# The operating point's hit position resolution (the paper's Sections 14-16:
# 50 um), from the paper's smearing settings; P_true is scaled to it as
# (sigma/SIGMA_REF)^ndof (Sections 5.1, 13.2). Versions up to 1.6 used P_true
# at SIGMA_REF unscaled - 2^8 = 256 times too high for the helix models.
import configparser
_cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
_cp.read(Path(__file__).resolve().parent.parent / "docs" / "paper" / "settings" / "__smearing_config.txt")
SIGMA_OP = _cp.getfloat("position", "sigma_u_mm", fallback=SIGMA_REF)
assert abs(SIGMA_OP - _cp.getfloat("position", "sigma_v_mm", fallback=SIGMA_OP)) < 1e-12
print(f"operating position resolution {SIGMA_OP} mm (P_true calibrated at {SIGMA_REF} mm)")

Y_true = np.array([164.0, 354.0, 554.0, 819.0, 1153.0, 1486.0])
W_true = 1000.0 * Y_true / Y_true[-1]
TABLE = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path(__file__).resolve().parent.parent / "docs" / "paper" / "data" / "section14_density_per_layer.csv"
_rho = {}
with open(TABLE) as f:
    for r in csv.DictReader(f):
        if r["system"] in ("3", "5"):          # IT barrel, OT barrel
            _rho[(int(r["system"]), int(r["layer"]))] = float(r["bib_mean_density_after"])
RHO_REAL = np.array([_rho[(3, 0)], _rho[(3, 1)], _rho[(3, 2)], _rho[(5, 0)], _rho[(5, 1)], _rho[(5, 2)]])
print(f"densities from {TABLE}:", ", ".join(f"{x:.6g}" for x in RHO_REAL))
N_REAL = RHO_REAL * W_true ** 2
PROD_N_REAL = float(np.prod(N_REAL))

MODELS = {
    "exact helix":       dict(ndof=8, p_true=7.3054e-26, color="#1b9e77", ls="-"),
    "conservative quad": dict(ndof=8, p_true=6.2186e-26, color="#d95f02", ls="--"),
    "line (B=0)":        dict(ndof=9, p_true=3.6895e-28, color="#7570b3", ls=":"),
}
for m in MODELS.values():          # P_true at the operating resolution
    m["p_ref"] = m["p_true"]
    m["p_true"] = m["p_ref"] * (SIGMA_OP / SIGMA_REF) ** m["ndof"]

# --- Plot 1: E[#fakes] vs overall density multiplier mu (today's real density = mu=1) ---
mu = np.geomspace(0.1, 3000, 400)
fig, ax = plt.subplots(figsize=(6.0, 4.6))
for name, m in MODELS.items():
    efakes = m["p_true"] * PROD_N_REAL * mu ** N_PLANES
    ax.plot(mu, efakes, label=name, color=m["color"], ls=m["ls"], lw=2)
ax.axhline(1.0, color="gray", lw=0.8)
ax.axhline(10.0, color="gray", lw=0.8, ls="--")
ax.axvline(1.0, color="black", lw=1.0, alpha=0.6)
ax.text(1.08, 0.04, "today's\nestimated\nBIB density", fontsize=8, va="bottom",
        transform=ax.get_xaxis_transform())
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"overall density multiplier $\mu$ (BIB hit density relative to today's estimate)")
ax.set_ylabel(r"$E[\#\mathrm{fakes}]$ (6-layer tower)")
ax.set_title("Real-geometry tower: fake rate vs. background level")
ax.legend(fontsize=9, loc="upper left")
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig("section14_density.png", dpi=150)
print("Saved section14_density.png")

# --- Plot 2: E[#fakes] vs hit resolution sigma (density fixed at today's real value) ---
sigma = np.geomspace(0.02, 10.0, 400)  # mm
fig2, ax2 = plt.subplots(figsize=(6.0, 4.6))
for name, m in MODELS.items():
    efakes = m["p_ref"] * PROD_N_REAL * (sigma / SIGMA_REF) ** m["ndof"]
    ax2.plot(sigma, efakes, label=name, color=m["color"], ls=m["ls"], lw=2)
ax2.axhline(1.0, color="gray", lw=0.8)
ax2.axhline(10.0, color="gray", lw=0.8, ls="--")
ax2.axvline(SIGMA_OP, color="black", lw=1.0, alpha=0.6)
ax2.text(SIGMA_OP * 1.1, 0.04, f"today's\nassumed\nresolution\n({SIGMA_OP*1000:.0f} μm)", fontsize=8, va="bottom",
         transform=ax2.get_xaxis_transform())
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlabel(r"hit resolution $\sigma$ (mm)")
ax2.set_ylabel(r"$E[\#\mathrm{fakes}]$ (6-layer tower)")
ax2.set_title("Real-geometry tower: fake rate vs. hit resolution")
ax2.legend(fontsize=9, loc="upper left")
ax2.grid(True, which="both", alpha=0.25)
fig2.tight_layout()
fig2.savefig("section14_resolution.png", dpi=150)
print("Saved section14_resolution.png")

# --- Critical values table ---
print("\nCritical multipliers/resolutions for E[#fakes]=1 and =10:")
for name, m in MODELS.items():
    base = m["p_true"] * PROD_N_REAL
    mu1 = (1.0 / base) ** (1.0 / N_PLANES)
    mu10 = (10.0 / base) ** (1.0 / N_PLANES)
    sig1 = SIGMA_OP * (1.0 / base) ** (1.0 / m["ndof"])
    sig10 = SIGMA_OP * (10.0 / base) ** (1.0 / m["ndof"])
    print(f"{name:20s}  E@mu=1(today)={base:.3e}   "
          f"mu(E=1)={mu1:6.1f}x  mu(E=10)={mu10:6.1f}x   "
          f"sigma(E=1)={sig1:5.2f}mm ({sig1/SIGMA_OP:.0f}x)  sigma(E=10)={sig10:5.2f}mm   "
          f"P_true={m['p_true']:.3e}")
