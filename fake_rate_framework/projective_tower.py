"""
projective_tower.py

Concrete example requested by Luciano: 6 detector layers at 20cm (200mm)
intervals, layer sizes forming a projective tower with vertex at the
origin, farthest layer 1m x 1m. Sigma_x=sigma_y=100 micron on every
layer (uniform resolution; only the layer SIZE varies with z, following
the cone). Starting hit count: 10 hits/layer.

Because sigma is uniform across all layers, the same one-parameter
"conformal" scaling trick still applies even though the plane SIZES are
different from each other: chi2 = lambda^2 * chi2_shrunk exactly, if we
shrink EVERY layer's width by the same common factor lambda (preserving
the cone's proportions and each layer's true sigma). So we:

  1. Define the TRUE geometry (z_i, W_i = 1000*i/6 mm for i=1..6, sigma
     = 0.1mm for all).
  2. Shrink every W_i by a common factor lambda, chosen so the SMALLEST
     (nearest) layer's shrunk width is a modest, tractable multiple of
     sigma (~10 sigma, matching what worked well before).
  3. Run high-statistics vectorized MC on the SHRUNK geometry (tracksim
     already supports per-plane W/sigma, no code changes needed).
  4. Fit the deep-tail power law there (fixed exponent = ndof/2 = 4 for
     N=6 planes), get K_shrunk.
  5. Rescale to the TRUE geometry: P_true(chi2<cut) = P_shrunk(chi2 <
     cut/lambda^2) = K_shrunk * (cut/lambda^2)^4.
  6. Combine with the combinatorial factor (10 hits/layer, 6 layers) to
     get the expected number of fake tracks vs chi2 cut, and show how it
     scales with an overall hit-density multiplier mu (expect mu^6).
"""

import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tracksim import Plane, Detector, batch_random_chi2

N_PLANES = 6
SIGMA = 0.100  # mm, same for every layer
Z_SPACING = 200.0  # mm

# True geometry
z_true = np.array([Z_SPACING * i for i in range(1, N_PLANES + 1)])
W_true = np.array([1000.0 * i / N_PLANES for i in range(1, N_PLANES + 1)])  # mm

NDOF = 2 * N_PLANES - 4  # = 8
POWER = NDOF / 2.0       # = 4

# Choose shrink factor so the SMALLEST (nearest) layer's shrunk width is
# ~10 sigma (the regime that worked well in earlier calibrations).
TARGET_RATIO_NEAREST = 10.0
LAMBDA = W_true[0] / (TARGET_RATIO_NEAREST * SIGMA)
W_shrunk = W_true / LAMBDA

print("True geometry:")
for i in range(N_PLANES):
    print(f"  layer {i+1}: z={z_true[i]:8.1f}mm  W={W_true[i]:9.3f}mm  W/sigma={W_true[i]/SIGMA:9.1f}")
print(f"\nShrink factor lambda = {LAMBDA:.4f}")
print("Shrunk geometry (used for direct MC):")
for i in range(N_PLANES):
    print(f"  layer {i+1}: z={z_true[i]:8.1f}mm  W={W_shrunk[i]:7.4f}mm  W/sigma={W_shrunk[i]/SIGMA:7.2f}")

detector_shrunk = Detector(planes=[
    Plane(z=z_true[i], size_x=W_shrunk[i], size_y=W_shrunk[i], sigma_x=SIGMA, sigma_y=SIGMA)
    for i in range(N_PLANES)
])

# ---------------------------------------------------------------------------
# High-statistics MC on the shrunk geometry
# ---------------------------------------------------------------------------
TOTAL_TRIALS = 150_000_000
CHUNK = 5_000_000
SEED = 424242

rng = np.random.default_rng(SEED)
n_chunks = TOTAL_TRIALS // CHUNK
all_chi2 = np.empty(TOTAL_TRIALS, dtype=np.float64)

t0 = time.time()
for i in range(n_chunks):
    all_chi2[i*CHUNK:(i+1)*CHUNK] = batch_random_chi2(detector_shrunk, CHUNK, rng)
print(f"\nGenerated {TOTAL_TRIALS:,} trials in {time.time()-t0:.0f}s")

sorted_chi2 = np.sort(all_chi2)
n = len(sorted_chi2)
del all_chi2
print(f"min chi2 observed (shrunk): {sorted_chi2[0]:.4f}   (ndof={NDOF})")

# ---------------------------------------------------------------------------
# Fixed-exponent (ndof/2=4) plateau fit for K_shrunk
# ---------------------------------------------------------------------------
event_counts = np.unique(np.geomspace(20, int(0.02 * n), 60).astype(int))
cuts = sorted_chi2[event_counts]
p_emp = event_counts / n
K_local = p_emp / cuts ** POWER
rel_err = 1.0 / np.sqrt(event_counts)

print(f"\n{'n_events':>10} {'cut(chi2)':>12} {'P_emp':>12} {'K_local':>14} {'rel.err':>8}")
for ne, c, p, k, re in zip(event_counts, cuts, p_emp, K_local, rel_err):
    print(f"{ne:>10d} {c:>12.4f} {p:>12.3e} {k:>14.4e} {re*100:>7.1f}%")

# Use the deepest quarter of the sampled range as the plateau region
deep = cuts < np.percentile(cuts, 25)
w = event_counts[deep].astype(float)
K_shrunk = np.sum(K_local[deep] * w) / np.sum(w)
print(f"\nPlateau K_shrunk (deep quartile, cut<{np.percentile(cuts,25):.2f}) = {K_shrunk:.4e}")

# also report shallower quartile for comparison (trend check)
deep2 = (cuts >= np.percentile(cuts, 25)) & (cuts < np.percentile(cuts, 50))
w2 = event_counts[deep2].astype(float)
K_shrunk_2 = np.sum(K_local[deep2] * w2) / np.sum(w2)
print(f"Next quartile K_shrunk (cut in [{np.percentile(cuts,25):.2f},{np.percentile(cuts,50):.2f}]) = {K_shrunk_2:.4e}"
      f"  (ratio to plateau: {K_shrunk_2/K_shrunk:.3f})")

# ---------------------------------------------------------------------------
# Rescale to TRUE geometry and combine with combinatorics
# ---------------------------------------------------------------------------
N_HITS_PER_LAYER = 10
C = N_HITS_PER_LAYER ** N_PLANES

def P_true(cut):
    cut_shrunk = cut / LAMBDA**2
    return K_shrunk * cut_shrunk ** POWER

print(f"\nCombinatorial factor C = {N_HITS_PER_LAYER}^{N_PLANES} = {C:,}")
print(f"\n{'chi2 cut':>10} {'P_true(chi2<cut)':>18} {'E[#fakes] (10 hits/layer)':>26}")
for cut in [8, 10, 15, 20, 25, 30, 40, 50]:
    p = P_true(cut)
    efakes = C * p
    print(f"{cut:>10.1f} {p:>18.3e} {efakes:>26.3e}")

print("\nScaling with an overall hit-density multiplier mu (E[#fakes] ~ mu^N):")
cut_example = 20.0
p20 = P_true(cut_example)
for mu in [0.5, 1, 2, 5, 10]:
    efakes = C * (mu ** N_PLANES) * p20
    print(f"  mu={mu:>4.1f}  ->  E[#fakes at cut={cut_example:.0f}] = {efakes:.3e}   (n_hits/layer={N_HITS_PER_LAYER*mu:.1f})")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

ax = axes[0]
ax.errorbar(cuts, K_local, yerr=K_local*rel_err, fmt="o", ms=4, color="#4C72B0",
            ecolor="#4C72B0", alpha=0.7, capsize=2)
ax.axhline(K_shrunk, color="r", ls="--", lw=1.5, label=f"plateau K={K_shrunk:.3e}")
ax.set_xscale("log")
ax.set_xlabel(r"$\chi^2$ cut (shrunk geometry)")
ax.set_ylabel(r"implied $K=P(\chi^2<\mathrm{cut})/\mathrm{cut}^4$")
ax.set_title(f"Projective tower (N={N_PLANES}, ndof={NDOF}):\nfixed-exponent plateau check")
ax.legend()
ax.grid(True, which="both", alpha=0.3)

ax2 = axes[1]
cuts_plot = np.linspace(5, 50, 200)
efakes_plot = [C * P_true(c) for c in cuts_plot]
ax2.semilogy(cuts_plot, efakes_plot, color="#4C72B0", lw=2)
ax2.axvline(NDOF, color="k", ls=":", alpha=0.6, label=f"ndof={NDOF}")
ax2.set_xlabel(r"$\chi^2$ cut")
ax2.set_ylabel("E[# fake tracks] (10 hits/layer)")
ax2.set_title("Expected fake tracks vs chi2 cut\n(true 6-layer projective tower geometry)")
ax2.legend()
ax2.grid(True, which="both", alpha=0.3)

fig.tight_layout()
fig.savefig("projective_tower_result.png", dpi=150)
print("\nSaved plot to projective_tower_result.png")
