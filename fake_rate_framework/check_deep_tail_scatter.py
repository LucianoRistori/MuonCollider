"""
check_deep_tail_scatter.py

Luciano asked: is the low K_local value seen at the very deepest chi2
tail (in the projective-tower plateau check) just statistical noise, or
a real systematic effect?

Test: run many INDEPENDENT replicas of the same (shrunk) projective-
tower MC, and for each replica compute K_local = P(chi2<cut)/cut^4 at a
FIXED, shallow order statistic (e.g. the 20th-smallest chi2 out of each
replica's trials). If the deep-tail dip were a real systematic, the
replica-to-replica values should cluster tightly (well below the
plateau, small scatter). If it's just Poisson counting noise, the
replica values should scatter with relative width ~1/sqrt(n_events)
(e.g. ~22% for n_events=20), consistent with the error bars already
used in the plateau plot.
"""

import time
import numpy as np
from tracksim import Plane, Detector, batch_random_chi2

N_PLANES = 6
SIGMA = 0.100
Z_SPACING = 200.0
NDOF = 2 * N_PLANES - 4
POWER = NDOF / 2.0

z_true = np.array([Z_SPACING * i for i in range(1, N_PLANES + 1)])
W_true = np.array([1000.0 * i / N_PLANES for i in range(1, N_PLANES + 1)])
LAMBDA = W_true[0] / (10.0 * SIGMA)
W_shrunk = W_true / LAMBDA

detector_shrunk = Detector(planes=[
    Plane(z=z_true[i], size_x=W_shrunk[i], size_y=W_shrunk[i], sigma_x=SIGMA, sigma_y=SIGMA)
    for i in range(N_PLANES)
])

N_REPLICAS = 25
TRIALS_PER_REPLICA = 10_000_000
N_ORDER = 20  # look at the 20th-smallest chi2 in each replica (matches the noisiest point before)

K_values = []
cuts_used = []
t0 = time.time()
for rep in range(N_REPLICAS):
    rng = np.random.default_rng(1000 + rep)  # different seed per replica
    chi2 = batch_random_chi2(detector_shrunk, TRIALS_PER_REPLICA, rng)
    chi2.sort()
    cut = chi2[N_ORDER - 1]  # 20th smallest (0-indexed: index 19)
    P = N_ORDER / TRIALS_PER_REPLICA
    K = P / cut ** POWER
    K_values.append(K)
    cuts_used.append(cut)
    del chi2

print(f"Total time: {time.time()-t0:.0f}s for {N_REPLICAS} replicas x {TRIALS_PER_REPLICA:,} trials")

K_values = np.array(K_values)
cuts_used = np.array(cuts_used)

print(f"\ncut values (20th-smallest chi2) across replicas: min={cuts_used.min():.3f} max={cuts_used.max():.3f}")
print(f"\nK_local values across {N_REPLICAS} independent replicas (n_events={N_ORDER} each):")
for i, (k, c) in enumerate(zip(K_values, cuts_used)):
    print(f"  replica {i:2d}: cut={c:8.4f}  K={k:.4e}")

mean_K = K_values.mean()
std_K = K_values.std(ddof=1)
rel_scatter = std_K / mean_K
theory_rel_scatter = 1.0 / np.sqrt(N_ORDER)

print(f"\nMean K across replicas   = {mean_K:.4e}")
print(f"Std  K across replicas   = {std_K:.4e}")
print(f"Observed relative scatter (std/mean) = {rel_scatter*100:.1f}%")
print(f"Theoretical (Poisson, 1/sqrt(n={N_ORDER})) = {theory_rel_scatter*100:.1f}%")
print(f"Ratio observed/theory = {rel_scatter/theory_rel_scatter:.2f}  "
      f"(close to 1.0 => consistent with pure statistics)")

# Also: min/max range check -- does the earlier single-run low value
# (~8.4e-12 from the 150M run) fall within the range we see here across
# independent replicas?
print(f"\nRange of K across replicas: [{K_values.min():.3e}, {K_values.max():.3e}]")
print("For comparison, the original 150M-trial run's n_events=20 point was K=8.38e-12,")
print("and the fitted plateau was K=1.13e-11.")
