"""
worked_example.py - paper 1, Section 6: the six-layer barrel tower inspired by a
Muon Collider tracker, with rounded radii and densities.

Expected fakes on at least 5 of 6 layers: E = sum_S P_S prod_{i in S} n_i over the
6-of-6 set and the six 5-of-6 sets, with P_S = K_S c^(ndof/2) and K_S from the
volume of the accepted track manifold (fake_rate_framework/geometric_K.py).
Then the window occupancy eta of each set (eta.py) and the distinct-fake estimate.

Usage: python3 worked_example.py      -> worked_example.json, printed table
"""
import itertools, json, math, sys
from pathlib import Path
import numpy as np
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "fake_rate_framework"))
from geometric_K import K_geom, P_3SIGMA

Y = np.array([160.0, 350.0, 550.0, 820.0, 1200.0, 1500.0])        # mm, two significant figures
W = 1000.0 * Y / Y[-1]                                              # projective, 1 m^2 outermost plane
RHO = np.array([0.039, 0.011, 0.0036, 0.0013, 0.00083, 0.00053])    # hits/mm^2 after cuts (100 um, 100 ps)
SIGMA = 0.100                                                       # mm, both coordinates
B, PT = 5.0, 10.0
KMAX = 0.2998 * B / PT / 1000.0


def pred_var(Yo, Yi, basis):
    A = np.stack([f(Yo) for f in basis], 1) / SIGMA
    a = np.array([f(np.array([Yi]))[0] for f in basis])
    return a @ np.linalg.inv(A.T @ A) @ a


BEND = [lambda y: y, lambda y: 0.5 * y ** 2]          # X = b Y + (kappa/2) Y^2, linearized at kappa = 0
DEPTH = [lambda y: np.ones_like(y), lambda y: y]      # Z = p + q r


def main():
    n = RHO * W ** 2
    out = dict(Y=Y.tolist(), W=W.tolist(), rho=RHO.tolist(), n=n.tolist(), sigma=SIGMA, kappa_max=KMAX, sets={})
    E_sum = E_dist = 0.0
    print(f"hits per plane n_i: {', '.join(f'{x:.0f}' for x in n)}")
    for S in [tuple(range(6))] + list(itertools.combinations(range(6), 5)):
        S = list(S); d = 2 * len(S) - 4; c = float(chi2dist.isf(P_3SIGMA, d))
        K = K_geom(Y[S], W[S], SIGMA, KMAX); P = K * c ** (d / 2); E = P * float(np.prod(n[S]))
        eta = 0.0
        for i in S:
            o = [j for j in S if j != i]
            Vb, Vd = pred_var(Y[o], Y[i], BEND), pred_var(Y[o], Y[i], DEPTH)
            eta += (n[i] - 1) * math.pi * (c / 2) * math.sqrt((SIGMA ** 2 + Vb) * (SIGMA ** 2 + Vd)) / W[i] ** 2
        out["sets"][",".join(map(str, S))] = dict(ndof=d, cut=c, K=K, P=P, E=E, eta=eta)
        E_sum += E; E_dist += E * math.exp(-eta / 2)
        miss = sorted(set(range(6)) - set(S))
        print(f"{'all 6' if not miss else 'without layer %d' % (miss[0] + 1):18s} ndof {d}  P {P:.3e}  E {E:.3e}  eta {eta:.3f}", flush=True)
    E6 = out["sets"]["0,1,2,3,4,5"]["E"]
    mu1 = math.exp(next(x for x in np.linspace(-3, 5, 80001) if E6 * math.exp(6 * x) + (E_sum - E6) * math.exp(5 * x) >= 1))
    out.update(E_sum=E_sum, E_distinct=E_dist, mu_one_fake=mu1)
    print(f"at least 5 of 6: E = {E_sum:.3e};  distinct ~ {E_dist:.3e};  one fake at {mu1:.2f}x the densities")
    json.dump(out, open(Path(__file__).with_name("worked_example.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
