"""How stable is K? One large Monte Carlo of the calibration detector
(5-plane subset 0-4 of the reference geometry, ndof 6), then K fitted with
P = K c^(d/2) (1 + b c + g c^2 [+ h c^3]) over different ranges c' <= c_max."""
import sys, numpy as np
from scipy.optimize import minimize
from validate_line import projective, cut_for
n_tr = int(float(sys.argv[1])) if len(sys.argv) > 1 else 300_000_000
rng = np.random.default_rng(99)
g = projective(6, 100, 600, 2 * 0.1 * 6, 0.1).subset((0, 1, 2, 3, 4))
d = g.ndof
edges = np.concatenate([[0.0], np.geomspace(1e-2, 4.0, 300)])
obs = np.zeros(len(edges) - 1, np.int64); done = 0
while done < n_tr:
    u, v = g.random_hits(rng, (2_000_000,)); c = np.sort(g.chi2(u, v))
    obs += np.diff(np.searchsorted(c, edges)); done += 2_000_000
print(f"{n_tr:.1e} trials, {obs.sum()} below c'=4")
cum = np.cumsum(obs)
for c in (0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
    j = int(np.searchsorted(edges, c)) - 1
    n = cum[j]; ce = edges[j + 1]
    print(f"  K_local(c'={ce:.3f}) = {n / n_tr / ce ** (d / 2):.6g} +- {np.sqrt(n) / n_tr / ce ** (d / 2):.2g}  ({n})")
for cmax in (0.5, 1.0, 2.0, 4.0):
    nb = int(np.searchsorted(edges, cmax)); e = edges[:nb + 1]; o = obs[:nb]
    for order in (1, 2, 3):
        def nll(p):
            K = np.exp(p[0]); co = list(p[1:]) + [0] * (3 - len(p[1:]))
            mu = np.diff(n_tr * K * e ** (d / 2) * (1 + co[0] * e + co[1] * e**2 + co[2] * e**3))
            return 1e30 if np.any(mu <= 0) else float(np.sum(mu - o * np.log(mu)))
        x = [np.log(o.sum() / n_tr / cmax ** (d / 2))] + [0.0] * order
        for _ in range(3):
            x = minimize(nll, x, method="Nelder-Mead", options=dict(xatol=1e-11, fatol=1e-11, maxiter=40000, maxfev=40000)).x
        print(f"  c_max {cmax:3.1f}  order {order}:  K = {np.exp(x[0]):.6g}   (entries {o.sum()})")
