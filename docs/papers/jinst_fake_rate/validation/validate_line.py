"""
validate_line.py - end-to-end validation of the fake-rate method for
straight-line tracks (B = 0), in small detectors where fakes are frequent
enough to be counted directly.

Detector: N planes at depths Y_i, each measuring two transverse
coordinates (u, v) with resolution sigma_i, square sensitive area of side
W_i. Track model: a free straight line in each view, u = a + b Y,
v = c + d Y (4 parameters, ndof = 2N - 4). Everything is done in whitened
coordinates (each coordinate divided by its resolution), where the fit is
an unweighted least-squares fit and chi2 = |Q u|^2 + |Q v|^2, with Q the
projector orthogonal to the track model.

For each configuration three numbers are compared, per event with n_i
noise hits on plane i:
  1. COUNTED: the mean number of hit combinations passing the chi2 cut,
     found by fitting ALL combinations in simulated noise-only events
     (and, for "at least k of N", the number of distinct fakes after
     merging candidates that share hits);
  2. DIRECT:  P_direct x prod n_i, with P_direct the passing probability
     of one random combination measured directly in the same detector
     (tests the counting identity E = P prod n_i);
  3. METHOD:  P_method x prod n_i, with P_method = K (c/lambda^2)^(ndof/2)
     and K calibrated by Monte Carlo on the detector shrunk by lambda
     (tests the small-chi2 power law and the rescaling).

Usage: python3 validate_line.py [quick]
Writes results_line.csv and prints a table.
"""
import itertools
import sys
import time

import numpy as np
from scipy.stats import chi2 as chi2dist

P_ONE_SIDED_3SIGMA = 1.3499e-3          # P(chi2 > cut) for genuine tracks


def cut_for(ndof):
    return float(chi2dist.isf(P_ONE_SIDED_3SIGMA, ndof))


# ---------------------------------------------------------------- geometry
class Geometry:
    def __init__(self, Y, W, sigma):
        self.Y = np.asarray(Y, float)            # plane depths (mm)
        self.W = np.asarray(W, float)            # plane sides (mm)
        self.sigma = np.asarray(sigma, float)    # resolution per coordinate (mm)
        self.w = self.W / self.sigma             # whitened plane sides
        A = np.stack([1.0 / self.sigma, self.Y / self.sigma], axis=1)   # whitened design
        self.Q = np.eye(len(self.Y)) - A @ np.linalg.solve(A.T @ A, A.T)
        self.ndof = 2 * (len(self.Y) - 2)

    def subset(self, idx):
        idx = list(idx)
        return Geometry(self.Y[idx], self.W[idx], self.sigma[idx])

    def shrunk(self, lam):
        """Plane sides divided by lam at fixed resolution (equivalently:
        every resolution multiplied by lam)."""
        return Geometry(self.Y, self.W / lam, self.sigma)

    def chi2(self, u, v):
        """u, v: whitened coordinates, shape (..., N)."""
        qu = u @ self.Q.T
        qv = v @ self.Q.T
        return np.einsum("...i,...i->...", qu, qu) + np.einsum("...i,...i->...", qv, qv)

    def random_hits(self, rng, shape):
        """Uniform noise hits in whitened coordinates, shape (*shape, N)."""
        h = self.w / 2.0
        u = rng.uniform(-h, h, size=tuple(shape) + (len(self.w),))
        v = rng.uniform(-h, h, size=tuple(shape) + (len(self.w),))
        return u, v


# ---------------------------------------------------------------- P direct
def p_direct(geo, cut, n_trials, rng, batch=2_000_000):
    passed, done = 0, 0
    while done < n_trials:
        m = min(batch, n_trials - done)
        u, v = geo.random_hits(rng, (m,))
        passed += int(np.count_nonzero(geo.chi2(u, v) < cut))
        done += m
    p = passed / n_trials
    return p, np.sqrt(passed) / n_trials if passed else np.nan, passed


# ---------------------------------------------------------------- the method
def calibrate_K(geo, lam, n_trials, rng, batch=2_000_000, n_bins=24, order=2):
    """Monte Carlo on the detector shrunk by lam. Near c' = 0 the residual
    density is even in r, so P(chi2 < c') = K c'^(ndof/2) (1 + beta c' + ...).
    K and beta are fitted (Poisson likelihood, binned in c'^(ndof/2)) over
    c' <= c_max = (smallest shrunk whitened plane)^2, and K is the
    c' -> 0 limit. Also returns the K_local curve for plots."""
    from scipy.optimize import minimize
    g = geo.shrunk(lam)
    d = g.ndof
    c_max = float(g.w.min() ** 2)
    # bins equally spaced in c'^(d/2), so each holds a similar number of entries
    edges = np.linspace(0.0, c_max ** (d / 2), n_bins + 1) ** (2.0 / d)
    grid = np.geomspace(c_max / 300, 4 * c_max, 30)
    obs = np.zeros(n_bins, np.int64)
    counts = np.zeros(len(grid), np.int64)
    done = 0
    while done < n_trials:                       # accumulate, never store
        m = min(batch, n_trials - done)
        u, v = g.random_hits(rng, (m,))
        c2 = np.sort(g.chi2(u, v))
        obs += np.diff(np.searchsorted(c2, edges))
        counts += np.searchsorted(c2, grid)
        done += m

    def model_cum(c, K, beta, gamma=0.0):
        return n_trials * K * c ** (d / 2) * (1 + beta * c + gamma * c * c)

    def nll(par):
        K, beta = np.exp(par[0]), par[1]
        gamma = par[2] if len(par) > 2 else 0.0
        mu = np.diff(model_cum(edges, K, beta, gamma))
        if np.any(mu <= 0):
            return 1e30
        return float(np.sum(mu - obs * np.log(mu)))
    K0 = max(obs.sum(), 1) / n_trials / c_max ** (d / 2)
    start = [np.log(K0), 0.0] + ([0.0] if order == 2 else [])
    res = minimize(nll, start, method="Nelder-Mead",
                   options=dict(xatol=1e-10, fatol=1e-10, maxiter=20000, maxfev=20000))
    res = minimize(nll, res.x, method="Nelder-Mead",
                   options=dict(xatol=1e-10, fatol=1e-10, maxiter=20000, maxfev=20000))
    K, beta = float(np.exp(res.x[0])), float(res.x[1])
    # statistical error of K: numerical Hessian of the likelihood
    x0 = res.x.copy()
    npar = len(x0)
    steps = np.array([1e-3] + [max(1e-3 * abs(v), 1e-6) for v in x0[1:]])
    H = np.zeros((npar, npar))
    f0 = nll(x0)
    for i in range(npar):
        for j in range(i, npar):
            ei = np.eye(npar)[i] * steps[i]
            ej = np.eye(npar)[j] * steps[j]
            H[i, j] = H[j, i] = (nll(x0 + ei + ej) - nll(x0 + ei - ej) - nll(x0 - ei + ej)
                                 + nll(x0 - ei - ej)) / (4 * steps[i] * steps[j])
    try:
        K_err = K * np.sqrt(max(np.linalg.inv(H)[0, 0], 0))
    except np.linalg.LinAlgError:
        K_err = np.nan
    return dict(grid=grid, counts=counts, k_local=counts / n_trials / grid ** (d / 2),
                K=K, K_err=K_err, beta=beta, gamma=(float(res.x[2]) if len(res.x) > 2 else 0.0),
                c_max=c_max, n_fit=int(obs.sum()))


def p_method(geo, cut, lam, n_trials, rng):
    cal = calibrate_K(geo, lam, n_trials, rng)
    d = geo.ndof
    P = cal["K"] * (cut / lam ** 2) ** (d / 2)
    return P, P * cal["K_err"] / cal["K"], cal


# ---------------------------------------------------------------- brute force
def count_all_combinations(geo, n_hits, cut, n_events, rng, k=None, merge_share=None):
    """Simulate noise-only events with n_hits[i] hits on plane i, fit every
    combination of one hit per plane (and, if k is given, every combination
    on every set of >= k planes, each with its own cut), and count those
    passing. Returns per-event arrays: passing combinations, and distinct
    fakes after merging candidates that share >= merge_share hits."""
    N = len(geo.Y)
    sets = [tuple(range(N))] if k is None else \
        [s for m in range(k, N + 1) for s in itertools.combinations(range(N), m)]
    subgeo = {s: geo.subset(s) for s in sets}
    cuts = {s: cut_for(subgeo[s].ndof) if k is not None else cut for s in sets}
    n_comb = np.zeros(n_events, int)
    n_distinct = np.zeros(n_events, int)
    for e in range(n_events):
        hits_u, hits_v = [], []
        for i in range(N):
            h = geo.w[i] / 2.0
            hits_u.append(rng.uniform(-h, h, n_hits[i]))
            hits_v.append(rng.uniform(-h, h, n_hits[i]))
        cands = []
        for s in sets:
            grids = np.meshgrid(*[np.arange(n_hits[i]) for i in s], indexing="ij")
            idx = np.stack([g.ravel() for g in grids], axis=1)        # (combos, |s|)
            u = np.stack([hits_u[i][idx[:, j]] for j, i in enumerate(s)], axis=1)
            v = np.stack([hits_v[i][idx[:, j]] for j, i in enumerate(s)], axis=1)
            ok = np.nonzero(subgeo[s].chi2(u, v) < cuts[s])[0]
            for r in ok:
                cands.append(frozenset((i, int(idx[r, j])) for j, i in enumerate(s)))
        n_comb[e] = len(cands)
        # distinct fakes: union-find over candidates sharing >= merge_share hits
        parent = list(range(len(cands)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        if merge_share is not None:
            for a in range(len(cands)):
                for b in range(a + 1, len(cands)):
                    if len(cands[a] & cands[b]) >= merge_share:
                        parent[find(a)] = find(b)
        n_distinct[e] = len({find(a) for a in range(len(cands))})
    return n_comb, n_distinct


# ---------------------------------------------------------------- configurations
def projective(N, Y0, Yout, W_out, sigma):
    Y = np.linspace(Y0, Yout, N)
    return Geometry(Y, W_out * Y / Yout, np.full(N, sigma))


def main():
    quick = len(sys.argv) > 1 and sys.argv[1] == "quick"
    rng = np.random.default_rng(20261009)
    n_ev = 200 if quick else 2000
    n_cal = 4_000_000 if quick else 40_000_000
    n_dir = 20_000_000 if quick else 200_000_000
    configs = [
        # name, geometry, hits per plane
        ("4 planes, w 40-160", projective(4, 100, 400, 16.0, 0.1), [12, 12, 12, 12]),
        ("6 planes, w 60-360", projective(6, 100, 600, 36.0, 0.1), [8] * 6),
        ("6 planes, w 30-180", projective(6, 100, 600, 18.0, 0.1), [6] * 6),
    ]
    rows = []
    for name, geo, n in configs:
        t0 = time.time()
        cut = cut_for(geo.ndof)
        prod_n = float(np.prod(n))
        lam = geo.w.min() / 2.0                  # smallest shrunk plane: 2 resolutions
        Pd, Pd_err, nd = p_direct(geo, cut, n_dir, rng)
        Pm, Pm_err, cal = p_method(geo, cut, lam, n_cal, rng)
        nc, _ = count_all_combinations(geo, n, cut, n_ev, rng)
        row = dict(config=name, N=len(geo.Y), ndof=geo.ndof, cut=cut, lam=lam,
                   hits=n[0], events=n_ev, counted=nc.mean(), counted_err=nc.std(ddof=1) / np.sqrt(n_ev),
                   direct=Pd * prod_n, direct_err=Pd_err * prod_n,
                   method=Pm * prod_n, method_err=Pm_err * prod_n, K=cal["K"], beta=cal["beta"], n_fit=cal["n_fit"])
        rows.append(row)
        print(f"{name:22s} ndof={geo.ndof} cut={cut:.2f} lam={lam:.2f}  "
              f"counted {row['counted']:.4g}+-{row['counted_err']:.2g}  "
              f"direct {row['direct']:.4g}+-{row['direct_err']:.2g}  "
              f"method {row['method']:.4g}+-{row['method_err']:.2g}   ({time.time()-t0:.0f}s)", flush=True)
    import csv
    with open("results_line.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
