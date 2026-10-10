"""
helix_validation.py - end-to-end validation of the fake-rate method for
the exact helix fit (tracks from the beam line, pT acceptance), the
counterpart of the straight-line tests in validate_line.py.

  calib  <N> <trials>              K for the reference tower (w_min = 2)
  extrap <N> <w_min list> <max>    test 1: P direct vs P method for towers
                                   of growing size (exact rescalings)
  count  <w_min> <hits> <events>   test 2: noise-only events, every
                                   combination on 6 planes and on every
                                   5-plane subset fitted; counted vs
                                   P_direct x prod n; distinct fakes; eta

Acceptance throughout: chi2 < one-sided-3-sigma cut for the ndof AND
|kappa_fit| <= xi / Y_N with xi = 0.2.
"""
import csv
import itertools
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2dist

from helix_core import Tower

P3 = 1.3499e-3


def cut_for(ndof):
    return float(chi2dist.isf(P3, ndof))


def calibrate(tower, n_trials, rng, batch=1_000_000, n_bins=24):
    d = tower.ndof
    c_max = float((tower.W.min() / tower.sigma) ** 2)
    edges = np.linspace(0.0, c_max ** (d / 2), n_bins + 1) ** (2.0 / d)
    grid = np.geomspace(c_max / 300, c_max, 25)
    obs = np.zeros(n_bins, np.int64)
    counts = np.zeros(len(grid), np.int64)
    done = 0
    while done < n_trials:
        m = min(batch, n_trials - done)
        X, Z = tower.random_hits(rng, m)
        c2, k = tower.chi2(X, Z, c_max)
        c2 = np.sort(c2[np.abs(k) <= tower.kappa_max])
        obs += np.diff(np.searchsorted(c2, edges))
        counts += np.searchsorted(c2, grid)
        done += m

    def nll(p):
        K = np.exp(p[0])
        mu = np.diff(n_trials * K * edges ** (d / 2) * (1 + p[1] * edges + p[2] * edges ** 2))
        return 1e30 if np.any(mu <= 0) else float(np.sum(mu - obs * np.log(mu)))
    x = [np.log(max(obs.sum(), 1) / n_trials / c_max ** (d / 2)), 0.0, 0.0]
    for _ in range(3):
        x = minimize(nll, x, method="Nelder-Mead",
                     options=dict(xatol=1e-11, fatol=1e-11, maxiter=40000, maxfev=40000)).x
    steps = np.array([1e-3, max(1e-3 * abs(x[1]), 1e-6), max(1e-3 * abs(x[2]), 1e-7)])
    H = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            ei, ej = np.eye(3)[i] * steps[i], np.eye(3)[j] * steps[j]
            H[i, j] = H[j, i] = (nll(x + ei + ej) - nll(x + ei - ej) - nll(x - ei + ej)
                                 + nll(x - ei - ej)) / (4 * steps[i] * steps[j])
    K = float(np.exp(x[0]))
    K_err = K * float(np.sqrt(max(np.linalg.inv(H)[0, 0], 0)))
    return dict(K=K, K_err=K_err, beta=float(x[1]), gamma=float(x[2]), c_max=c_max,
                n_fit=int(obs.sum()), grid=grid.tolist(),
                k_local=(counts / n_trials / grid ** (d / 2)).tolist(), counts=counts.tolist())


def p_direct(tower, cut, n_max, rng, want=400, batch=1_000_000, min_trials=2_000_000):
    passed, done = 0, 0
    while done < n_max and (passed < want or done < min_trials):
        m = min(batch, n_max - done)
        X, Z = tower.random_hits(rng, m)
        ok, _, _ = tower.accepted(X, Z, cut)
        passed += int(ok.sum())
        done += m
    return passed / done, np.sqrt(passed) / done, passed, done


def cmd_calib(N, n_trials):
    rng = np.random.default_rng(5000 + N)
    t0 = time.time()
    cal = calibrate(Tower(N, 2.0), n_trials, rng)
    cal.update(N=N, ndof=2 * N - 4, cut=cut_for(2 * N - 4), trials=n_trials, seconds=time.time() - t0)
    json.dump(cal, open(f"helix_calib_N{N}.json", "w"), indent=1)
    print(f"N={N}: K={cal['K']:.5g} +- {cal['K_err']:.2g}  beta={cal['beta']:.4g} gamma={cal['gamma']:.3g} "
          f"({cal['n_fit']} accepted below c'={cal['c_max']:g}; {cal['seconds']:.0f}s)", flush=True)


def cmd_extrap(N, w_list, n_max):
    cal = json.load(open(f"helix_calib_N{N}.json"))
    rng = np.random.default_rng(6000 + N)
    d, cut = 2 * N - 4, cut_for(2 * N - 4)
    rows = []
    for wm in w_list:
        lam = wm / 2.0
        cp = cut / lam ** 2
        Pm = cal["K"] * cp ** (d / 2)
        t0 = time.time()
        Pd, Pde, npass, ntr = p_direct(Tower(N, wm), cut, n_max, rng)
        r = dict(N=N, ndof=d, w_min=wm, c_prime=cp, P_method=Pm, P_method_err=Pm * cal["K_err"] / cal["K"],
                 P_direct=Pd, P_direct_err=Pde, n_pass=npass, n_trials=ntr,
                 ratio=Pd / Pm, ratio_err=Pde / Pm, K=cal["K"], beta=cal["beta"], gamma=cal["gamma"])
        rows.append(r)
        print(f"  w_min {wm:6.1f} c'={cp:8.4g}  direct {Pd:.4e} +- {Pde:.1e} ({npass} of {ntr:.2e})  "
              f"method {Pm:.4e}  ratio {r['ratio']:.4f} +- {r['ratio_err']:.4f}  ({time.time()-t0:.0f}s)", flush=True)
        with open(f"helix_extrapolation_N{N}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


def cmd_count(w_min, n, n_ev, n_p=4_000_000, seed=0):
    N = 6
    rng = np.random.default_rng(int(100 * w_min) + n + 777 + 1000003 * seed)
    sets = [tuple(range(N))] + list(itertools.combinations(range(N), 5))
    tow = {s: Tower(N, w_min, planes=s) for s in sets}
    cut = {s: cut_for(tow[s].ndof) for s in sets}
    full = Tower(N, w_min)
    t0 = time.time()
    P = {s: p_direct(tow[s], cut[s], n_p, rng, want=4000)[:2] for s in sets}
    print(f"helix count: w_min {w_min}, {n} hits/plane: direct P done ({time.time()-t0:.0f}s)", flush=True)
    idx = {s: np.stack([g.ravel() for g in np.meshgrid(*[np.arange(n)] * len(s), indexing="ij")], 1)
           for s in sets}
    count = {s: np.zeros(n_ev, int) for s in sets}
    distinct = np.zeros(n_ev, int)
    eta = []
    t0 = time.time()
    for e in range(n_ev):
        hx = np.stack([rng.uniform(-0.5, 0.5, n) * W for W in full.W])
        hz = np.stack([rng.uniform(-0.5, 0.5, n) * W for W in full.W])
        cands = []
        for s in sets:
            ix = idx[s]
            X = np.stack([hx[i][ix[:, j]] for j, i in enumerate(s)], 1)
            Z = np.stack([hz[i][ix[:, j]] for j, i in enumerate(s)], 1)
            ok, _, _ = tow[s].accepted(X, Z, cut[s])
            rows = np.nonzero(ok)[0]
            count[s][e] = len(rows)
            for r in rows:
                cands.append((s, frozenset((i, int(ix[r, j])) for j, i in enumerate(s))))
        # window occupancy: other hits that could replace one of a candidate's
        for s, c in cands:
            hit = dict(c)
            tot = 0
            for i in s:
                X = np.tile([hx[p][hit[p]] for p in s], (n, 1))
                Z = np.tile([hz[p][hit[p]] for p in s], (n, 1))
                j = s.index(i)
                X[:, j], Z[:, j] = hx[i], hz[i]
                ok, _, _ = tow[s].accepted(X, Z, cut[s])
                tot += int(ok.sum()) - 1
            eta.append(tot)
        parent = list(range(len(cands)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        for a in range(len(cands)):
            for b in range(a + 1, len(cands)):
                if len(cands[a][1] & cands[b][1]) >= 4:
                    parent[find(a)] = find(b)
        distinct[e] = len({find(a) for a in range(len(cands))})
        if (e + 1) % max(1, n_ev // 10) == 0:
            print(f"  {e+1} events ({time.time()-t0:.0f}s)", flush=True)
    se = lambda x: x.std(ddof=1) / np.sqrt(n_ev)
    out = []
    tot = np.zeros(n_ev, int)
    pred, pvar = 0.0, 0.0
    for s in sets:
        f = float(n) ** len(s)
        label = "all 6" if len(s) == 6 else "no plane " + str(sorted(set(range(N)) - set(s))[0] + 1)
        out.append(dict(set=label, counted=count[s].mean(), counted_err=se(count[s]),
                        predicted=P[s][0] * f, predicted_err=P[s][1] * f))
        tot += count[s]
        pred += P[s][0] * f
        pvar += (P[s][1] * f) ** 2
    out.append(dict(set="sum", counted=tot.mean(), counted_err=se(tot), predicted=pred, predicted_err=np.sqrt(pvar)))
    out.append(dict(set="distinct (>= 4 shared hits)", counted=distinct.mean(), counted_err=se(distinct),
                    predicted=np.nan, predicted_err=np.nan))
    out.append(dict(set="eta", counted=float(np.mean(eta)) if eta else np.nan, counted_err=np.nan,
                    predicted=np.nan, predicted_err=np.nan))
    for r in out:
        r.update(w_min=w_min, hits=n, events=n_ev)
        print(f"  {r['set']:28s} counted {r['counted']:.4g} +- {r['counted_err']:.2g}   predicted {r['predicted']:.4g} +- {r['predicted_err']:.2g}")
    with open(f"helix_count_w{w_min:g}_n{n}" + (f"_s{seed}" if seed else "") + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "calib":
        cmd_calib(int(sys.argv[2]), int(float(sys.argv[3])))
    elif mode == "extrap":
        cmd_extrap(int(sys.argv[2]), [float(x) for x in sys.argv[3].split(",")], int(float(sys.argv[4])))
    elif mode == "count":
        cmd_count(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), seed=int(sys.argv[5]) if len(sys.argv) > 5 else 0)
