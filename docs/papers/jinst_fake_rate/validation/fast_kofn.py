"""
fast_kofn.py - "at least 5 of 6 planes" with realistic numbers of hits.

For each simulated noise-only event (n hits on each of 6 projective
planes) every passing combination is found exactly with fastcount.search:
on all 6 planes (ndof 8) and on each 5-plane subset (ndof 6), each with
its own 3-sigma cut. Recorded per configuration:
  - the passing combinations per set, and their sum;
  - the METHOD prediction P_S prod n for each set, with K_S from
    calib_sets.json (an end-to-end test of the method deep in the
    extrapolation: P_S down to ~1e-15);
  - distinct fakes, after merging candidates sharing >= 4 (also >= 3) hits;
  - the window occupancy eta: for a passing candidate, the number of OTHER
    hits, summed over its planes, that could replace its own hit on that
    plane and still pass. Distinct fakes -> sum of combinations as eta -> 0.
The plane sides grow with n so that the expected sum stays around 2-3 per
event: w_min = 15 (n/8)^(5/6) resolutions.

Usage: python3 fast_kofn.py <hits per plane> <events>   (appends to kofn_results.csv)
"""
import csv
import itertools
import json
import os
import sys
import time

import numpy as np

from fastcount import Plane, search
from validate_line import cut_for, projective


def fit_rest(planes, u, v):
    """Straight-line fit (both views) to the given points; returns chi2 and
    a function giving (u_hat, v_hat, V) at depth Y."""
    w = np.array([1 / p.sigma ** 2 for p in planes])
    Y = np.array([p.Y for p in planes])
    S0, S1, S2 = w.sum(), (w * Y).sum(), (w * Y * Y).sum()
    D = S0 * S2 - S1 ** 2
    out = []
    chi2 = 0.0
    for x in (u, v):
        Sx, SxY = (w * x).sum(), (w * Y * x).sum()
        b = (S0 * SxY - S1 * Sx) / D
        a = (Sx - b * S1) / S0
        chi2 += float((w * (x - a - b * Y) ** 2).sum())
        out.append((a, b))

    def pred(Yk):
        V = (S2 - 2 * Yk * S1 + Yk ** 2 * S0) / D
        return out[0][0] + out[0][1] * Yk, out[1][0] + out[1][1] * Yk, V
    return chi2, pred


def main():
    n = int(sys.argv[1])
    n_ev = int(sys.argv[2])
    N, sigma = 6, 0.1
    w_min = 15.0 * (n / 8.0) ** (5.0 / 6.0)
    geo = projective(N, 100, 100 * N, w_min * sigma * N, sigma)
    planes = [Plane(y, W, s) for y, W, s in zip(geo.Y, geo.W, geo.sigma)]
    cal = json.load(open("calib_sets.json"))
    lam = w_min / 2.0                       # calibration geometry has w_min = 2
    sets = [tuple(range(N))] + list(itertools.combinations(range(N), 5))
    pred, pred_err = {}, {}
    for s in sets:
        c = cal[",".join(map(str, s))]
        P = c["K"] * (c["cut"] / lam ** 2) ** (c["ndof"] / 2)
        pred[s] = P * float(n) ** len(s)
        pred_err[s] = pred[s] * c["K_err"] / c["K"]
    rng = np.random.default_rng(31 * n + 7 + (1000003 * int(sys.argv[3]) if len(sys.argv) > 3 else 0))
    count = {s: np.zeros(n_ev, int) for s in sets}
    distinct = {4: np.zeros(n_ev, int), 3: np.zeros(n_ev, int)}
    eta = []
    t0 = time.time()
    for e in range(n_ev):
        hu = [rng.uniform(-W / 2, W / 2, n) for W in geo.W]
        hv = [rng.uniform(-W / 2, W / 2, n) for W in geo.W]
        cands = []
        for s in sets:
            res, _ = search([planes[i] for i in s], [(hu[i], hv[i]) for i in s], cut_for(2 * len(s) - 4))
            count[s][e] = len(res)
            for r in res:
                cands.append((s, frozenset((i, int(r[j])) for j, i in enumerate(s))))
        # window occupancy of each passing candidate
        for s, c in cands:
            cut = cut_for(2 * len(s) - 4)
            hit = dict(c)
            tot = 0
            for i in s:
                rest = [p for p in s if p != i]
                chi2_rest, pr = fit_rest([planes[p] for p in rest],
                                         np.array([hu[p][hit[p]] for p in rest]),
                                         np.array([hv[p][hit[p]] for p in rest]))
                uh, vh, V = pr(planes[i].Y)
                d2 = ((hu[i] - uh) ** 2 + (hv[i] - vh) ** 2) / (planes[i].sigma ** 2 + V)
                tot += int(np.count_nonzero(chi2_rest + d2 < cut)) - 1
            eta.append(tot)
        for m in (4, 3):
            parent = list(range(len(cands)))

            def find(a):
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a
            for a in range(len(cands)):
                for b in range(a + 1, len(cands)):
                    if len(cands[a][1] & cands[b][1]) >= m:
                        parent[find(a)] = find(b)
            distinct[m][e] = len({find(a) for a in range(len(cands))})
    total = sum(count[s] for s in sets)
    se = lambda x: x.std(ddof=1) / np.sqrt(n_ev)
    row = dict(hits=n, w_min=w_min, events=n_ev,
               E6_counted=count[sets[0]].mean(), E6_counted_err=se(count[sets[0]]), E6_pred=pred[sets[0]],
               E5_counted=sum(count[s] for s in sets[1:]).mean(),
               E5_counted_err=se(sum(count[s] for s in sets[1:])),
               E5_pred=sum(pred[s] for s in sets[1:]),
               sum_counted=total.mean(), sum_counted_err=se(total),
               sum_pred=sum(pred.values()), sum_pred_err=np.sqrt(sum(v ** 2 for v in pred_err.values())),
               distinct4=distinct[4].mean(), distinct4_err=se(distinct[4]),
               distinct3=distinct[3].mean(), distinct3_err=se(distinct[3]),
               eta=float(np.mean(eta)) if eta else np.nan, seconds=time.time() - t0)
    print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in row.items()), flush=True)
    new = not os.path.exists("kofn_results.csv")
    with open("kofn_results.csv", "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        if new:
            w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    main()
