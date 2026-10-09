"""
validate_counting.py - the counting identity and "at least k of N planes",
checked by brute force in small detectors.

In each simulated noise-only event (n hits on each of 6 planes) EVERY
combination is fitted: one hit on each of the 6 planes (ndof 8), and one
hit on each plane of every 5-plane subset (ndof 6), each with its own
3-sigma cut. Counted per event:
  - E6, and E5(S) for each subset S: the passing combinations;
  - distinct fakes: passing candidates merged when they share >= m hits
    (m = 4 = k-1: a 6-hit fake and its passing 5-hit subsets, or two
    5-hit fakes that differ in one hit, count once; also m = 3).
These are compared with the prediction P_S x prod_{i in S} n, with P_S
measured directly for one random combination in the same geometry.

Usage: python3 validate_counting.py <w_min> <hits per plane> <events> [P trials]
Writes counting_w<w_min>_n<n>.csv (one row per set, plus sum and distinct).
"""
import csv
import itertools
import sys
import time

import numpy as np

from validate_line import cut_for, p_direct, projective


def main():
    w_min = float(sys.argv[1])
    n = int(sys.argv[2])
    n_ev = int(sys.argv[3])
    n_p = int(float(sys.argv[4])) if len(sys.argv) > 4 else 50_000_000
    N, k, sigma = 6, 5, 0.1
    rng = np.random.default_rng(int(100 * w_min) + n)
    geo = projective(N, 100, 100 * N, w_min * sigma * N, sigma)
    sets = [tuple(range(N))] + list(itertools.combinations(range(N), k))
    sub = {s: geo.subset(s) for s in sets}
    cut = {s: cut_for(sub[s].ndof) for s in sets}
    t0 = time.time()
    P = {s: p_direct(sub[s], cut[s], n_p, rng)[:2] for s in sets}
    print(f"w_min {w_min}, {n} hits/plane: direct P measured ({time.time()-t0:.0f}s)", flush=True)

    # all index combinations of each set, precomputed once
    idx = {s: np.stack([g.ravel() for g in np.meshgrid(*[np.arange(n)] * len(s), indexing="ij")], 1)
           for s in sets}
    count = {s: np.zeros(n_ev, int) for s in sets}
    distinct = {4: np.zeros(n_ev, int), 3: np.zeros(n_ev, int)}
    sub6_passing = 0              # 5-hit subsets of passing 6-hit fakes that also pass
    t0 = time.time()
    for e in range(n_ev):
        hu = np.stack([rng.uniform(-w / 2, w / 2, n) for w in geo.w])      # (N, n)
        hv = np.stack([rng.uniform(-w / 2, w / 2, n) for w in geo.w])
        cands = []
        for s in sets:
            ix = idx[s]
            u = np.stack([hu[i][ix[:, j]] for j, i in enumerate(s)], 1)
            v = np.stack([hv[i][ix[:, j]] for j, i in enumerate(s)], 1)
            ok = np.nonzero(sub[s].chi2(u, v) < cut[s])[0]
            count[s][e] = len(ok)
            for r in ok:
                cands.append(frozenset((i, int(ix[r, j])) for j, i in enumerate(s)))
        five = {c for c in cands if len(c) == 5}
        for c in cands:
            if len(c) == 6:
                sub6_passing += sum(1 for h in c if (c - {h}) in five)
        for m in (4, 3):
            parent = list(range(len(cands)))

            def find(a):
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a
            for a in range(len(cands)):
                for b in range(a + 1, len(cands)):
                    if len(cands[a] & cands[b]) >= m:
                        parent[find(a)] = find(b)
            distinct[m][e] = len({find(a) for a in range(len(cands))})
        if (e + 1) % max(1, n_ev // 10) == 0:
            print(f"  {e+1} events ({time.time()-t0:.0f}s)", flush=True)

    rows = []
    def add(label, counted, pred, pred_err):
        rows.append(dict(w_min=w_min, hits=n, events=n_ev, set=label, counted=counted.mean(),
                         counted_err=counted.std(ddof=1) / np.sqrt(n_ev), predicted=pred, predicted_err=pred_err))
    total = np.zeros(n_ev, int)
    pred_sum, pred_var = 0.0, 0.0
    for s in sets:
        p, pe = P[s]
        f = float(n) ** len(s)
        add("all 6" if len(s) == 6 else "5: no plane " + str(sorted(set(range(N)) - set(s))[0] + 1),
            count[s], p * f, pe * f)
        total += count[s]
        pred_sum += p * f
        pred_var += (pe * f) ** 2
    add("sum (6 of 6 + 5 of 6)", total, pred_sum, np.sqrt(pred_var))
    for m in (4, 3):
        add(f"distinct (merge if >= {m} shared hits)", distinct[m], np.nan, np.nan)
    print(f"\nw_min {w_min}, {n} hits/plane, {n_ev} events "
          f"(5-hit subsets of the 6-hit fakes that also pass: {sub6_passing / max(count[sets[0]].sum(), 1):.2f} per 6-hit fake)")
    for r in rows:
        print(f"  {r['set']:38s} counted {r['counted']:9.4g} +- {r['counted_err']:.2g}   "
              f"predicted {r['predicted']:9.4g} +- {r['predicted_err']:.2g}")
    with open(f"counting_w{w_min:g}_n{n}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
