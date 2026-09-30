"""
Step 2 of rebuilding docs/efficiency_pt_infinity.pdf (see README.txt): every
number quoted in the note, from fitdata.npz, written to numbers.json.

The fit is the one of bib_common.efficiency_at_infinite_pt (least squares of
eps_inf + c/pT^2 to the individual trials above 10 GeV/c, sandwich errors,
hits grouped by muon), repeated here with the extra outputs the note shows.
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PT_MIN = 10.0
REGIONS = {**{str(k): [k] for k in range(1, 7)}, "VXD": [1, 2], "IT+OT": [3, 4, 5, 6],
           "ALL": [1, 2, 3, 4, 5, 6]}


def fit(pt, y, pmin=PT_MIN, groups=None):
    pt = np.asarray(pt, dtype=float)
    y = np.asarray(y, dtype=float)
    sel = pt > pmin
    n = int(sel.sum())
    yy = y[sel]
    X = np.column_stack([np.ones(n), (1.0 / pt[sel]) ** 2])
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ (X.T @ yy)
    r = yy - X @ beta
    sc = X * r[:, None]
    cov_ungrouped = xtx_inv @ (sc.T @ sc) @ xtx_inv
    cov = cov_ungrouped
    if groups is not None:
        g = np.unique(np.asarray(groups)[sel], return_inverse=True)[1]
        scg = np.column_stack([np.bincount(g, weights=X[:, j] * r) for j in range(2)])
        cov = xtx_inv @ (scg.T @ scg) @ xtx_inv
    x = X[:, 1]
    return dict(a=beta[0], c=beta[1], cov=cov.tolist(), sa=np.sqrt(cov[0, 0]), sc=np.sqrt(cov[1, 1]),
                sa_ungrouped=np.sqrt(cov_ungrouped[0, 0]), n=n, xmean=x.mean(), xrms=x.std(),
                avg_above=yy.mean(), avg_all=y.mean(),
                n_groups=(len(np.unique(np.asarray(groups)[sel])) if groups is not None else n))


def binned(pt, y, edges, groups=None):
    """Efficiency in bins of generated pT; errors per muon when grouped."""
    pt = np.asarray(pt, dtype=float)
    y = np.asarray(y, dtype=float)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (pt > lo) & (pt <= hi)
        n = int(s.sum())
        e = y[s].mean()
        if groups is None:
            err, ng = np.sqrt(e * (1 - e) / n), n
        else:
            g = np.unique(groups[s], return_inverse=True)[1]
            ng = int(g.max() + 1)
            err = np.sqrt(np.sum((np.bincount(g, weights=y[s]) - e * np.bincount(g)) ** 2)) / n
        x = 1.0 / pt[s] ** 2
        out.append(dict(lo=lo, hi=hi, n=n, ngroups=ng, eff=e, err=err, xmean=x.mean(),
                        xlo=1.0 / hi ** 2, xhi=1.0 / lo ** 2))
    return out


def main():
    d = np.load(HERE / "fitdata.npz")
    res = {"label": str(d["label"]), "date": str(d["date"])}
    tpt, found = d["trk_pt"], d["trk_found"]
    res["tracks"] = fit(tpt, found)
    res["tracks_thresholds"] = {str(p): fit(tpt, found, p) for p in (5, 7, 10, 15, 20, 40)}
    res["tracks_bins"] = binned(tpt, found, [1.5, 2, 2.5, 5, 10, 15, 20, 30, 50, 100, 1e9])
    q = d["trk_q"]
    res["tracks_frac"] = dict(above10=float((tpt > PT_MIN).mean()),
                              above30_of10=float((tpt > 30).sum() / (tpt > PT_MIN).sum()),
                              above100_of10=float((tpt > 100).sum() / (tpt > PT_MIN).sum()),
                              ptmax=float(tpt.max()), n=int(len(tpt)),
                              n_plus=int((q > 0).sum()), n_minus=int((q < 0).sum()))
    hpt, hp, hs, hev = d["hit_pt"], d["hit_passed"], d["hit_system"], d["hit_ev"]
    res["hits"] = {k: fit(hpt[np.isin(hs, ids)], hp[np.isin(hs, ids)], PT_MIN, hev[np.isin(hs, ids)])
                   for k, ids in REGIONS.items()}
    res["ot_barrel_turnon"] = binned(hpt[hs == 5], hp[hs == 5], [1.5, 2, 2.5, 3, 4, 10], hev[hs == 5])
    res["vxd_barrel_bins"] = binned(hpt[hs == 1], hp[hs == 1], [1.5, 2, 3, 5, 10, 20, 50, 1e9], hev[hs == 1])
    json.dump(res, open(HERE / "numbers.json", "w"), indent=1, default=float)

    t = res["tracks"]
    print(f"track-finding: eps_inf = {t['a']*100:.2f} +- {t['sa']*100:.2f} %, c = {t['c']:+.2f} +- "
          f"{t['sc']:.2f} (GeV/c)^2, {t['n']:,} muons")
    for k, f in res["hits"].items():
        print(f"hits {k:>5}: eps_inf = {f['a']*100:.2f} +- {f['sa']*100:.2f} %, c = {f['c']:+.2f} +- {f['sc']:.2f}")
    print(f"-> {HERE / 'numbers.json'}")


if __name__ == "__main__":
    main()
