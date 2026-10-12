import sys, time, json, numpy as np, gsurf as g
from scipy.stats import chi2
from geoms import GEOMS
name = sys.argv[1]; Ls = [float(x) for x in sys.argv[2].split(",")]; n_max = int(float(sys.argv[3]))
cut = chi2.isf(1.3499e-3, 4); rng = np.random.default_rng(abs(hash(name)) % 2**31)
out = []
for L in Ls:
    surfs, kmax, box = GEOMS[name](L)
    t0 = time.time()
    Pt, Pt_err, edge = g.tube_P_fast(surfs, box, kmax, cut, 20_000_000, rng)
    t1 = time.time()
    npass = done = lost = 0
    while done < n_max and (npass < 400 or done < 200_000):
        m = 20_000
        UV = [S.random_hits(rng, m) for S in surfs]
        U = np.stack([a for a, b in UV], 1); V = np.stack([b for a, b in UV], 1)
        c2q, thq = g.fit_fast(surfs, U, V, kmax, n_iter=6, k_seeds=(0.0,))
        sel = np.nonzero(c2q < 20 * cut)[0]
        c2 = c2q.copy(); th = thq.copy()
        if len(sel):
            c2s, ths = g.fit_fast(surfs, U[sel], V[sel], kmax)
            b = c2s < c2[sel]; c2[sel[b]] = c2s[b]; th[sel[b]] = ths[b]
        npass += int(np.sum((c2 < cut) & (np.abs(th[:, 1]) <= kmax)))
        done += m
    Pd = npass / done; Pd_err = np.sqrt(npass) / done
    r = dict(geom=name, L=L, P_tube=Pt, P_tube_err=Pt_err, edge=edge, P_direct=Pd, P_direct_err=Pd_err,
             n_pass=npass, n_trials=done, ratio=Pd / Pt, ratio_err=Pd_err / Pt, t_tube=t1 - t0, t_direct=time.time() - t1)
    out.append(r)
    print(json.dumps(r), flush=True)
    json.dump(out, open(f"val_{name}_{sys.argv[2].replace(",","_")}.json", "w"), indent=1)
