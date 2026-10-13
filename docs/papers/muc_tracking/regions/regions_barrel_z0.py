"""barrel-only fakes (sets of barrel layers only) with the track z0 restricted to |z0| <= Z0W,
the window that keeps 98% of the fitted real barrel tracks (sig_z0.py)."""
import sys, math, json, numpy as np
import regions as R
Z0W = float(sys.argv[1]); n_samp = int(float(sys.argv[2])); seed = int(sys.argv[3])
surfs = [S for S in R.load_geometry() if S["kind"] == "barrel"]
rng = np.random.default_rng(seed)
DPHI = 0.2
box = [(0.0, DPHI), (-R.KMAX, R.KMAX), (-Z0W, Z0W), (0.0, 3.0)]
lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
vol = np.prod(hi - lo) * (2 * np.pi / DPHI) * 2
s1 = s2 = 0.0; accm = {}; done = 0; B = 100_000
tb = np.linspace(0, 3, 13); ht = np.zeros(12)
while done < n_samp:
    th = lo + (hi - lo) * rng.random((B, 4))
    res, by_m, _ = R.weights(surfs, th)
    w = res["barrel"]; s1 += w.sum(); s2 += (w ** 2).sum()
    ht += np.histogram(th[:, 3], tb, weights=w)[0]
    for m, x in by_m.items(): accm[m] = accm.get(m, 0) + x.sum()
    done += B
mean = s1 / done; err = math.sqrt(max(s2 / done - mean ** 2, 0) / done)
print(f"|z0|<={Z0W} mm: barrel E = {vol*mean:.3e} +- {vol*err:.1e} per BX;",
      {m: f"{vol*v/done:.2e}" for m, v in accm.items()})
print("t hist:", " ".join(f"{vol*x/done:.1e}" for x in ht))
json.dump(dict(z0w=Z0W, E=vol*mean, E_err=vol*err, by_m={m: vol*v/done for m, v in accm.items()},
               t_bins=tb.tolist(), t_hist=(vol*ht/done).tolist()), open(f"barrel_z0_{Z0W}_s{seed}.json", "w"))
