"""Calibrate K for the 6-plane set and its six 5-plane subsets, on the
reference projective geometry (planes at 100..600 mm, sides prop. to Y,
smallest whitened side 2). Writes calib_sets.json."""
import itertools, json, sys, time
import numpy as np
from validate_line import calibrate_K, cut_for, projective
n_cal = int(float(sys.argv[1])) if len(sys.argv) > 1 else 200_000_000
rng = np.random.default_rng(424242)
ref = projective(6, 100, 600, 2 * 0.1 * 6, 0.1)
out = {}
for s in [tuple(range(6))] + list(itertools.combinations(range(6), 5)):
    g = ref.subset(s)
    t0 = time.time()
    cal = calibrate_K(g, 1.0, n_cal, rng)
    out[",".join(map(str, s))] = dict(ndof=g.ndof, cut=cut_for(g.ndof), K=cal["K"], K_err=cal["K_err"],
                                       beta=cal["beta"], c_max=cal["c_max"], n_fit=cal["n_fit"])
    print(s, out[",".join(map(str, s))], f"{time.time()-t0:.0f}s", flush=True)
json.dump(out, open("calib_sets.json", "w"), indent=1)
