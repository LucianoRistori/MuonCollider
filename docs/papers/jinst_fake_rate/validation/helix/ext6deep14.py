import csv, json, time, numpy as np
from helix_core import Tower
from helix_validation import p_direct, cut_for
g = json.load(open("helix_kgeom.json"))["0,1,2,3,4,5"]
rng = np.random.default_rng(6262); cut = cut_for(8); rows = []
for wm in (14,):
    cp = cut / (wm / 2) ** 2; Pg = g["K"] * cp ** 4; t0 = time.time()
    Pd, Pde, npass, ntr = p_direct(Tower(6, wm), cut, 6_000_000_000, rng, want=300)
    rows.append(dict(N=6, ndof=8, w_min=wm, c_prime=cp, P_direct=Pd, P_direct_err=Pde, n_pass=npass, n_trials=ntr, P_geom=Pg, ratio_geom=Pd/Pg, ratio_geom_err=Pde/Pg))
    print(rows[-1], f"{time.time()-t0:.0f}s", flush=True)
    with open("helix_extrap_deep_N6_w14.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
