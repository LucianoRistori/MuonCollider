import csv, json, time, numpy as np
from helix_core import Tower
from helix_validation import p_direct, cut_for
g = json.load(open("helix_kgeom.json"))["N4"]
rng = np.random.default_rng(4242); cut = cut_for(4); rows = []
for wm in (16, 24, 32, 48, 64, 96):
    cp = cut / (wm / 2) ** 2; Pg = g["K"] * cp ** 2; t0 = time.time()
    Pd, Pde, npass, ntr = p_direct(Tower(4, wm), cut, 3_000_000_000, rng, want=2000)
    rows.append(dict(N=4, ndof=4, w_min=wm, c_prime=cp, P_direct=Pd, P_direct_err=Pde, n_pass=npass, n_trials=ntr, P_geom=Pg, ratio_geom=Pd/Pg, ratio_geom_err=Pde/Pg))
    print(rows[-1], f"{time.time()-t0:.0f}s", flush=True)
    with open("helix_extrap_deep_N4.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
