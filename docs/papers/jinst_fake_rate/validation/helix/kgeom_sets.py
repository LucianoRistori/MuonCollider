import itertools, json
from helix_core import Tower
from kgeom import K_geom
from helix_validation import cut_for
out = {}
t = Tower(4, 2.0); out["N4"] = dict(K=K_geom(t, 400, 400, 600), ndof=4, cut=cut_for(4))
print(out["N4"], flush=True)
for S in [tuple(range(6))] + list(itertools.combinations(range(6), 5)):
    t = Tower(6, 2.0, planes=S)
    out[",".join(map(str, S))] = dict(K=K_geom(t, 400, 400, 600), ndof=t.ndof, cut=cut_for(t.ndof))
    print(S, out[",".join(map(str, S))], flush=True)
json.dump(out, open("helix_kgeom.json", "w"), indent=1)
