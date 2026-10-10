import sys, itertools, json, numpy as np
from helix_core import Tower
from helix_validation import p_direct, cut_for
w = float(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2]))
out = {}
for S in [tuple(range(6))] + list(itertools.combinations(range(6), 5)):
    t = Tower(6, w, planes=S)
    P, Pe, npass, ntr = p_direct(t, cut_for(t.ndof), 1_500_000_000, rng, want=1500)
    out[",".join(map(str, S))] = dict(P=P, P_err=Pe, n_pass=npass, n_trials=ntr)
    print(S, out[",".join(map(str, S))], flush=True)
json.dump(out, open(f"helix_pS_w{w:g}.json", "w"), indent=1)
