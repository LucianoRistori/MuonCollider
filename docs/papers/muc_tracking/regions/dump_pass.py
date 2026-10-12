import sys, numpy as np, time
sys.path.insert(0, "/sessions/" ) if False else None
import os
sys.path.insert(0, os.path.expanduser("~/mnt/MuonCollider"))
from bib_common import (load_hits, add_incidence_angles, add_time_of_flight, load_cuts, apply_cuts,
                        load_smearing_config, smearing_rng, apply_position_time_smearing, apply_angle_smearing)
run = sys.argv[1]; out = sys.argv[2]
t0 = time.time()
smear = load_smearing_config(f"{run}/__smearing_config.txt"); rng = smearing_rng(smear)
h = load_hits(os.path.expanduser("~/mnt/MuonColliderSimulation/Data/ntu_bib_ipp_1bx.root"))
apply_position_time_smearing(h, smear, rng); add_incidence_angles(h); apply_angle_smearing(h, smear, rng); add_time_of_flight(h)
cuts = load_cuts(f"{run}/__cuts_config.txt")
m, _ = apply_cuts(h, cuts)
keys = [k for k in ("system", "side", "layer", "x", "y", "z") if k in h]
np.savez_compressed(out, pass_mask_count=int(m.sum()), **{k: np.asarray(h[k])[m] for k in keys},
                    **{f"all_{k}": np.asarray(h[k]) for k in ("system", "side", "layer")})
print("passing", int(m.sum()), "of", len(m), f"{time.time()-t0:.0f}s", keys)
