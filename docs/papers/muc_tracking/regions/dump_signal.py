"""dump_signal.py <run folder> <out.npz>: the muon-gun hits (the muon's own),
smeared and cut exactly as in the track-efficiency step of that run, IT+OT
systems only, with the pass mask and the generated muon's parameters."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.expanduser("~/mnt/MuonCollider"))
from bib_common import (load_hits, add_incidence_angles, add_time_of_flight, load_cuts, apply_cuts,
                        load_smearing_config, smearing_rng, apply_position_time_smearing,
                        apply_angle_smearing, signal_muon_hits_only)
run, out = sys.argv[1], sys.argv[2]; t0 = time.time()
cfg = f"{run}/__cuts_config.txt"
smear = load_smearing_config(f"{run}/__smearing_config.txt"); rng = smearing_rng(smear)
h = load_hits(os.path.expanduser("~/mnt/MuonColliderSimulation/Data/ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root"),
              muon_hits_only=signal_muon_hits_only(cfg))
apply_position_time_smearing(h, smear, rng); add_incidence_angles(h); apply_angle_smearing(h, smear, rng); add_time_of_flight(h)
m, _ = apply_cuts(h, load_cuts(cfg))
sel = np.isin(h["system"], [3, 4, 5, 6])
keys = ("event_id", "system", "side", "layer", "x", "y", "z")
P = h["_primary"]
np.savez_compressed(out, passed=np.asarray(m)[sel], **{k: np.asarray(h[k])[sel] for k in keys},
                    **{f"gen_{k}": np.asarray(P[k], float) for k in ("px", "py", "pz", "vx", "vy", "vz", "q")},
                    n_events=h["_n_events"])
print("IT+OT hits", int(sel.sum()), "passing", int(m[sel].sum()), f"{time.time()-t0:.0f}s")
