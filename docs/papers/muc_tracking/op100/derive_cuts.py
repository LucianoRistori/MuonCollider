"""Derive the z0 and time cuts of every subsystem (98% signal containment, N-1),
as resolution_scan.py does, for the smearing in argv[1]; pT cut from the cuts config."""
import sys, os, numpy as np
sys.path.insert(0, os.path.expanduser("~/mnt/MuonCollider"))
import bib_common as bc
import resolution_scan as rs
smear = bc.load_smearing_config(sys.argv[1]); cuts_cfg = bc.load_cuts(sys.argv[2])
sig = bc.load_hits(os.path.expanduser("~/mnt/MuonColliderSimulation/Data/ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root"), muon_hits_only=True)
rng = bc.smearing_rng(smear)
bc.apply_position_time_smearing(sig, smear, rng); bc.add_incidence_angles(sig)
bc.apply_angle_smearing(sig, smear, rng); bc.add_time_of_flight(sig, B_FIELD_T=rs.B_FIELD_T)
pt_only = {n: dict(cuts_cfg[n]) for n in bc.CUT_NAMES}
for n in ("t_corrected_ns", "z_axis_intercept_mm"):
    pt_only[n]["enabled"] = False
sig_pt, _ = bc.apply_cuts(sig, pt_only, B_FIELD_T=rs.B_FIELD_T)
z0, tl, it = rs.n1_cuts(sig, sig_pt, 98.0)
for s in rs.SYSTEMS:
    print(f"{bc.SYSTEM_NAMES[s]:12s} time {tl[s]:.3f} ns   z0 {z0[s]:.0f} mm")
print("iterations", it)
