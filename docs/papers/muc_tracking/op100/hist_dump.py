"""Histograms for paper 2 figures (IT and OT barrel): corrected time, z-axis intercept and
curvature pT, for BIB and signal, all hits and N-1 (passing the other two cuts), at the
smearing and cuts of the given run."""
import sys, os, numpy as np
sys.path.insert(0, os.path.expanduser("~/mnt/MuonCollider"))
import bib_common as bc
run = sys.argv[1]; out = sys.argv[2]
smear = bc.load_smearing_config(f"{run}/__smearing_config.txt"); cuts = bc.load_cuts(f"{run}/__cuts_config.txt")
D = os.path.expanduser("~/mnt/MuonColliderSimulation/Data")
res = {}
for name, f, mu in (("bib", "ntu_bib_ipp_1bx.root", False), ("sig", "ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root", True)):
    rng = bc.smearing_rng(smear)
    h = bc.load_hits(f"{D}/{f}", muon_hits_only=mu)
    bc.apply_position_time_smearing(h, smear, rng); bc.add_incidence_angles(h); bc.apply_angle_smearing(h, smear, rng)
    bc.add_time_of_flight(h)
    _, per = bc.apply_cuts(h, cuts)
    keys = list(per.keys())
    sy = np.asarray(h["system"])
    t = np.asarray(h["t_corrected_ns"]); z = np.asarray(h["z_axis_intercept_mm"]); k = np.asarray(h["inv_radius_per_mm"]) * 1000.0
    for s, lab in ((3, "IT"), (5, "OT"), (4, "ITE"), (6, "OTE")):
        m = sy == s
        pt, pz, pm = (np.asarray(per[c])[m] for c in ("t_corrected_ns", "z_axis_intercept_mm", "momentum_gev"))
        for var, x, edges, other in (("t", t[m], np.linspace(-1, 3, 161), pz & pm), ("z0", z[m], np.linspace(-1000, 1000, 201) if s in (4, 6) else np.linspace(-400, 400, 161), pt & pm),
                                     ("k", k[m], np.linspace(-3, 3, 121), pt & pz)):
            res[f"{name}_{lab}_{var}_all"] = np.histogram(x, edges)[0]
            res[f"{name}_{lab}_{var}_n1"] = np.histogram(x, edges)[0] if False else np.histogram(x[other], edges)[0]
            res[f"edges_{var}"] = edges; res[f"edges_{lab}_{var}"] = edges
        res[f"{name}_{lab}_n"] = int(m.sum()); res[f"{name}_{lab}_pass"] = int((pt & pz & pm).sum())
    print(name, keys)
np.savez_compressed(out, **res)
