"""IT and OT BIB hits (systems 3-6) with positions and the per-cut pass bits (time, z0, pT),
at the smearing and cuts of the given run, for fake-track estimates at intermediate stages."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.expanduser("~/mnt/MuonCollider"))
import bib_common as bc
run = sys.argv[1]; out = sys.argv[2]; t0 = time.time()
smear = bc.load_smearing_config(f"{run}/__smearing_config.txt"); rng = bc.smearing_rng(smear)
h = bc.load_hits(os.path.expanduser("~/mnt/MuonColliderSimulation/Data/ntu_bib_ipp_1bx.root"))
bc.apply_position_time_smearing(h, smear, rng); bc.add_incidence_angles(h); bc.apply_angle_smearing(h, smear, rng); bc.add_time_of_flight(h)
m_all, per = bc.apply_cuts(h, bc.load_cuts(f"{run}/__cuts_config.txt"))
sy = np.asarray(h["system"]); sel = sy >= 3
bits = (np.asarray(per["t_corrected_ns"]).astype(np.uint8) | (np.asarray(per["z_axis_intercept_mm"]).astype(np.uint8) << 1)
        | (np.asarray(per["momentum_gev"]).astype(np.uint8) << 2))
np.savez_compressed(out, system=sy[sel].astype(np.int8), side=np.asarray(h["side"])[sel].astype(np.int8),
                    layer=np.asarray(h["layer"])[sel].astype(np.int8),
                    x=np.asarray(h["x"])[sel].astype(np.float32), y=np.asarray(h["y"])[sel].astype(np.float32),
                    z=np.asarray(h["z"])[sel].astype(np.float32), bits=bits[sel])
print("IT+OT hits", int(sel.sum()), "all pass", int((m_all & sel).sum()), f"{time.time()-t0:.0f}s")
