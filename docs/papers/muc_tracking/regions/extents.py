import sys, os, numpy as np, csv
sys.path.insert(0, os.path.expanduser("~/mnt/MuonCollider"))
from bib_common import load_hits
h = load_hits(os.path.expanduser("~/mnt/MuonColliderSimulation/Data/ntu_bib_ipp_1bx.root"))
s, sd, l = np.asarray(h["system"]), np.asarray(h["side"]), np.asarray(h["layer"])
r = np.hypot(h["x"], h["y"]); z = np.asarray(h["z"])
w = csv.writer(open("claude_scratch/regions/extents.csv", "w")); w.writerow(["system", "side", "layer", "n", "rmin", "rmax", "r_p1", "r_p99", "zmin", "zmax", "absz_p1", "absz_p99"])
for key in sorted(set(zip(s.tolist(), sd.tolist(), l.tolist()))):
    m = (s == key[0]) & (sd == key[1]) & (l == key[2])
    rr, zz = r[m], z[m]
    w.writerow(list(key) + [int(m.sum()), rr.min(), rr.max(), np.percentile(rr, 1), np.percentile(rr, 99), zz.min(), zz.max(), np.percentile(abs(zz), 1), np.percentile(abs(zz), 99)])
