"""
regions.py - expected fake tracks per bunch crossing in the IT+OT tracker of
the MuSIC_v2 geometry, split into barrel / intermediate / endcap, with the
general-surface formula:

  E = Int dtheta  sum over sets S (>= K of the surfaces crossed by theta)
        V_d c_S^(d/2) sqrt(det J_S^T J_S) prod_{i in S} rho_i(x_i) sigma^2

Surfaces: IT and OT barrel layers as full coaxial cylinders (|z| <= half
length), IT and OT endcap disks as planes z = +-z_d with annular active
regions. Densities after the cuts: rho(|z|) per barrel layer, rho(r) per disk
(both ends averaged), from the BIB hits passing the cuts of a run.
Track parameters sampled uniformly: phi0 in a narrow window (phi symmetry),
kappa within the p_T threshold, z0 free, t = cot(theta) > 0 (the t < 0 half is
the mirror image). A fake on set S is "barrel" if S has only barrel layers,
"endcap" if only disks, "intermediate" if both.
"""
import csv, itertools, json, math, os, sys
import numpy as np
from scipy.stats import chi2 as chi2dist
import gsurf as g

P3 = 1.3499e-3
SIG = 0.1
K_MIN = 5
KMAX = 0.2998 * 5.0 / 10.0 / 1000.0
UP = "/mnt/user-data/uploads/MuonColliderSimulation/Analysis"


def load_geometry():
    ext = {}
    for r in csv.DictReader(open(f"{UP}/claude_scratch/regions/extents.csv")):
        ext[(int(r["system"]), int(r["side"]), int(r["layer"]))] = {k: float(v) for k, v in r.items()}
    area = {}
    for r in csv.DictReader(open(os.environ.get("REG_DENS", f"{UP}/runs/2026-10-10_112338_pdf/step4_cuts/density_per_layer_before_after_cuts.csv"))):
        area[(int(r["system"]), int(r["layer"]))] = float(r["area_mm2"])
    hits = np.load(os.environ.get("REG_BIB", f"{UP}/claude_scratch/regions/bib_pass.npz"))
    s, sd, l = hits["system"], hits["side"], hits["layer"]
    r = np.hypot(hits["x"], hits["y"]); z = hits["z"]
    surfs = []
    for sysid, nl in ((3, 3), (5, 3)):                       # barrel cylinders
        for L in range(nl):
            e = ext[(sysid, 0, L)]
            R = 0.5 * (e["rmin"] + e["rmax"]); Z = e["zmax"]
            cov = area[(sysid, L)] / (2 * np.pi * R * 2 * Z)
            edges = np.linspace(0, Z, 13)
            m = (s == sysid) & (l == L)
            cnt, _ = np.histogram(np.abs(z[m]), edges)
            rho = cnt / (cov * 2 * np.pi * R * 2 * np.diff(edges))
            surfs.append(dict(kind="barrel", name=f"{'IT' if sysid == 3 else 'OT'}B{L}", R=R, Z=Z,
                              edges=edges, rho=rho, n_after=int(m.sum()), cov=cov,
                              obj=g.Cylinder(R, 0.0, np.pi, -Z, Z, SIG)))
    for sysid, nl in ((4, 7), (6, 4)):                       # endcap disks (+z side; -z mirrored)
        for L in range(nl):
            e = ext[(sysid, 1, L)]
            zd = 0.5 * (e["zmin"] + e["zmax"]); r0, r1 = e["rmin"], e["rmax"]
            cov = area[(sysid, L)] / (2 * np.pi * (r1 ** 2 - r0 ** 2))
            edges = np.linspace(r0, r1, 9)
            m = (s == sysid) & (l == L)
            cnt, _ = np.histogram(r[m], edges)
            rho = cnt / (cov * 2 * np.pi * np.diff(edges ** 2))      # both ends together
            surfs.append(dict(kind="disk", name=f"{'IT' if sysid == 4 else 'OT'}D{L}", zd=zd, r0=r0, r1=r1,
                              edges=edges, rho=rho, n_after=int(m.sum()), cov=cov,
                              obj=g.Plane([0, 0, zd], [0, 0, 1], [1, 0, 0], -r1, r1, -r1, r1, SIG)))
    return surfs


def crossing_info(surfs, th):
    """for each surface: inside mask, J (n,2,4) whitened, density at crossing"""
    out = []
    for S in surfs:
        u, v, ok, _, s, du, dv = g.cross_jac(S["obj"], th)
        if S["kind"] == "barrel":
            inside = ok & (np.abs(v) <= S["Z"])
            x = np.abs(v)
        else:
            rr = np.hypot(u, v)
            inside = ok & (rr >= S["r0"]) & (rr <= S["r1"])
            x = rr
        idx = np.clip(np.searchsorted(S["edges"], x) - 1, 0, len(S["rho"]) - 1)
        rho = S["rho"][idx]
        J = np.stack([du, dv], 1) / SIG
        out.append((inside, J, rho))
    return out


def weights(surfs, th):
    """per-sample fake density, split by region of the layer set, and by number of layers"""
    n = len(th); M = len(surfs)
    info = crossing_info(surfs, th)
    crossed = np.stack([i[0] for i in info], 1)
    kinds = np.array([S["kind"] == "barrel" for S in surfs])
    res = {k: np.zeros(n) for k in ("barrel", "intermediate", "endcap")}
    by_m = {}
    pats, inv = np.unique(crossed, axis=0, return_inverse=True)
    for p_i, pat in enumerate(pats):
        idx_s = np.nonzero(pat)[0]
        if len(idx_s) < K_MIN:
            continue
        rows = np.nonzero(inv.ravel() == p_i)[0]
        Jall = {i: info[i][1][rows] for i in idx_s}
        rhoall = {i: info[i][2][rows] for i in idx_s}
        for k in range(K_MIN, len(idx_s) + 1):
            d = 2 * k - 4; c = chi2dist.isf(P3, d)
            pref = math.pi ** (d / 2) / math.gamma(d / 2 + 1) * c ** (d / 2)
            for S in itertools.combinations(idx_s, k):
                J = np.concatenate([Jall[i] for i in S], 1)
                det = np.linalg.det(np.einsum("nia,nib->nab", J, J))
                w = pref * np.sqrt(np.maximum(det, 0)) * np.prod([rhoall[i] * SIG ** 2 for i in S], 0)
                nb = kinds[list(S)].sum()
                reg = "barrel" if nb == k else ("endcap" if nb == 0 else "intermediate")
                res[reg][rows] += w
                by_m.setdefault(k, np.zeros(n))[rows] += w
    return res, by_m, crossed


def main():
    n_samp = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2_000_000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    surfs = load_geometry()
    rng = np.random.default_rng(seed)
    DPHI = 0.2
    box = [(0.0, DPHI), (-KMAX, KMAX), (-1500.0, 1500.0), (0.0, 8.0)]
    lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
    vol = np.prod(hi - lo) * (2 * np.pi / DPHI) * 2          # full phi, and both signs of t
    acc = {k: [0.0, 0.0] for k in ("barrel", "intermediate", "endcap")}
    accm = {}
    z0_bins = np.linspace(-1500, 1500, 61); t_bins = np.linspace(0, 8, 33)
    hz = {k: np.zeros(60) for k in acc}; ht = {k: np.zeros(32) for k in acc}
    done = 0; B = 100_000
    while done < n_samp:
        th = lo + (hi - lo) * rng.random((B, 4))
        res, by_m, _ = weights(surfs, th)
        for k, w in res.items():
            acc[k][0] += w.sum(); acc[k][1] += (w ** 2).sum()
            hz[k] += np.histogram(th[:, 2], z0_bins, weights=w)[0]
            ht[k] += np.histogram(th[:, 3], t_bins, weights=w)[0]
        for m, w in by_m.items():
            accm[m] = accm.get(m, 0.0) + w.sum()
        done += B
    out = {"n_samples": done, "volume": vol, "surfaces": [{k: (v.tolist() if isinstance(v, np.ndarray) else v)
           for k, v in S.items() if k != "obj"} for S in surfs]}
    for k, (s1, s2) in acc.items():
        mean = s1 / done; err = math.sqrt(max(s2 / done - mean ** 2, 0) / done)
        out[k] = dict(E=vol * mean, E_err=vol * err,
                      z0_hist=(vol * hz[k] / done).tolist(), t_hist=(vol * ht[k] / done).tolist())
        print(f"{k:13s} E = {vol*mean:.3e} +- {vol*err:.1e} per bunch crossing", flush=True)
    out["by_n_layers"] = {m: vol * v / done for m, v in accm.items()}
    out["z0_bins"] = z0_bins.tolist(); out["t_bins"] = t_bins.tolist()
    print("by number of layers:", {m: f"{v:.2e}" for m, v in out["by_n_layers"].items()})
    json.dump(out, open(f"regions_result_s{seed}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
