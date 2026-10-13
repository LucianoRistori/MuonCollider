"""
sig_z0.py - exact-helix fit of the real signal muons (their own IT+OT hits that
pass the hit cuts, at least 5), to get the distribution of the fitted track z0
and the windows that keep 98% of the signal, per region of the hit set:
  barrel        - barrel layers only
  intermediate  - barrel layers and disks
  endcap        - disks only
The model is the one used for the fakes: a helix from (0, 0, z0) with
theta = (phi0, kappa, z0, t = cot theta). Barrel hits are fitted on a
cylinder of the hit's own radius (residuals R dphi and z), disk hits on the
plane z = z_hit (residuals x and y); sigma = 0.1 mm for every coordinate.
"""
import json, sys
import numpy as np

UP = "/mnt/user-data/uploads/MuonColliderSimulation/Analysis/claude_scratch/regions"
SIG = 0.1
K_MIN = 5
KMAX10 = 0.2998 * 5.0 / 10.0 / 1000.0          # kappa at pT = 10 GeV/c


def pad_events(d):
    p = d["passed"]
    ev = d["event_id"][p]; o = np.argsort(ev, kind="stable")
    ev = ev[o]
    x, y, z = d["x"][p][o].astype(float), d["y"][p][o].astype(float), d["z"][p][o].astype(float)
    sysid, side, lay = d["system"][p][o], d["side"][p][o], d["layer"][p][o]
    ue, start, cnt = np.unique(ev, return_index=True, return_counts=True)
    keep = cnt >= K_MIN
    ue, start, cnt = ue[keep], start[keep], cnt[keep]
    N = cnt.max(); n = len(ue)
    idx = start[:, None] + np.arange(N)[None]
    valid = np.arange(N)[None] < cnt[:, None]
    idx = np.where(valid, idx, start[:, None])
    barrel = np.isin(sysid, [3, 5])
    A = dict(x=x[idx], y=y[idx], z=z[idx], barrel=barrel[idx], valid=valid,
             sys=sysid[idx], side=side[idx], lay=lay[idx])
    return ue, A


def model(th, A):
    """residuals (n, 2N) in units of sigma; invalid entries zero"""
    phi0, k, z0, t = (th[:, i:i + 1] for i in range(4))
    x, y, z, bar, val = A["x"], A["y"], A["z"], A["barrel"], A["valid"]
    R = np.hypot(x, y)
    # barrel: crossing of the cylinder of radius R
    a = np.clip(k * R / 2, -0.999999, 0.999999)
    sb = np.where(np.abs(k) > 1e-12, 2 * np.arcsin(a) / np.where(np.abs(k) > 1e-12, k, 1), R)
    phib = phi0 + np.arcsin(a)
    dphi = np.angle(np.exp(1j * (np.arctan2(y, x) - phib)))
    rb1 = R * dphi; rb2 = z - (z0 + t * sb)
    # disk: crossing of the plane z = z_hit
    sd = (z - z0) / np.where(np.abs(t) > 1e-9, t, 1e-9)
    h = k * sd / 2
    chord = sd * np.where(np.abs(h) > 1e-9, np.sin(h) / np.where(np.abs(h) > 1e-9, h, 1), 1)
    ang = phi0 + h
    rd1 = x - chord * np.cos(ang); rd2 = y - chord * np.sin(ang)
    r1 = np.where(bar, rb1, rd1); r2 = np.where(bar, rb2, rd2)
    bad = (~bar) & (sd < 0)                              # disk behind the origin
    r1 = np.where(bad, 1e3, r1)
    return np.concatenate([r1 * val, r2 * val], 1) / SIG


def fit(A, n_iter=30, k_seeds=(0.0, -1.0, 1.0)):
    x, y, z, val = A["x"], A["y"], A["z"], A["valid"]
    r = np.hypot(x, y); w = val.astype(float)
    # straight-line starting values in the r-z view and the inner-hit azimuth
    Sw = w.sum(1); Sr = (w * r).sum(1); Srr = (w * r * r).sum(1); Sz = (w * z).sum(1); Srz = (w * r * z).sum(1)
    det = Sw * Srr - Sr ** 2
    t0 = (Sw * Srz - Sr * Sz) / det; z00 = (Sz - t0 * Sr) / Sw
    inner = np.argmin(np.where(val, r, 1e9), 1)
    ph = np.arctan2(y[np.arange(len(x)), inner], x[np.arange(len(x)), inner])
    n = len(x)
    best = np.full(n, np.inf); bth = np.zeros((n, 4))
    eps = np.array([1e-7, 1e-9, 1e-5, 1e-7])
    for ks in k_seeds:
        th = np.stack([ph, np.full(n, ks * 1e-3), z00, t0], 1)
        Rr = model(th, A); c2 = (Rr ** 2).sum(1); lam = np.full(n, 1e-3)
        for _ in range(n_iter):
            J = np.stack([(model(th + eps[i] * np.eye(4)[i], A) - Rr) / eps[i] for i in range(4)], 2)
            JTJ = np.einsum("nia,nib->nab", J, J); g = -np.einsum("nia,ni->na", J, Rr)
            D = np.diagonal(JTJ, axis1=1, axis2=2)
            step = np.linalg.solve(JTJ + lam[:, None, None] * (np.eye(4)[None] * D[:, :, None]) + 1e-15 * np.eye(4)[None],
                                   g[..., None])[..., 0]
            tn = th + step; Rn = model(tn, A); c2n = (Rn ** 2).sum(1)
            b = c2n < c2
            th = np.where(b[:, None], tn, th); Rr = np.where(b[:, None], Rn, Rr); c2 = np.where(b, c2n, c2)
            lam = np.where(b, lam * 0.3, lam * 10)
        u = c2 < best; best = np.where(u, c2, best); bth = np.where(u[:, None], th, bth)
    return best, bth


def window98(v):
    """shortest interval containing 98% of v, and the symmetric |z0| <= a one"""
    v = np.sort(v); n = len(v); m = int(np.ceil(0.98 * n))
    w = v[m - 1:] - v[:n - m + 1]; i = np.argmin(w)
    return (float(v[i]), float(v[i + m - 1])), float(np.quantile(np.abs(v), 0.98))


def main():
    d = np.load(f"{UP}/sig_pass.npz")
    ev, A = pad_events(d)
    c2, th = fit(A)
    nh = A["valid"].sum(1)
    nb = (A["barrel"] & A["valid"]).sum(1)
    reg = np.where(nb == nh, "barrel", np.where(nb == 0, "endcap", "intermediate"))
    px, py, pz, vz, q = (d[f"gen_{k}"][ev] for k in ("px", "py", "pz", "vz", "q"))
    pt = np.hypot(px, py); tgen = pz / pt
    ndof = 2 * nh - 4
    out = {"n_tracks": int(len(ev)), "regions": {}}
    print(f"{len(ev)} muons with >= {K_MIN} passing IT+OT hits")
    print(f"chi2/ndof median {np.median(c2 / ndof):.2f}; fitted z0 - true vz: "
          f"rms {np.std(th[:, 2] - vz):.2f} mm, t - true t: rms {np.std(th[:, 3] - tgen):.4f}")
    for sel_name, sel in (("all pT", np.ones(len(ev), bool)), ("pT > 10", pt > 10)):
        print(f"--- {sel_name}")
        for R in ("barrel", "intermediate", "endcap"):
            m = sel & (reg == R)
            z0 = th[m, 2]
            (lo, hi), a = window98(z0)
            dz = z0 - vz[m]
            res = dict(n=int(m.sum()), z0_lo=lo, z0_hi=hi, z0_abs98=a,
                       resid_rms=float(np.std(dz)), resid_98abs=float(np.quantile(np.abs(dz), 0.98)),
                       t_range=[float(np.quantile(np.abs(th[m, 3]), q_)) for q_ in (0.01, 0.5, 0.99)])
            out["regions"].setdefault(sel_name, {})[R] = res
            print(f"{R:13s} n={m.sum():6d}  98% window [{lo:7.1f},{hi:7.1f}] mm  |z0|<{a:6.1f}  "
                  f"(z0-vz: rms {res['resid_rms']:.2f}, 98% {res['resid_98abs']:.1f})  |t| 1/50/99%: "
                  + " ".join(f"{v:.2f}" for v in res["t_range"]))
    np.savez_compressed("sig_fit.npz", ev=ev, th=th, c2=c2, nh=nh, nb=nb, pt=pt, vz=vz, tgen=tgen, q=q)
    json.dump(out, open("sig_z0_result.json", "w"), indent=1)


if __name__ == "__main__":
    main()
