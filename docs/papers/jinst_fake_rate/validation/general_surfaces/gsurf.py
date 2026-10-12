"""
gsurf.py - the fake-probability formula for arbitrary detector surfaces.

Track model: helix from the beam line (transverse circle through the origin),
parameters theta = (phi0, kappa, z0, t):  phi0 initial azimuth, kappa signed
curvature, z0 longitudinal origin, t = cot(polar angle). At transverse arc
length s:
    x = s cos(phi0 + kappa s/2) sinc(kappa s/2)
    y = s sin(phi0 + kappa s/2) sinc(kappa s/2)
    z = z0 + t s
Surfaces: Cylinder (coaxial, local coords u = R*phi, v = z) and Plane (any
position and orientation: local coords along two in-plane unit vectors).
Each surface has an active region (a rectangle in its local coordinates),
resolution sigma in both local coordinates.

Leading-order formula (tube around the track manifold):
    P = V_d c^(d/2) * prod_i sigma_i^2 / A_i * Int_{Theta_B} sqrt(det J^T J) dtheta
with J the derivatives of the whitened local crossing coordinates.
"""
import math
import numpy as np


def sinc(u):
    out = np.ones_like(u)
    m = np.abs(u) > 1e-6
    out[m] = np.sin(u[m]) / u[m]
    out[~m] = 1 - u[~m] ** 2 / 6
    return out


def helix_point(th, s):
    """th: (n,4) phi0, kappa, z0, t; s: (n,) -> (n,3)"""
    phi0, k, z0, t = th.T
    h = k * s / 2
    a = phi0 + h
    sc = s * sinc(h)
    return np.stack([sc * np.cos(a), sc * np.sin(a), z0 + t * s], 1)


def helix_tangent(th, s):
    phi0, k, z0, t = th.T
    a = phi0 + k * s
    return np.stack([np.cos(a), np.sin(a), t], 1)


def asin_over(x):
    """arcsin(x)/x, stable at 0"""
    out = np.ones_like(x)
    m = np.abs(x) > 1e-6
    out[m] = np.arcsin(np.clip(x[m], -1, 1)) / x[m]
    out[~m] = 1 + x[~m] ** 2 / 6
    return out


class Cylinder:
    def __init__(self, R, phi_c, half_dphi, z_lo, z_hi, sigma):
        self.R, self.phi_c, self.half_dphi = R, phi_c, half_dphi
        self.z_lo, self.z_hi, self.sigma = z_lo, z_hi, sigma
        self.area = 2 * half_dphi * R * (z_hi - z_lo)

    def cross(self, th):
        k = th[:, 1]
        x = k * self.R / 2
        ok = np.abs(x) < 1
        s = self.R * asin_over(np.where(ok, x, 0.0))      # first crossing
        p = helix_point(th, s)
        phi = np.arctan2(p[:, 1], p[:, 0])
        dphi = (phi - self.phi_c + np.pi) % (2 * np.pi) - np.pi
        u, v = self.R * dphi, p[:, 2]
        inside = ok & (np.abs(dphi) <= self.half_dphi) & (v >= self.z_lo) & (v <= self.z_hi)
        return u, v, ok, inside

    def random_hits(self, rng, n):
        u = rng.uniform(-self.half_dphi, self.half_dphi, n) * self.R
        v = rng.uniform(self.z_lo, self.z_hi, n)
        return u, v


class Plane:
    """Plane through P0 with unit normal n; local axes eu, ev (orthonormal, in plane);
    active rectangle u in [u_lo, u_hi], v in [v_lo, v_hi]."""
    def __init__(self, P0, n, eu, u_lo, u_hi, v_lo, v_hi, sigma):
        self.P0 = np.asarray(P0, float)
        self.n = np.asarray(n, float) / np.linalg.norm(n)
        eu = np.asarray(eu, float); eu = eu - (eu @ self.n) * self.n
        self.eu = eu / np.linalg.norm(eu)
        self.ev = np.cross(self.n, self.eu)
        self.u_lo, self.u_hi, self.v_lo, self.v_hi, self.sigma = u_lo, u_hi, v_lo, v_hi, sigma
        self.area = (u_hi - u_lo) * (v_hi - v_lo)

    def cross(self, th, n_iter=8):
        phi0, k, z0, t = th.T
        d = np.stack([np.cos(phi0), np.sin(phi0), t], 1)
        o = np.stack([np.zeros_like(z0), np.zeros_like(z0), z0], 1)
        den = d @ self.n
        s = np.where(np.abs(den) > 1e-12, ((self.P0 - o) @ self.n) / np.where(np.abs(den) > 1e-12, den, 1), -1.0)
        for _ in range(n_iter):
            f = (helix_point(th, s) - self.P0) @ self.n
            fp = helix_tangent(th, s) @ self.n
            fp = np.where(np.abs(fp) > 1e-12, fp, 1e-12)
            step = f / fp
            s = s - np.clip(step, -abs(1.0) * (1 + np.abs(s)), (1 + np.abs(s)))
        p = helix_point(th, s)
        res = np.abs((p - self.P0) @ self.n)
        ok = (s > 0) & (res < 1e-9 * (1 + np.abs(s))) & (np.abs(k * s) < np.pi)   # before the half turn
        q = p - self.P0
        u, v = q @ self.eu, q @ self.ev
        inside = ok & (u >= self.u_lo) & (u <= self.u_hi) & (v >= self.v_lo) & (v <= self.v_hi)
        return u, v, ok, inside

    def random_hits(self, rng, n):
        return rng.uniform(self.u_lo, self.u_hi, n), rng.uniform(self.v_lo, self.v_hi, n)


def crossings(surfs, th):
    us, vs, oks, ins = zip(*(S.cross(th) for S in surfs))
    return np.stack(us, 1), np.stack(vs, 1), np.all(np.stack(oks, 1), 1), np.all(np.stack(ins, 1), 1)


def whitened_jacobian(surfs, th, h=(1e-6, 1e-7, 1e-6, 1e-6)):
    """(n, 2N, 4) derivatives of whitened local coordinates"""
    n, N = len(th), len(surfs)
    sig = np.array([S.sigma for S in surfs])
    J = np.zeros((n, 2 * N, 4))
    for a in range(4):
        dp = th.copy(); dm = th.copy()
        dp[:, a] += h[a]; dm[:, a] -= h[a]
        up, vp, _, _ = crossings(surfs, dp)
        um, vm, _, _ = crossings(surfs, dm)
        J[:, :N, a] = (up - um) / (2 * h[a]) / sig
        J[:, N:, a] = (vp - vm) / (2 * h[a]) / sig
    return J


def tube_P(surfs, box, kmax, cut, n_samples, rng, batch=200_000):
    """P = V_d c^(d/2) prod sigma^2/A * Int sqrt(det J^T J) over theta in box with
    |kappa| <= kmax and all crossings inside. box: [(lo, hi)] for phi0, kappa, z0, t.
    Returns P, its MC (integration) error, and the fraction of accepted samples on the box edges."""
    N = len(surfs); d = 2 * N - 4
    lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
    vol_box = np.prod(hi - lo)
    tot = tot2 = 0.0; edge = 0; acc = 0; done = 0
    while done < n_samples:
        m = min(batch, n_samples - done)
        th = lo + (hi - lo) * rng.random((m, 4))
        _, _, ok, inside = crossings(surfs, th)
        sel = inside & (np.abs(th[:, 1]) <= kmax)
        g = np.zeros(m)
        if sel.any():
            J = whitened_jacobian(surfs, th[sel])
            g[sel] = np.sqrt(np.maximum(np.linalg.det(np.einsum("nia,nib->nab", J, J)), 0))
            e = ((th[sel] - lo) / (hi - lo))[:, [0, 2, 3]]
            edge += int(np.sum(np.any((e < 0.002) | (e > 0.998), axis=1)))
            acc += int(sel.sum())
        tot += g.sum(); tot2 += (g ** 2).sum(); done += m
    mean = tot / done
    err = math.sqrt(max(tot2 / done - mean ** 2, 0) / done)
    Vd = math.pi ** (d / 2) / math.gamma(d / 2 + 1)
    pref = Vd * cut ** (d / 2) * np.prod([S.sigma ** 2 / S.area for S in surfs]) * vol_box
    return pref * mean, pref * err, edge / max(acc, 1)


# ------------------------------------------------------------------ the fit
def fit_two_stage(surfs, U, V, kmax, cut, factor=20.0):
    """Quick single-start fit for all combinations; the full multi-start fit
    only where the quick chi2 is below factor*cut. Returns chi2, theta, and
    the quick chi2 (for checking that the prefilter loses nothing)."""
    c2q, thq = fit(surfs, U, V, kmax, n_iter=6, k_seeds=(0.0,))
    c2 = c2q.copy(); th = thq.copy()
    sel = np.nonzero(c2q < factor * cut)[0]
    if len(sel):
        c2[sel], th[sel] = fit(surfs, U[sel], V[sel], kmax)
        c2[sel] = np.minimum(c2[sel], c2q[sel])
    return c2, th, c2q


def fit(surfs, U, V, kmax, n_iter=20, k_seeds=(0.0, -0.6, 0.6)):
    """Least-squares helix fit to hits (U, V: (n, N) local coordinates).
    Multi-start Levenberg-Marquardt. Returns chi2 (n,), theta (n,4)."""
    n, N = U.shape
    sig = np.array([S.sigma for S in surfs])
    # starting values: global positions of the hits, a straight line from the beam line
    P = np.stack([local_to_global(S, U[:, i], V[:, i]) for i, S in enumerate(surfs)], 1)   # (n,N,3)
    r = np.hypot(P[..., 0], P[..., 1])
    phi = np.arctan2(P[..., 1], P[..., 0])
    phi0 = np.angle(np.mean(np.exp(1j * phi), 1))
    Sr, Srr, Sz, Srz = r.sum(1), (r * r).sum(1), P[..., 2].sum(1), (r * P[..., 2]).sum(1)
    det = N * Srr - Sr ** 2
    t0 = (N * Srz - Sr * Sz) / det
    z00 = (Sz - t0 * Sr) / N
    best_c2 = np.full(n, np.inf); best_th = np.zeros((n, 4))

    def resid(th):
        u, v, ok, _ = crossings(surfs, th)
        R = np.concatenate([(U - u) / sig, (V - v) / sig], 1)
        R[~ok] = 1e6
        return R

    for ks in k_seeds:
        th = np.stack([phi0, np.full(n, ks * kmax), z00, t0], 1)
        lam = np.full(n, 1e-3)
        R = resid(th); c2 = (R ** 2).sum(1)
        for _ in range(n_iter):
            J = -whitened_jacobian(surfs, th) * 1.0
            # residual derivative = -(d model / d theta)/sigma -> J as above (sign absorbed)
            JTJ = np.einsum("nia,nib->nab", J, J)
            g = np.einsum("nia,ni->na", J, R)
            Ad = JTJ + lam[:, None, None] * (np.eye(4)[None] * np.diagonal(JTJ, axis1=1, axis2=2)[:, :, None] + 1e-12 * np.eye(4)[None])
            try:
                step = -np.linalg.solve(Ad, g[..., None])[..., 0]
            except np.linalg.LinAlgError:
                break
            tn = th + step
            Rn = resid(tn); c2n = (Rn ** 2).sum(1)
            better = c2n < c2
            th = np.where(better[:, None], tn, th)
            R = np.where(better[:, None], Rn, R)
            c2 = np.where(better, c2n, c2)
            lam = np.where(better, lam * 0.3, lam * 10)
        upd = c2 < best_c2
        best_c2 = np.where(upd, c2, best_c2)
        best_th = np.where(upd[:, None], th, best_th)
    return best_c2, best_th


def local_to_global(S, u, v):
    if isinstance(S, Cylinder):
        ph = S.phi_c + u / S.R
        return np.stack([S.R * np.cos(ph), S.R * np.sin(ph), v], 1)
    return S.P0 + u[:, None] * S.eu + v[:, None] * S.ev


# ================================================================ analytic version
def _dC_dk(s, k):
    """d/dkappa of C = s*sinc(kappa s/2) = (2/kappa) sin(kappa s/2), at fixed s"""
    h = k * s / 2
    out = np.empty_like(h)
    m = np.abs(h) > 1e-4
    hm = h[m]
    out[m] = (s[m] ** 2 / 2) * (np.cos(hm) - np.sin(hm) / hm) / hm
    out[~m] = (s[~m] ** 2 / 2) * (-h[~m] / 3)
    return out


def helix_dtheta(th, s):
    """partial derivatives of the helix point at fixed s: (n,3,4)"""
    phi0, k, z0, t = th.T
    h = k * s / 2; a = phi0 + h
    C = s * sinc(h)
    x, y = C * np.cos(a), C * np.sin(a)
    dC = _dC_dk(s, k)
    D = np.zeros((len(th), 3, 4))
    D[:, 0, 0], D[:, 1, 0] = -y, x
    D[:, 0, 1] = dC * np.cos(a) - y * s / 2
    D[:, 1, 1] = dC * np.sin(a) + x * s / 2
    D[:, 2, 2] = 1.0
    D[:, 2, 3] = s
    return D


def _plane_cross(S, th, s0=None, n_iter=8):
    phi0, k, z0, t = th.T
    if s0 is None:
        d = np.stack([np.cos(phi0), np.sin(phi0), t], 1)
        o = np.stack([np.zeros_like(z0), np.zeros_like(z0), z0], 1)
        den = d @ S.n
        s = np.where(np.abs(den) > 1e-12, ((S.P0 - o) @ S.n) / np.where(np.abs(den) > 1e-12, den, 1), -1.0)
    else:
        s = s0.copy()
    for _ in range(n_iter):
        f = (helix_point(th, s) - S.P0) @ S.n
        fp = helix_tangent(th, s) @ S.n
        fp = np.where(np.abs(fp) > 1e-12, fp, 1e-12)
        s = s - np.clip(f / fp, -(1 + np.abs(s)), 1 + np.abs(s))
    return s


def cross_jac(S, th, s0=None):
    """local (u, v), validity, inside, s, and du/dtheta, dv/dtheta (n,4) - analytic."""
    k = th[:, 1]
    if isinstance(S, Cylinder):
        x = k * S.R / 2
        ok = np.abs(x) < 1
        s = S.R * asin_over(np.where(ok, x, 0.0))
    else:
        s = _plane_cross(S, th, s0)
    p = helix_point(th, s)
    T = helix_tangent(th, s)
    D = helix_dtheta(th, s)
    if isinstance(S, Cylinder):
        pt = p[:, :2]
        dsd = -np.einsum("ni,nia->na", pt, D[:, :2, :]) / np.einsum("ni,ni->n", pt, T[:, :2])[:, None]
        dp = D + T[:, :, None] * dsd[:, None, :]
        phi = np.arctan2(p[:, 1], p[:, 0])
        dphi = (phi - S.phi_c + np.pi) % (2 * np.pi) - np.pi
        u, v = S.R * dphi, p[:, 2]
        du = (p[:, 0, None] * dp[:, 1, :] - p[:, 1, None] * dp[:, 0, :]) / S.R
        dv = dp[:, 2, :]
        inside = ok & (np.abs(dphi) <= S.half_dphi) & (v >= S.z_lo) & (v <= S.z_hi)
    else:
        res = np.abs((p - S.P0) @ S.n)
        ok = (s > 0) & (res < 1e-9 * (1 + np.abs(s))) & (np.abs(k * s) < np.pi)
        dsd = -np.einsum("i,nia->na", S.n, D) / (T @ S.n)[:, None]
        dp = D + T[:, :, None] * dsd[:, None, :]
        q = p - S.P0
        u, v = q @ S.eu, q @ S.ev
        du = np.einsum("i,nia->na", S.eu, dp)
        dv = np.einsum("i,nia->na", S.ev, dp)
        inside = ok & (u >= S.u_lo) & (u <= S.u_hi) & (v >= S.v_lo) & (v <= S.v_hi)
    return u, v, ok, inside, s, du, dv


def model_jac(surfs, th, s0=None):
    """whitened model coordinates (n,2N), validity, inside, s list, J (n,2N,4)"""
    N = len(surfs)
    outs = [cross_jac(S, th, None if s0 is None else s0[i]) for i, S in enumerate(surfs)]
    sig = np.array([S.sigma for S in surfs])
    m = np.concatenate([np.stack([o[0] for o in outs], 1) / sig, np.stack([o[1] for o in outs], 1) / sig], 1)
    J = np.concatenate([np.stack([o[5] for o in outs], 1) / sig[None, :, None],
                        np.stack([o[6] for o in outs], 1) / sig[None, :, None]], 1)
    ok = np.all(np.stack([o[2] for o in outs], 1), 1)
    ins = np.all(np.stack([o[3] for o in outs], 1), 1)
    return m, ok, ins, [o[4] for o in outs], J


def tube_P_fast(surfs, box, kmax, cut, n_samples, rng, batch=400_000):
    N = len(surfs); d = 2 * N - 4
    lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
    tot = tot2 = 0.0; edge = acc = done = 0
    while done < n_samples:
        m = min(batch, n_samples - done)
        th = lo + (hi - lo) * rng.random((m, 4))
        _, ok, ins, _, J = model_jac(surfs, th)
        sel = ins & (np.abs(th[:, 1]) <= kmax)
        gval = np.zeros(m)
        if sel.any():
            Js = J[sel]
            gval[sel] = np.sqrt(np.maximum(np.linalg.det(np.einsum("nia,nib->nab", Js, Js)), 0))
            e = ((th[sel] - lo) / (hi - lo))[:, [0, 2, 3]]
            edge += int(np.sum(np.any((e < 0.002) | (e > 0.998), axis=1))); acc += int(sel.sum())
        tot += gval.sum(); tot2 += (gval ** 2).sum(); done += m
    mean = tot / done; err = math.sqrt(max(tot2 / done - mean ** 2, 0) / done)
    Vd = math.pi ** (d / 2) / math.gamma(d / 2 + 1)
    pref = Vd * cut ** (d / 2) * np.prod([S.sigma ** 2 / S.area for S in surfs]) * np.prod(hi - lo)
    return pref * mean, pref * err, edge / max(acc, 1)


def fit_fast(surfs, U, V, kmax, n_iter=15, k_seeds=(0.0, -0.6, 0.6), th_init=None):
    """Multi-start Levenberg-Marquardt with analytic Jacobian and warm-started crossings."""
    n, N = U.shape
    sig = np.array([S.sigma for S in surfs])
    meas = np.concatenate([U / sig, V / sig], 1)
    if th_init is None:
        P = np.stack([local_to_global(S, U[:, i], V[:, i]) for i, S in enumerate(surfs)], 1)
        r = np.hypot(P[..., 0], P[..., 1]); phi = np.arctan2(P[..., 1], P[..., 0])
        phi0 = np.angle(np.mean(np.exp(1j * phi), 1))
        Sr, Srr, Sz, Srz = r.sum(1), (r * r).sum(1), P[..., 2].sum(1), (r * P[..., 2]).sum(1)
        det = N * Srr - Sr ** 2
        t0 = (N * Srz - Sr * Sz) / det; z00 = (Sz - t0 * Sr) / N
    best_c2 = np.full(n, np.inf); best_th = np.zeros((n, 4))
    for ks in k_seeds:
        th = np.stack([phi0, np.full(n, ks * kmax), z00, t0], 1) if th_init is None else th_init.copy()
        m, ok, _, s, J = model_jac(surfs, th)
        R = meas - m; R[~ok] = 1e6
        c2 = (R ** 2).sum(1); lam = np.full(n, 1e-3)
        for _ in range(n_iter):
            JTJ = np.einsum("nia,nib->nab", J, J)
            g = np.einsum("nia,ni->na", J, R)
            A = JTJ + lam[:, None, None] * (np.eye(4)[None] * np.diagonal(JTJ, axis1=1, axis2=2)[:, :, None]) + 1e-12 * np.eye(4)[None]
            step = np.linalg.solve(A, g[..., None])[..., 0]
            tn = th + step
            mn, okn, _, sn, Jn = model_jac(surfs, tn, s)
            Rn = meas - mn; Rn[~okn] = 1e6
            c2n = (Rn ** 2).sum(1)
            b = c2n < c2
            th = np.where(b[:, None], tn, th); R = np.where(b[:, None], Rn, R)
            J = np.where(b[:, None, None], Jn, J); c2 = np.where(b, c2n, c2)
            s = [np.where(b, a1, a0) for a0, a1 in zip(s, sn)]
            lam = np.where(b, lam * 0.3, lam * 10)
        upd = c2 < best_c2
        best_c2 = np.where(upd, c2, best_c2); best_th = np.where(upd[:, None], th, best_th)
        if th_init is not None:
            break
    return best_c2, best_th
