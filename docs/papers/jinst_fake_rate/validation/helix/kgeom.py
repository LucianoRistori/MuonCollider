"""K from geometry: K = V_d * Vol(track manifold in acceptance, whitened) / prod (W_i/sigma)^2.
Bending (b, kappa) on a grid; depth (p, q) polygon area exact per q-grid."""
import sys, math, numpy as np
from exact_helix_fit_core import model_X
from helix_core import Tower

def K_geom(tower, nb=400, nk=400, nq=600):
    Y, W, s = tower.Y, tower.W, tower.sigma
    N = len(Y); d = 2 * N - 4
    km = tower.kappa_max
    bmax = 0.7
    b = np.linspace(-bmax, bmax, nb + 1); b = 0.5 * (b[1:] + b[:-1])
    k = np.linspace(-km, km, nk + 1); k = 0.5 * (k[1:] + k[:-1])
    k[np.abs(k) < 1e-12] = 1e-12
    B, Kp = [a.ravel() for a in np.meshgrid(b, k, indexing="ij")]
    X, valid = model_X(Y, B, Kp)
    inside = valid & np.all(np.abs(X) <= W / 2, axis=1)
    B, Kp, X = B[inside], Kp[inside], X[inside]
    hb, hk = 1e-6, km * 1e-4
    dXb = (model_X(Y, B + hb, Kp)[0] - model_X(Y, B - hb, Kp)[0]) / (2 * hb) / s
    dXk = (model_X(Y, B, Kp + hk)[0] - model_X(Y, B, Kp - hk)[0]) / (2 * hk) / s
    detA = np.sqrt(np.maximum((dXb**2).sum(1) * (dXk**2).sum(1) - ((dXb*dXk).sum(1))**2, 0))
    r = np.sqrt(X**2 + Y**2)
    detD = np.sqrt(N * (r**2).sum(1) - r.sum(1)**2) / s**2
    # area of {(p,q): |p + q r_i| <= W_i/2 all i}
    qmax = 1.2 * (W / r.min(0)).max()
    area = np.zeros(len(B))
    q = np.linspace(-qmax, qmax, nq + 1)
    for chunk in np.array_split(np.arange(len(B)), max(1, len(B) // 20000)):
        rr = r[chunk][:, None, :]
        hi = np.min(W / 2 - q[None, :, None] * rr, axis=2)
        lo = np.max(-W / 2 - q[None, :, None] * rr, axis=2)
        L = np.maximum(hi - lo, 0)
        area[chunk] = np.trapezoid(L, q, axis=1)
    cell = (b[1] - b[0]) * (k[1] - k[0])
    vol = np.sum(detA * detD * area) * cell
    Vd = math.pi ** (d / 2) / math.gamma(d / 2 + 1)
    return Vd * vol / np.prod((W / s) ** 2)

if __name__ == "__main__":
    N = int(sys.argv[1]); planes = None
    if len(sys.argv) > 2: planes = [int(x) for x in sys.argv[2].split(",")]
    t = Tower(N, 2.0, planes=planes)
    for res in ((200, 200, 300), (400, 400, 600), (800, 800, 1200)):
        print(res, f"{K_geom(t, *res):.6g}", flush=True)
