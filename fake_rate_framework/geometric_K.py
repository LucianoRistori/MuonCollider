"""
geometric_K.py - fake probability P(chi2 < cut) of the exact helix fit for the
6-layer IT+OT barrel projective tower, computed from the geometry (no Monte Carlo).

For random (uniform) hits, the combinations with chi2 < c fill a thin tube of
radius sqrt(c) around the track manifold M (the ideal hit points of all accepted
tracks) in whitened hit space. To leading order in c (volume-of-tubes formula,
Hotelling 1939, Weyl 1939):

    P(chi2 < c) = K c^(ndof/2),   K = V_ndof * Vol(M) / prod_i (W_i/sigma)^2,

V_d = pi^(d/2)/Gamma(d/2+1), Vol(M) = integral of sqrt(det J^T J) over the
accepted track parameters. Exact helix from the beam line: bending view circle
through the origin (b, kappa), |kappa| <= kappa_max (p_T threshold); depth view
Z = p + q r. Vol(M) factorizes into a (b, kappa) grid integral times the area of
the (p, q) polygon whose hits lie inside every plane. Validated against direct
Monte Carlo and brute-force counting in docs/papers/jinst_fake_rate/validation/helix.

Computed for 6 of 6 and for every subset of 5 and of 4 planes, at sigma = SIGMA_REF in both
coordinates; P scales as (sigma_u/SIGMA_REF)^nb (sigma_v/SIGMA_REF)^nd with
nb = nd = N - 2 for N planes.

Usage: python3 geometric_K.py      (prints the constants for fake_rate.py)
"""
import itertools
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2dist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exact_helix_fit_core import model_X

P_3SIGMA = 1.3499e-3          # one-sided 3 sigma: P(chi2 > cut) for genuine tracks
SIGMA_REF = 0.100             # mm
Y_TOWER = np.array([164.0, 354.0, 554.0, 819.0, 1153.0, 1486.0])   # IT L0-L2, OT L0-L2 (mm)
W_TOWER = 1000.0 * Y_TOWER / Y_TOWER[-1]                            # 1 m^2 outermost plane
B_TESLA, PT_MIN = 5.0, 10.0
KAPPA_MAX = 0.2998 * B_TESLA / PT_MIN / 1000.0                     # 1/mm


def K_geom(Y, W, sigma, kappa_max, nb=400, nk=400, nq=600, bmax=0.7):
    Y, W = np.asarray(Y, float), np.asarray(W, float)
    N = len(Y); d = 2 * N - 4
    b = np.linspace(-bmax, bmax, nb + 1); b = 0.5 * (b[1:] + b[:-1])
    k = np.linspace(-kappa_max, kappa_max, nk + 1); k = 0.5 * (k[1:] + k[:-1])
    k[np.abs(k) < 1e-15] = 1e-15
    B, Kp = [a.ravel() for a in np.meshgrid(b, k, indexing="ij")]
    X, valid = model_X(Y, B, Kp)
    inside = valid & np.all(np.abs(X) <= W / 2, axis=1)
    B, Kp, X = B[inside], Kp[inside], X[inside]
    hb, hk = 1e-6, kappa_max * 1e-4
    dXb = (model_X(Y, B + hb, Kp)[0] - model_X(Y, B - hb, Kp)[0]) / (2 * hb) / sigma
    dXk = (model_X(Y, B, Kp + hk)[0] - model_X(Y, B, Kp - hk)[0]) / (2 * hk) / sigma
    detA = np.sqrt(np.maximum((dXb ** 2).sum(1) * (dXk ** 2).sum(1) - ((dXb * dXk).sum(1)) ** 2, 0))
    r = np.sqrt(X ** 2 + Y ** 2)
    detD = np.sqrt(N * (r ** 2).sum(1) - r.sum(1) ** 2) / sigma ** 2
    qmax = 1.2 * (W / r.min(0)).max()
    q = np.linspace(-qmax, qmax, nq + 1)
    area = np.zeros(len(B))
    for chunk in np.array_split(np.arange(len(B)), max(1, len(B) // 20000)):
        rr = r[chunk][:, None, :]
        hi = np.min(W / 2 - q[None, :, None] * rr, axis=2)
        lo = np.max(-W / 2 - q[None, :, None] * rr, axis=2)
        area[chunk] = np.trapezoid(np.maximum(hi - lo, 0), q, axis=1)
    vol = np.sum(detA * detD * area) * (b[1] - b[0]) * (k[1] - k[0])
    Vd = math.pi ** (d / 2) / math.gamma(d / 2 + 1)
    return Vd * vol / np.prod((W / sigma) ** 2)


def plane_sets(n=6, k_min=4):
    """All sets of at least k_min of the n planes, largest first."""
    return [S for k in range(n, k_min - 1, -1) for S in itertools.combinations(range(n), k)]


def main():
    print(f"kappa_max = {KAPPA_MAX:.6e} /mm (p_T > {PT_MIN} GeV/c, B = {B_TESLA} T), sigma = {SIGMA_REF} mm")
    print("P_REF_SETS = {")
    for S in plane_sets():
        d = 2 * len(S) - 4
        cut = float(chi2dist.isf(P_3SIGMA, d))
        K = K_geom(Y_TOWER[list(S)], W_TOWER[list(S)], SIGMA_REF, KAPPA_MAX)
        print(f"    {S}: {K * cut ** (d / 2):.5e},   # ndof {d}, K {K:.5e}", flush=True)
    print("}")


if __name__ == "__main__":
    main()
