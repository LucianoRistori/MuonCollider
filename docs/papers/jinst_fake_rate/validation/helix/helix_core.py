"""
helix_core.py - chi2 of the exact helix fit for tracks from the beam line,
vectorized over many hit combinations (straight-line counterpart:
validate_line.py).

Planes at depths (radii) Y_i, flat, measuring X (bending direction) and
Z (along the field), resolution sigma in both. Track model:
  bending view: circle through the origin, X(Y; b, kappa), exact
                (fit: multi-start Levenberg-Marquardt, exact_helix_fit_core)
  depth view:   Z = p + q r,  r = sqrt(X^2 + Y^2)  (linear least squares)
ndof = (N - 2) + (N - 2) = 2N - 4. Acceptance: chi2 < cut and
|kappa_fit| <= kappa_max, with kappa_max given through the dimensionless
xi = kappa_max * Y_N, which is invariant when the detector is shrunk.
Everything is scale-equivariant: shrinking all lengths by lambda at fixed
sigma multiplies chi2 by lambda^2 and kappa by lambda.
"""
import numpy as np
from exact_helix_fit_core import fit_exact_bend

SEEDS = np.array([1e-6, -0.3, -0.12, 0.12, 0.3])   # x 7.2/Y_max, as in the paper's calibration


def chi2_helix(X, Z, Y, sigma, n_iter=22):
    """X, Z: (m, N) hit coordinates (mm); Y: (N,) plane depths (mm).
    Returns chi2_total, kappa_fit (m,)."""
    Y = np.asarray(Y, float)
    inv_s2 = 1.0 / sigma ** 2
    Ymax = np.max(np.abs(Y))
    seeds = SEEDS * (7.2 / Ymax)
    _, kappa, chi2_b = fit_exact_bend(X, Y, inv_s2, n_iter=n_iter, kappa_seeds=seeds)
    r = np.sqrt(X ** 2 + Y[None, :] ** 2)
    n = Y.size
    Sz, Szz = Z.sum(1), (Z * r).sum(1)
    Sr, Srr = r.sum(1), (r * r).sum(1)
    det = n * Srr - Sr ** 2
    q = (n * Szz - Sr * Sz) / det
    p = (Sz - q * Sr) / n
    chi2_d = (((Z - p[:, None] - q[:, None] * r) ** 2).sum(1)) * inv_s2
    return chi2_b + chi2_d, kappa


def chi2_depth(X, Z, Y, sigma):
    r = np.sqrt(X ** 2 + Y[None, :] ** 2)
    n = Y.size
    Sz, Szr = Z.sum(1), (Z * r).sum(1)
    Sr, Srr = r.sum(1), (r * r).sum(1)
    det = n * Srr - Sr ** 2
    q = (n * Szr - Sr * Sz) / det
    p = (Sz - q * Sr) / n
    return (((Z - p[:, None] - q[:, None] * r) ** 2).sum(1)) / sigma ** 2


class Tower:
    """Projective tower of N planes at Y = Y1 * (1, 2, ..., N), sides
    W_i = shape * Y_i (shape = 1000/1486: a 1 m outermost plane at
    1.486 m, as in the paper's example), resolution sigma. Its size is
    set by the smallest plane in resolutions, w_min = W_1 / sigma; all
    lengths scale together, so towers of different w_min are exact
    rescalings of each other. kappa_max = xi / Y_N (xi = 0.2: about
    10 GeV/c at 1.5 m in 5 T)."""
    SHAPE = 1000.0 / 1486.0

    def __init__(self, N, w_min, sigma=0.1, xi=0.2, planes=None):
        self.sigma = sigma
        Y1 = w_min * sigma / self.SHAPE
        Yall = Y1 * np.arange(1, N + 1)
        self.planes = list(range(N)) if planes is None else list(planes)
        self.Y = Yall[self.planes]
        self.W = self.SHAPE * self.Y
        self.kappa_max = xi / Yall[-1]
        self.ndof = 2 * len(self.planes) - 4

    def random_hits(self, rng, m):
        X = rng.uniform(-0.5, 0.5, (m, len(self.Y))) * self.W
        Z = rng.uniform(-0.5, 0.5, (m, len(self.Y))) * self.W
        return X, Z

    def chi2(self, X, Z, cut_max):
        """chi2 and kappa of every combination; the helix fit is run only
        where the depth-view chi2 alone is below cut_max (a necessary
        condition for chi2 < cut_max); elsewhere chi2 = inf."""
        c2 = np.full(X.shape[0], np.inf)
        kap = np.full(X.shape[0], np.nan)
        cd = chi2_depth(X, Z, self.Y, self.sigma)
        sel = np.nonzero(cd < cut_max)[0]
        if len(sel):
            c2[sel], kap[sel] = chi2_helix(X[sel], Z[sel], self.Y, self.sigma)
        return c2, kap

    def accepted(self, X, Z, cut):
        c2, k = self.chi2(X, Z, cut)
        return (c2 < cut) & (np.abs(k) <= self.kappa_max), c2, k
