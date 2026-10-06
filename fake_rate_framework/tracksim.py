"""
tracksim.py

Simplified simulation of track fitting in a particle detector made of
N parallel planes (perpendicular to Z), each measuring (X, Y) of a hit.

Phase 1 (current): no magnetic field -> straight 3D lines.
Later phases (planned, not yet implemented): B field -> helical tracks.

Units: millimeters (mm) for all lengths, unless noted otherwise.
"""

from dataclasses import dataclass, field
from typing import List, Tuple, Optional
import numpy as np


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

@dataclass
class Plane:
    """A single measurement plane, perpendicular to Z, centered on the beam axis."""
    z: float          # position along Z [mm]
    size_x: float      # full width in X [mm]
    size_y: float      # full width in Y [mm]
    sigma_x: float      # measurement resolution in X [mm]
    sigma_y: float      # measurement resolution in Y [mm]

    def sample_uniform_xy(self, rng: np.random.Generator) -> Tuple[float, float]:
        """Draw a single (x, y) uniformly over the plane's active area."""
        x = rng.uniform(-self.size_x / 2.0, self.size_x / 2.0)
        y = rng.uniform(-self.size_y / 2.0, self.size_y / 2.0)
        return x, y


@dataclass
class Detector:
    """A stack of parallel planes. Z positions are hardwired here."""
    planes: List[Plane] = field(default_factory=list)

    @property
    def n_planes(self) -> int:
        return len(self.planes)

    @classmethod
    def make_uniform(
        cls,
        n_planes: int,
        z_positions: Optional[List[float]] = None,
        size_x: float = 1000.0,
        size_y: float = 1000.0,
        sigma_x: float = 0.010,   # 10 micron, in mm
        sigma_y: float = 0.010,
        z_spacing: float = 100.0,
    ) -> "Detector":
        """
        Convenience constructor: N identical planes.
        If z_positions is not given, planes are hardwired equally spaced
        by z_spacing, starting at z=0.
        """
        if z_positions is None:
            z_positions = [i * z_spacing for i in range(n_planes)]
        assert len(z_positions) == n_planes

        planes = [
            Plane(z=z, size_x=size_x, size_y=size_y, sigma_x=sigma_x, sigma_y=sigma_y)
            for z in z_positions
        ]
        return cls(planes=planes)


# ---------------------------------------------------------------------------
# Hit generation
# ---------------------------------------------------------------------------

def generate_one_random_hit_per_plane(
    detector: Detector, rng: np.random.Generator
) -> np.ndarray:
    """
    Pick exactly one uniformly-random (x, y) hit on each plane.
    Returns an array of shape (n_planes, 2): columns are (x, y).
    """
    hits = np.empty((detector.n_planes, 2))
    for i, plane in enumerate(detector.planes):
        hits[i, 0], hits[i, 1] = plane.sample_uniform_xy(rng)
    return hits


# ---------------------------------------------------------------------------
# Straight-line fit (no B field): x(z) = a + b*z, y(z) = c + d*z
# ---------------------------------------------------------------------------

def _weighted_linear_fit(z: np.ndarray, y: np.ndarray, w: np.ndarray) -> Tuple[float, float]:
    """
    Weighted least-squares fit of y = a + b*z.
    w is the per-point weight, typically 1/sigma**2.
    Returns (a, b).
    """
    Sw = np.sum(w)
    Swz = np.sum(w * z)
    Swy = np.sum(w * y)
    Swzz = np.sum(w * z * z)
    Swzy = np.sum(w * z * y)

    denom = Sw * Swzz - Swz ** 2
    b = (Sw * Swzy - Swz * Swy) / denom
    a = (Swy - b * Swz) / Sw
    return a, b


@dataclass
class LineFitResult:
    a: float   # x-intercept
    b: float   # dx/dz slope
    c: float   # y-intercept
    d: float   # dy/dz slope
    chi2: float
    ndof: int

    def x_at(self, z):
        return self.a + self.b * z

    def y_at(self, z):
        return self.c + self.d * z


def fit_straight_line_3d(detector: Detector, hits: np.ndarray) -> LineFitResult:
    """
    Fit a 3D straight line to one hit per plane (no B field).
    x(z) and y(z) are fit independently by weighted least squares.

    hits: array of shape (n_planes, 2), columns (x, y), one row per plane
          (same order as detector.planes).
    """
    n = detector.n_planes
    z = np.array([p.z for p in detector.planes])
    sigma_x = np.array([p.sigma_x for p in detector.planes])
    sigma_y = np.array([p.sigma_y for p in detector.planes])

    x = hits[:, 0]
    y = hits[:, 1]

    wx = 1.0 / sigma_x ** 2
    wy = 1.0 / sigma_y ** 2

    a, b = _weighted_linear_fit(z, x, wx)
    c, d = _weighted_linear_fit(z, y, wy)

    resid_x = x - (a + b * z)
    resid_y = y - (c + d * z)
    chi2 = np.sum(wx * resid_x ** 2) + np.sum(wy * resid_y ** 2)

    # 2*n measurements (x_i, y_i for each of n planes), 4 free parameters (a,b,c,d)
    ndof = 2 * n - 4

    return LineFitResult(a=a, b=b, c=c, d=d, chi2=chi2, ndof=ndof)


# ---------------------------------------------------------------------------
# Convenience: one full trial (generate + fit)
# ---------------------------------------------------------------------------

def one_random_combination_chi2(detector: Detector, rng: np.random.Generator) -> float:
    """Pick one random hit per plane and return the chi2 of the straight-line fit."""
    hits = generate_one_random_hit_per_plane(detector, rng)
    result = fit_straight_line_3d(detector, hits)
    return result.chi2


# ---------------------------------------------------------------------------
# Vectorized batch version (fast: no Python loop over trials).
# Needed to reach the many-millions-to-billions of trials required to
# populate the low-chi2 tail well enough to calibrate its power-law shape.
# ---------------------------------------------------------------------------

def batch_random_chi2(detector: Detector, n_trials: int, rng: np.random.Generator) -> np.ndarray:
    """
    Vectorized equivalent of calling one_random_combination_chi2 n_trials times.
    Generates n_trials independent random combinations (one uniform hit per
    plane each) and fits/returns their chi2 values as an array of length
    n_trials. Much faster than a Python loop because everything is done with
    numpy array operations across all trials at once.
    """
    n = detector.n_planes
    z = np.array([p.z for p in detector.planes])
    sigma_x = np.array([p.sigma_x for p in detector.planes])
    sigma_y = np.array([p.sigma_y for p in detector.planes])
    size_x = np.array([p.size_x for p in detector.planes])
    size_y = np.array([p.size_y for p in detector.planes])

    # x, y shape: (n_trials, n_planes)
    x = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * size_x[None, :]
    y = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * size_y[None, :]

    wx = 1.0 / sigma_x ** 2  # shape (n,)
    wy = 1.0 / sigma_y ** 2

    def batch_fit_chi2(coord: np.ndarray, w: np.ndarray) -> np.ndarray:
        # coord: (n_trials, n), w: (n,), z: (n,) -- all planes share the same z, w here
        Sw = np.sum(w)
        Swz = np.sum(w * z)
        Swzz = np.sum(w * z * z)
        Swc = np.sum(w[None, :] * coord, axis=1)         # (n_trials,)
        Swzc = np.sum(w[None, :] * z[None, :] * coord, axis=1)  # (n_trials,)

        denom = Sw * Swzz - Swz ** 2
        b = (Sw * Swzc - Swz * Swc) / denom
        a = (Swc - b * Swz) / Sw

        resid = coord - (a[:, None] + b[:, None] * z[None, :])
        chi2 = np.sum(w[None, :] * resid ** 2, axis=1)
        return chi2

    chi2_x = batch_fit_chi2(x, wx)
    chi2_y = batch_fit_chi2(y, wy)
    return chi2_x + chi2_y


# ---------------------------------------------------------------------------
# General weighted polynomial batch fit (degree >= 1).
#
# Phase 2: low-curvature ("first order in 1/R") approximation to a helical
# track. To first order in curvature, a helix looks like a straight line
# plus a small quadratic-in-z correction:
#   x(z) ~ a + b*z + c*z^2,   y(z) ~ d + e*z + f*z^2
# with c,f related to the curvature (kappa = 1/R) and the bend direction.
# Here we fit x(z) and y(z) as INDEPENDENT quadratics (c and f free, not
# constrained to a single shared kappa/direction) -- a 6-parameter model,
# one parameter more than the true 5-parameter helix, but still exactly
# LINEAR in the fit coefficients, so it keeps the closed-form / vectorized
# machinery below (design matrix [1, z, z^2, ...] instead of [1, z]).
#
# degree=1 reproduces batch_random_chi2 above (ndof = 2N-4).
# degree=2 is the first-order-in-curvature ("parabolic helix") model
# (ndof = 2N-6, since each coordinate now has 3 free parameters).
# ---------------------------------------------------------------------------

def batch_random_poly_fit(
    detector: Detector,
    n_trials: int,
    rng: np.random.Generator,
    degree: int = 2,
    return_coeffs: bool = False,
):
    """
    Vectorized weighted polynomial fit of degree `degree` to x(z) and y(z)
    independently, for n_trials independent random hit combinations (one
    uniform hit per plane each).

    Returns chi2 (array, shape (n_trials,)) if return_coeffs is False.
    If return_coeffs is True, returns (chi2, coeffs_x, coeffs_y) where
    coeffs_x, coeffs_y have shape (n_trials, degree+1), lowest order first
    (coeffs_x[:,0]=a, coeffs_x[:,1]=b, coeffs_x[:,2]=c, ...).
    """
    n = detector.n_planes
    ndeg = degree + 1
    z = np.array([p.z for p in detector.planes])
    sigma_x = np.array([p.sigma_x for p in detector.planes])
    sigma_y = np.array([p.sigma_y for p in detector.planes])
    size_x = np.array([p.size_x for p in detector.planes])
    size_y = np.array([p.size_y for p in detector.planes])

    x = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * size_x[None, :]
    y = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * size_y[None, :]

    wx = 1.0 / sigma_x ** 2
    wy = 1.0 / sigma_y ** 2

    # Design matrix powers of z: Z[k, i] = z_i^k, k=0..degree
    Zpow = np.stack([z ** k for k in range(ndeg)], axis=0)  # (ndeg, n)

    def batch_fit(coord: np.ndarray, w: np.ndarray):
        # Normal-equation matrix M[j,k] = sum_i w_i z_i^j z_i^k -- same for
        # every trial (depends only on geometry), so build/invert it once.
        M = np.einsum("i,ji,ki->jk", w, Zpow, Zpow)  # (ndeg, ndeg)
        Minv = np.linalg.inv(M)

        # RHS vector per trial: R[trial, j] = sum_i w_i z_i^j coord[trial,i]
        R = np.einsum("i,ji,ti->tj", w, Zpow, coord)  # (n_trials, ndeg)

        coeffs = R @ Minv.T  # (n_trials, ndeg), since M symmetric Minv.T=Minv

        # model[trial, i] = sum_j coeffs[trial,j] * z_i^j
        model = coeffs @ Zpow  # (n_trials, n)
        resid = coord - model
        chi2 = np.sum(w[None, :] * resid ** 2, axis=1)
        return chi2, coeffs

    chi2_x, coeffs_x = batch_fit(x, wx)
    chi2_y, coeffs_y = batch_fit(y, wy)
    chi2 = chi2_x + chi2_y

    if return_coeffs:
        return chi2, coeffs_x, coeffs_y
    return chi2


# ---------------------------------------------------------------------------
# Vertex-constrained (origin-through) straight-line fit:
#   x(z) = b*z,   y(z) = e*z          (NO intercept -- track forced through
# the origin, i.e. d0=0 exactly). Only 2 free parameters total (b, e), vs.
# 4 for the general straight-line fit -- so ndof = 2*N - 2, not 2*N - 4.
#
# Still exactly LINEAR in the data (b = sum(w*z*x)/sum(w*z^2), a single
# weighted regression coefficient with no intercept), so the same "WLS
# residual is an exactly linear function of the data" argument applies,
# and the exact universal scaling law chi2/(W/sigma)^2 = s survives
# unchanged for this model too.
# ---------------------------------------------------------------------------

def batch_random_origin_chi2(
    detector: Detector,
    n_trials: int,
    rng: np.random.Generator,
    return_coeffs: bool = False,
):
    """
    Vectorized origin-constrained (vertex-constrained, d0=0) straight-line
    fit: x(z)=b*z, y(z)=e*z, weighted least squares, NO intercept.
    ndof = 2*n_planes - 2.
    """
    n = detector.n_planes
    z = np.array([p.z for p in detector.planes])
    sigma_x = np.array([p.sigma_x for p in detector.planes])
    sigma_y = np.array([p.sigma_y for p in detector.planes])
    size_x = np.array([p.size_x for p in detector.planes])
    size_y = np.array([p.size_y for p in detector.planes])

    x = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * size_x[None, :]
    y = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * size_y[None, :]

    wx = 1.0 / sigma_x ** 2
    wy = 1.0 / sigma_y ** 2

    def batch_fit_origin(coord, w):
        Swzz = np.sum(w * z * z)  # scalar, same for every trial
        Swzc = np.sum(w[None, :] * z[None, :] * coord, axis=1)  # (n_trials,)
        b = Swzc / Swzz
        resid = coord - b[:, None] * z[None, :]
        chi2 = np.sum(w[None, :] * resid ** 2, axis=1)
        return chi2, b

    chi2_x, b_coeff = batch_fit_origin(x, wx)
    chi2_y, e_coeff = batch_fit_origin(y, wy)
    chi2 = chi2_x + chi2_y

    if return_coeffs:
        return chi2, b_coeff, e_coeff
    return chi2
