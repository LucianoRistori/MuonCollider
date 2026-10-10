"""
exact_helix_fit_core.py

Vectorized nonlinear fit of the EXACT (non-Taylor-truncated) beam-origin-
constrained circular-arc bending view X(Y;b,kappa), for the "true"
4-parameter helix model (b, kappa, p, q) -- replacing the Section 10.4
"conservative" first-order quadratic stand-in with the actual physics.

Exact circle through the origin with initial slope b at Y=0 and signed
curvature kappa=1/R:
    R  = 1/kappa                       (signed)
    Xc = R/sqrt(1+b^2)
    Yc = -b*Xc
    X(Y;b,kappa) = Xc - sign(kappa)*sqrt(R^2 - (Y-Yc)^2)
Verified (sympy + numeric) to reduce to X(Y)=bY+(kappa/2)(1+b^2)^1.5 Y^2
+O(kappa^2) in the small-kappa limit (Section 10.3), and to satisfy
X(0)=0, dX/dY(0)=b for both signs of kappa.

No closed-form minimum exists. A first attempt at alternating 1D
golden-section coordinate descent converged poorly -- b and kappa are
strongly correlated near the minimum (a narrow curved valley), the
classic pathology coordinate descent handles badly. Replaced with
vectorized Levenberg-Marquardt (Gauss-Newton with damping), numerical
Jacobian via central differences, 2x2 normal-equations solved in closed
form per trial -- the right tool for a smooth nonlinear least-squares
problem with 2 correlated parameters, and still fully vectorized across
the trial axis (no per-trial Python loop).
"""

import numpy as np

BIG = 1.0e18


def model_X(Y, b, kappa):
    """Y: (n_planes,), b, kappa: (n_trials,). Returns (n_trials, n_planes),
    and a validity mask (n_trials,) (True where circle reaches every plane)."""
    R = 1.0 / kappa
    Xc = R / np.sqrt(1 + b ** 2)
    Yc = -b * Xc
    s = np.sign(kappa)
    arg = R[:, None] ** 2 - (Y[None, :] - Yc[:, None]) ** 2
    valid = np.all(arg >= 0, axis=1)
    arg_safe = np.clip(arg, 0.0, None)
    model = Xc[:, None] - s[:, None] * np.sqrt(arg_safe)
    return model, valid


def chi2_bend_exact(X, Y, b, kappa, inv_sigma2):
    model, valid = model_X(Y, b, kappa)
    chi2 = np.sum((X - model) ** 2, axis=1) * inv_sigma2
    return np.where(valid, chi2, BIG)


def _lm_single_start(X, Y, inv_sigma2, b_init, kappa_init, kappa_clip,
                      n_iter, lm_lambda0, h_b, h_k):
    n_trials = X.shape[0]
    b = b_init.copy()
    kappa = kappa_init.copy()
    lm_lambda = np.full(n_trials, lm_lambda0)

    def resid(b_, k_):
        model, valid = model_X(Y, b_, k_)
        r = (X - model) * np.sqrt(inv_sigma2)
        return r, valid

    r0, valid0 = resid(b, kappa)
    chi2_0 = np.sum(r0 ** 2, axis=1)

    for it in range(n_iter):
        # numeric Jacobian via central differences (n_trials, n_planes, 2)
        rb_p, _ = resid(b + h_b, kappa)
        rb_m, _ = resid(b - h_b, kappa)
        Jb = (rb_p - rb_m) / (2 * h_b)

        rk_p, _ = resid(b, kappa + h_k)
        rk_m, _ = resid(b, kappa - h_k)
        Jk = (rk_p - rk_m) / (2 * h_k)

        # normal equations: (J^T J + lambda*diag) delta = -J^T r
        Jbb = np.sum(Jb * Jb, axis=1)
        Jkk = np.sum(Jk * Jk, axis=1)
        Jbk = np.sum(Jb * Jk, axis=1)
        gb = np.sum(Jb * r0, axis=1)
        gk = np.sum(Jk * r0, axis=1)

        A11 = Jbb * (1 + lm_lambda)
        A22 = Jkk * (1 + lm_lambda)
        A12 = Jbk
        det = A11 * A22 - A12 * A12
        det = np.where(np.abs(det) < 1e-300, 1e-300, det)
        delta_b = -(A22 * gb - A12 * gk) / det
        delta_k = -(A11 * gk - A12 * gb) / det

        # cap step size in kappa to avoid wild jumps across the domain boundary
        max_step_k = 0.25 * (kappa_clip[1] - kappa_clip[0])
        delta_k = np.clip(delta_k, -max_step_k, max_step_k)
        delta_b = np.clip(delta_b, -2.0, 2.0)

        b_new = b + delta_b
        kappa_new = np.clip(kappa + delta_k, kappa_clip[0], kappa_clip[1])
        kappa_new = np.where(np.abs(kappa_new) < 1e-8,
                              1e-8 * np.sign(kappa_new + 1e-30), kappa_new)

        r_new, valid_new = resid(b_new, kappa_new)
        chi2_new = np.where(valid_new, np.sum(r_new ** 2, axis=1), BIG)

        improved = chi2_new < chi2_0
        b = np.where(improved, b_new, b)
        kappa = np.where(improved, kappa_new, kappa)
        chi2_0 = np.where(improved, chi2_new, chi2_0)
        r0 = np.where(improved[:, None], r_new, r0)
        lm_lambda = np.where(improved, lm_lambda * 0.5, lm_lambda * 2.0)
        lm_lambda = np.clip(lm_lambda, 1e-8, 1e8)

    return b, kappa, chi2_0


def fit_exact_bend(X, Y, inv_sigma2, n_iter=40, lm_lambda0=1e-2,
                    kappa_clip=None, h_b=1e-4, h_k=None, kappa_seeds=None):
    """
    Multi-start vectorized Levenberg-Marquardt fit of (b,kappa) minimizing
    chi2_bend_exact(X,Y,b,kappa). The chi2 surface has local minima (a
    single LM run from one starting point gets trapped, verified by
    noiseless-recovery tests), so this runs LM from several kappa seeds
    per trial and keeps the global-best result elementwise -- essential
    for a MC fake-rate calibration, where an under-converged (too-high)
    chi2 minimum would make the fake rate look artificially low.

    X: (n_trials, n_planes), Y: (n_planes,).
    Returns b_fit, kappa_fit, chi2_fit: each (n_trials,).
    """
    n_trials = X.shape[0]
    Y_max = np.max(np.abs(Y))
    if kappa_clip is None:
        kmax = 3.0 / Y_max
        kappa_clip = (-kmax, kmax)
    if h_k is None:
        h_k = 1e-4 * (kappa_clip[1] - kappa_clip[0])
    if kappa_seeds is None:
        # span the valid kappa range plus the near-zero (straight-line) case
        kappa_seeds = np.concatenate([
            [1e-6],
            np.linspace(kappa_clip[0] * 0.85, kappa_clip[1] * 0.85, 7),
        ])
        kappa_seeds = kappa_seeds[np.abs(kappa_seeds) > 1e-9]

    b_init0 = np.sum(X * Y[None, :], axis=1) / np.sum(Y * Y)

    best_b = np.full(n_trials, np.nan)
    best_k = np.full(n_trials, np.nan)
    best_chi2 = np.full(n_trials, np.inf)

    for k_seed in kappa_seeds:
        b0 = np.clip(b_init0, -7.5, 7.5)
        k0 = np.full(n_trials, k_seed)
        b_fit, k_fit, chi2_fit = _lm_single_start(
            X, Y, inv_sigma2, b0, k0, kappa_clip, n_iter, lm_lambda0, h_b, h_k)
        improved = chi2_fit < best_chi2
        best_b = np.where(improved, b_fit, best_b)
        best_k = np.where(improved, k_fit, best_k)
        best_chi2 = np.where(improved, chi2_fit, best_chi2)

    return best_b, best_k, best_chi2
