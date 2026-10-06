"""
real_tower_calibration.py

Fresh fake-rate calibration (all three models: origin-constrained B=0
line, Section-10.4 "conservative" quadratic, and the exact 4-parameter
helix fit) for the NEW real-detector-informed projective tower
requested by Luciano:

  - Radii Y_i taken from the actual Muon Collider tracker geometry
    (IT barrel L0-L2, OT barrel L0-L2):
        Y_true = [164, 354, 554, 819, 1153, 1486] mm
  - Outermost plane rescaled to 1 m x 1 m, all other planes scaled
    down proportionally so the tower is projective from the origin:
        W_true = 1000 * Y_true / Y_true[-1] mm
  - sigma = 100 micron, B = 5 T, p_T > 10 GeV/c, one-sided-3-sigma
    matched chi2 cut -- same defaults as Sections 11-13.

This geometry is NOT evenly spaced (unlike the synthetic tower of
Sections 7-13), so none of the earlier K/C_kappa constants apply --
a completely fresh calibration is required for each model. The
machinery (whole-geometry LAMBDA shrink, deep-left-tail power-law
plateau fit, kappa_max linear suppression law) is otherwise identical
and reused verbatim from barrel_line_origin_comparison.py,
barrel_helix_tower.py and exact_helix_mc.py.

Final step: combine with the REAL, non-uniform per-layer hit
densities (page 5 of the BIB highlights PDF, after cuts) to get
E[#fakes] for this actual configuration, using the general
(non-uniform-n_i) form of the combinatorial identity,
E[#fakes] = P_true(cut) * prod_i(n_i).
"""

import time
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Geometry: real radii, 1 m^2 outermost plane, projective scaling
# ---------------------------------------------------------------------------
N_PLANES = 6
SIGMA = 0.100  # mm

Y_true = np.array([164.0, 354.0, 554.0, 819.0, 1153.0, 1486.0])
W_true = 1000.0 * Y_true / Y_true[-1]

print("Y_true (mm):", np.round(Y_true, 2))
print("W_true (mm):", np.round(W_true, 4))

# Real per-layer hit densities, hits/mm^2, after cuts (page 5 highlights PDF)
# order: IT-L0, IT-L1, IT-L2, OT-L0, OT-L1, OT-L2
RHO_REAL = np.array([0.02758201712269548, 0.00698052377137699, 0.00306547303960387,
                      0.00199638791184234, 0.00108495481364917, 0.000496270894482882])
N_REAL = RHO_REAL * W_true ** 2
print("n_i (real, non-uniform hit counts on rescaled planes):", np.round(N_REAL, 2))

# Whole-geometry shrink (Section 11.3 fix): shrink X,Y,Z together by LAMBDA,
# holding sigma fixed, so chi2_shrunk*LAMBDA^2 is LAMBDA-independent.
TARGET_RATIO_NEAREST = 10.0
LAMBDA = W_true[0] / (TARGET_RATIO_NEAREST * SIGMA)
W_shrunk = W_true / LAMBDA
Y_shrunk = Y_true / LAMBDA
print(f"LAMBDA = {LAMBDA:.4f}")
print("W_shrunk:", np.round(W_shrunk, 4))
print("Y_shrunk:", np.round(Y_shrunk, 4))

B_TESLA = 5.0
PT_MIN_GEV = 10.0
R_MIN_MM = (PT_MIN_GEV / (0.2998 * B_TESLA)) * 1000.0
KAPPA_MAX_TRUE = 1.0 / R_MIN_MM
print(f"\nB={B_TESLA}T, pT_min={PT_MIN_GEV}GeV/c -> R_min={R_MIN_MM/1000:.4f} m, "
      f"kappa_max(true)={KAPPA_MAX_TRUE:.6e} /mm")


# ===========================================================================
# MODEL 1: origin-constrained straight line (B=0 control)
# ===========================================================================
def batch_fit_line(n_trials, rng, W, Y):
    n = len(Y)
    X = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    Z = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    w = 1.0 / SIGMA ** 2

    b = np.sum(X * Y[None, :], axis=1) / np.sum(Y * Y)
    resid_bend = X - b[:, None] * Y[None, :]
    chi2_bend = np.sum(resid_bend ** 2, axis=1) * w

    r = np.sqrt(X ** 2 + Y[None, :] ** 2)
    Sz = np.sum(Z, axis=1); Szz = np.sum(Z * Z, axis=1)
    Sr = np.sum(r, axis=1); Srz = np.sum(r * Z, axis=1)
    det = n * Szz - Sz ** 2
    q = (n * Srz - Sz * Sr) / det
    p = (Sr - q * Sz) / n
    resid_depth = r - (p[:, None] + q[:, None] * Z)
    chi2_depth = np.sum(resid_depth ** 2, axis=1) * w
    return chi2_bend + chi2_depth


def calibrate_line():
    NDOF = 2 * N_PLANES - 3
    POWER = NDOF / 2.0
    print(f"\n=== MODEL 1: origin-constrained line, ndof={NDOF}, power={POWER} ===")

    for lam_test in [LAMBDA, LAMBDA * 7.0]:
        Wt, Yt = W_true / lam_test, Y_true / lam_test
        c2 = batch_fit_line(500_000, np.random.default_rng(1), Wt, Yt)
        s = np.sort(c2)[5000] * lam_test ** 2
        print(f"  scaling check lambda={lam_test:9.2f}  chi2[1pct]*lambda^2 = {s:.6f}")

    N_TRIALS = 40_000_000
    BATCH = 4_000_000
    rng = np.random.default_rng(20261002)
    chi2_all = np.empty(N_TRIALS, dtype=np.float64)
    done = 0
    t0 = time.time()
    while done < N_TRIALS:
        nb = min(BATCH, N_TRIALS - done)
        chi2_all[done:done + nb] = batch_fit_line(nb, rng, W_shrunk, Y_shrunk)
        done += nb
    print(f"  generated {N_TRIALS:,} trials in {time.time()-t0:.0f}s")

    chi2_sorted = np.sort(chi2_all)
    n = len(chi2_sorted)
    event_counts = np.unique(np.geomspace(20, int(0.02 * n), 50).astype(int))
    cuts = chi2_sorted[event_counts]
    p_emp = event_counts / n
    K_local = p_emp / cuts ** POWER
    deep = cuts < np.percentile(cuts, 25)
    w_ = event_counts[deep].astype(float)
    K_LINE = float(np.sum(K_local[deep] * w_) / np.sum(w_))
    print(f"  K_LINE (new geometry) = {K_LINE:.6e}")

    CUT = stats.chi2.isf(stats.norm.sf(3), NDOF)
    p_true = K_LINE * (CUT / LAMBDA ** 2) ** POWER
    print(f"  ndof={NDOF}, matched 3sigma cut={CUT:.4f}, P_true={p_true:.4e}")
    return dict(ndof=NDOF, power=POWER, K=K_LINE, cut=CUT, p_true=p_true, is_line=True)


# ===========================================================================
# MODEL 2: conservative quadratic (independent b,c), with kappa_max cut
# ===========================================================================
def batch_fit_quad(n_trials, rng, W, Y):
    n = len(Y)
    X = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    Z = rng.uniform(-0.5, 0.5, size=(n_trials, n)) * W[None, :]
    w = 1.0 / SIGMA ** 2

    Y1, Y2 = Y, Y ** 2
    S11 = np.sum(Y1 * Y1) * w
    S12 = np.sum(Y1 * Y2) * w
    S22 = np.sum(Y2 * Y2) * w
    RX1 = np.sum(X * Y1[None, :], axis=1) * w
    RX2 = np.sum(X * Y2[None, :], axis=1) * w
    det = S11 * S22 - S12 * S12
    b = (RX1 * S22 - RX2 * S12) / det
    c = (RX2 * S11 - RX1 * S12) / det
    model_bend = b[:, None] * Y1[None, :] + c[:, None] * Y2[None, :]
    chi2_bend = w * np.sum((X - model_bend) ** 2, axis=1)

    r = np.sqrt(X ** 2 + Y[None, :] ** 2)
    Sw = n * w
    Swz = w * np.sum(Z, axis=1)
    Swzz = w * np.sum(Z * Z, axis=1)
    Swr = w * np.sum(r, axis=1)
    Swzr = w * np.sum(Z * r, axis=1)
    denom = Sw * Swzz - Swz ** 2
    q = (Sw * Swzr - Swz * Swr) / denom
    p = (Swr - q * Swz) / Sw
    model_depth = p[:, None] + q[:, None] * Z
    chi2_depth = w * np.sum((r - model_depth) ** 2, axis=1)
    return chi2_bend + chi2_depth, c


def calibrate_quad():
    NDOF = 2 * N_PLANES - 4
    POWER = NDOF / 2.0
    print(f"\n=== MODEL 2: conservative quadratic, ndof={NDOF}, power={POWER} ===")

    for lam_test in [LAMBDA, LAMBDA * 7.0]:
        Wt, Yt = W_true / lam_test, Y_true / lam_test
        c2, _ = batch_fit_quad(500_000, np.random.default_rng(1), Wt, Yt)
        s = np.sort(c2)[5000] * lam_test ** 2
        print(f"  scaling check lambda={lam_test:9.2f}  chi2[1pct]*lambda^2 = {s:.6f}")

    N_TRIALS = 40_000_000
    BATCH = 4_000_000
    rng = np.random.default_rng(20261003)
    all_chi2 = np.empty(N_TRIALS, dtype=np.float64)
    all_kappa_proxy = np.empty(N_TRIALS, dtype=np.float64)
    done = 0
    t0 = time.time()
    while done < N_TRIALS:
        nb = min(BATCH, N_TRIALS - done)
        c2, c = batch_fit_quad(nb, rng, W_shrunk, Y_shrunk)
        all_chi2[done:done + nb] = c2
        all_kappa_proxy[done:done + nb] = 2.0 * c
        done += nb
    print(f"  generated {N_TRIALS:,} trials in {time.time()-t0:.0f}s")

    order = np.argsort(all_chi2)
    sorted_chi2 = all_chi2[order]
    sorted_kappa = all_kappa_proxy[order]
    n = len(sorted_chi2)

    event_counts = np.unique(np.geomspace(20, int(0.02 * n), 50).astype(int))
    cuts = sorted_chi2[event_counts]
    p_emp = event_counts / n
    K_local = p_emp / cuts ** POWER
    deep = cuts < np.percentile(cuts, 25)
    w_ = event_counts[deep].astype(float)
    K_QUAD = float(np.sum(K_local[deep] * w_) / np.sum(w_))
    print(f"  K_QUAD (new geometry) = {K_QUAD:.6e}")

    kappa_test_values = np.percentile(np.abs(all_kappa_proxy), [0.5, 1, 2, 5, 10, 20])
    chi2_cuts_test = [np.percentile(all_chi2, p) for p in [1, 5, 20, 50]]
    f_by_kappa = {}
    for kmax in kappa_test_values:
        ratios = []
        for ccut in chi2_cuts_test:
            mask_k = np.abs(all_kappa_proxy) < kmax
            p_chi2_only = np.mean(all_chi2 < ccut)
            p_joint = np.mean((all_chi2 < ccut) & mask_k)
            ratios.append(p_joint / p_chi2_only)
        f_by_kappa[kmax] = np.mean(ratios)
    kmax_arr = np.array(list(f_by_kappa.keys()))
    f_arr = np.array(list(f_by_kappa.values()))
    C_KAPPA_QUAD = f_arr[0] / kmax_arr[0]
    print(f"  C_KAPPA_QUAD (new geometry) = {C_KAPPA_QUAD:.6e}")

    KAPPA_SHRUNK_MAX = KAPPA_MAX_TRUE * LAMBDA
    F_KAPPA = C_KAPPA_QUAD * KAPPA_SHRUNK_MAX
    CUT = stats.chi2.isf(stats.norm.sf(3), NDOF)
    p_true = F_KAPPA * K_QUAD * (CUT / LAMBDA ** 2) ** POWER
    print(f"  F_kappa={F_KAPPA:.4e}, ndof={NDOF}, cut={CUT:.4f}, P_true={p_true:.4e}")
    return dict(ndof=NDOF, power=POWER, K=K_QUAD, C_kappa=C_KAPPA_QUAD, cut=CUT, p_true=p_true)


# ===========================================================================
# MODEL 3: exact 4-parameter helix fit
# ===========================================================================
def calibrate_exact():
    from exact_helix_fit_core import fit_exact_bend
    NDOF = 2 * N_PLANES - 4
    POWER = NDOF / 2.0
    print(f"\n=== MODEL 3: exact 4-param helix, ndof={NDOF}, power={POWER} ===")
    inv_sigma2 = 1.0 / SIGMA ** 2
    SEEDS5 = np.array([1e-6, -0.3, -0.12, 0.12, 0.3])

    def batch_fit_exact(n_trials, rng, W, Y, n_iter=22):
        X = rng.uniform(-0.5, 0.5, size=(n_trials, len(Y))) * W[None, :]
        Z = rng.uniform(-0.5, 0.5, size=(n_trials, len(Y))) * W[None, :]
        Y_max = np.max(np.abs(Y))
        kmax = 3.0 / Y_max
        kappa_seeds = SEEDS5 * (kmax / (3.0 / 7.2))
        b_fit, k_fit, chi2_bend = fit_exact_bend(X, Y, inv_sigma2, n_iter=n_iter, kappa_seeds=kappa_seeds)

        n_pts = len(Y)
        r = np.sqrt(X ** 2 + Y[None, :] ** 2)
        Sz = np.sum(Z, axis=1); Szz = np.sum(Z * Z, axis=1)
        Sr = np.sum(r, axis=1); Srz = np.sum(r * Z, axis=1)
        det = n_pts * Szz - Sz ** 2
        q = (n_pts * Srz - Sz * Sr) / det
        p = (Sr - q * Sz) / n_pts
        resid_depth = r - (p[:, None] + q[:, None] * Z)
        chi2_depth = np.sum(resid_depth ** 2, axis=1) * inv_sigma2
        return chi2_bend + chi2_depth, chi2_bend, chi2_depth, b_fit, k_fit

    for lam_test in [LAMBDA, LAMBDA * 7.0]:
        Wt, Yt = W_true / lam_test, Y_true / lam_test
        c2, *_ = batch_fit_exact(3000, np.random.default_rng(1), Wt, Yt)
        s = np.sort(c2)[30] * lam_test ** 2
        print(f"  scaling check lambda={lam_test:9.2f}  chi2[1pct]*lambda^2 = {s:.6f}")

    N_TRIALS = 3_000_000
    BATCH = 200_000
    rng = np.random.default_rng(20261004)
    all_chi2 = np.empty(N_TRIALS, dtype=np.float64)
    all_kappa = np.empty(N_TRIALS, dtype=np.float64)
    done = 0
    t0 = time.time()
    while done < N_TRIALS:
        nb = min(BATCH, N_TRIALS - done)
        c2, cb, cd, bf, kf = batch_fit_exact(nb, rng, W_shrunk, Y_shrunk)
        all_chi2[done:done + nb] = c2
        all_kappa[done:done + nb] = kf
        done += nb
        if done % 500_000 < BATCH:
            print(f"    ... {done:,}/{N_TRIALS:,} trials, {time.time()-t0:.0f}s elapsed")
    print(f"  generated {N_TRIALS:,} trials in {time.time()-t0:.0f}s")

    order = np.argsort(all_chi2)
    sorted_chi2 = all_chi2[order]
    sorted_kappa = all_kappa[order]
    n = len(sorted_chi2)

    event_counts = np.unique(np.geomspace(20, int(0.02 * n), 40).astype(int))
    cuts = sorted_chi2[event_counts]
    p_emp = event_counts / n
    K_local = p_emp / cuts ** POWER
    deep = cuts < np.percentile(cuts, 25)
    w_ = event_counts[deep].astype(float)
    K_EXACT = float(np.sum(K_local[deep] * w_) / np.sum(w_))
    print(f"  K_EXACT (new geometry) = {K_EXACT:.6e}")

    kappa_abs = np.abs(all_kappa)
    kappa_pcts = np.percentile(kappa_abs, [0.5, 1, 2, 5, 10, 20])
    chi2_test_cuts_n = [200, 500, 1200, 3000]
    f_results = {}
    for kmax_s in kappa_pcts:
        fs = []
        for ne in chi2_test_cuts_n:
            cut = sorted_chi2[ne]
            mask_joint = (sorted_chi2 <= cut) & (np.abs(sorted_kappa) <= kmax_s)
            p_joint = np.sum(mask_joint) / n
            p_chi2only = ne / n
            fs.append(p_joint / p_chi2only)
        f_results[kmax_s] = np.mean(fs)
    kmax_arr = np.array(list(f_results.keys()))
    f_arr = np.array(list(f_results.values()))
    C_KAPPA_EXACT = f_arr[0] / kmax_arr[0]
    print(f"  C_KAPPA_EXACT (new geometry) = {C_KAPPA_EXACT:.6e}")

    KAPPA_SHRUNK_MAX = KAPPA_MAX_TRUE * LAMBDA
    F_KAPPA = C_KAPPA_EXACT * KAPPA_SHRUNK_MAX
    CUT = stats.chi2.isf(stats.norm.sf(3), NDOF)
    p_true = F_KAPPA * K_EXACT * (CUT / LAMBDA ** 2) ** POWER
    print(f"  F_kappa={F_KAPPA:.4e}, ndof={NDOF}, cut={CUT:.4f}, P_true={p_true:.4e}")
    return dict(ndof=NDOF, power=POWER, K=K_EXACT, C_kappa=C_KAPPA_EXACT, cut=CUT, p_true=p_true)


if __name__ == "__main__":
    res_line = calibrate_line()
    res_quad = calibrate_quad()
    res_exact = calibrate_exact()

    print("\n\n==================== SUMMARY: real-radii tower ====================")
    print(f"Y_true (mm) = {np.round(Y_true,1).tolist()}")
    print(f"W_true (mm) = {np.round(W_true,2).tolist()}")
    print(f"n_i (real density) = {np.round(N_REAL,2).tolist()}")
    prod_n_real = np.prod(N_REAL)
    print(f"prod(n_i) real-density = {prod_n_real:.6e}")

    for name, res in [("line(B=0)", res_line), ("conservative-quad", res_quad), ("exact-helix", res_exact)]:
        p_true = res["p_true"]
        efakes_real = p_true * prod_n_real
        # also report uniform-density equivalent n(E=1), n(E=10) for reference
        n1 = (1.0 / p_true) ** (1.0 / N_PLANES)
        n10 = (10.0 / p_true) ** (1.0 / N_PLANES)
        print(f"\n[{name}] ndof={res['ndof']} P_true={p_true:.4e}")
        print(f"   uniform-density reference: n(E=1)={n1:.0f}  n(E=10)={n10:.0f} hits/layer")
        print(f"   E[#fakes] at REAL non-uniform density = {efakes_real:.4e}")

    np.save('/tmp/real_tower_results.npy', np.array([res_line, res_quad, res_exact], dtype=object))
    print("\nSaved results.")
