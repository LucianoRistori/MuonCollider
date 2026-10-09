"""
validate_extrapolation.py - does the calibrated power law predict the
passing probability of a real (unshrunk) detector?

For a family of detectors with the same shape (projective, planes at
Y = 100, 200, ..., N*100 mm, sides proportional to Y) and growing size,
the probability that one random combination passes the 3-sigma chi2 cut
is measured DIRECTLY, and compared with the METHOD: K calibrated once on
the detector shrunk until its smallest plane is 2 resolutions wide, then
P = K (c / lambda^2)^(ndof/2). As the detector grows, c/lambda^2 moves
further below the calibration range, so the test probes the
extrapolation that the method relies on.

Usage: python3 validate_extrapolation.py <N planes> <w_min list> [cal trials] [max direct trials]
Writes extrapolation_N<N>.csv.
"""
import csv
import sys
import time

import numpy as np

from validate_line import calibrate_K, cut_for, p_direct, projective


def main():
    N = int(sys.argv[1])
    w_mins = [float(x) for x in sys.argv[2].split(",")]
    n_cal = int(float(sys.argv[3])) if len(sys.argv) > 3 else 400_000_000
    max_dir = int(float(sys.argv[4])) if len(sys.argv) > 4 else 1_000_000_000
    rng = np.random.default_rng(1000 + N)
    sigma = 0.1
    ratio = N                                  # outermost / innermost plane side
    ref = projective(N, 100, 100 * N, 2.0 * sigma * ratio, sigma)   # w_min = 2: the calibration detector
    cut = cut_for(ref.ndof)
    t0 = time.time()
    cal = calibrate_K(ref, 1.0, n_cal, rng)
    print(f"N={N} ndof={ref.ndof} cut={cut:.3f}: K = {cal['K']:.5g} +- {cal['K_err']:.2g} "
          f"(beta = {cal['beta']:.4g}, {cal['n_fit']} entries below c'={cal['c_max']:g}; {time.time()-t0:.0f}s)", flush=True)
    rows = []
    for wm in w_mins:
        geo = projective(N, 100, 100 * N, wm * sigma * ratio, sigma)
        lam = wm / 2.0
        c_prime = cut / lam ** 2
        P_m = cal["K"] * c_prime ** (ref.ndof / 2)
        n_dir = int(min(max(400 / P_m, 2e6), max_dir))
        t0 = time.time()
        P_d, P_d_err, n_pass = p_direct(geo, cut, n_dir, rng)
        r = dict(N=N, ndof=ref.ndof, w_min=wm, lam=lam, c_prime=c_prime, P_method=P_m,
                 P_method_err=P_m * cal["K_err"] / cal["K"], P_direct=P_d, P_direct_err=P_d_err,
                 n_pass=n_pass, n_trials=n_dir, ratio=P_d / P_m if P_m > 0 else np.nan,
                 ratio_err=(P_d_err / P_m) if P_m > 0 else np.nan, K=cal["K"], beta=cal["beta"])
        rows.append(r)
        print(f"  w_min {wm:6.1f}  c'={c_prime:9.4g}  direct {P_d:.4e} +- {P_d_err:.1e} ({n_pass} of {n_dir:.2e})"
              f"  method {P_m:.4e}   direct/method = {r['ratio']:.4f} +- {r['ratio_err']:.4f}  ({time.time()-t0:.0f}s)",
              flush=True)
    with open(f"extrapolation_N{N}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
