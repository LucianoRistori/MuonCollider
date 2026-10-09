End-to-end validation of the fake-rate method, straight lines (B = 0)
=====================================================================
Paper 1 (JINST), Section 4. Pure synthetic Monte Carlo: no simulation data.

Detector: N planes at Y = 100, 200, ... mm, sides proportional to Y
(a projective tower), two measured coordinates per plane, resolution
0.1 mm; free straight line per view (ndof = 2N - 4); one-sided 3-sigma
chi2 cuts (P(chi2 > cut) = 1.35e-3 for genuine tracks).

Code
  validate_line.py           geometry, whitened chi2, direct P, calibration of K
                             (fit P = K c^(d/2) (1 + beta c + gamma c^2) on the
                             detector shrunk to a smallest plane of 2 resolutions,
                             c' <= 4, Poisson likelihood)
  validate_extrapolation.py  test 1: P measured directly vs the method, for
                             detectors of growing size (P from 2e-2 to 3e-7)
  fastcount.py               exact enumeration of all passing combinations
                             (prunes with the monotonicity of chi2; identical
                             to brute force: check_fastcount.py)
  calib_sets.py              K for the 6-plane set and its six 5-plane subsets
  fast_kofn.py               test 2: "at least 5 of 6", 8 to 1000 hits per plane,
                             counted combinations vs the method, distinct fakes,
                             window occupancy eta
  validate_counting.py       brute-force version of test 2 (small detectors)
  k_systematics.py           stability of K vs fit range and order (about +-1%)
  make_validation_figures.py validation_line.png

Results
  extrapolation_N4.csv, extrapolation_N6.csv   test 1
  kofn_results.csv                              test 2 (150 hits: three runs)
  calib_sets.json, ref_corrections.json         calibration constants
  counting_w*_n*.csv                            brute-force runs

Findings (9-10 October 2026)
  - The counting identity holds: passing combinations counted in noise-only
    events agree with P x prod n for every set of planes; over all
    configurations counted/method = 1.006 +- 0.008, with P per 5-plane
    combination down to ~5e-16.
  - The calibrated power law matches direct Monte Carlo below the
    calibration range (ratio 1 within ~5% statistics down to c' = 0.004);
    the corrections (1 + beta c' + gamma c'^2) describe the approach.
    K needs the quadratic term: a linear one underestimates it by 7%.
  - Distinct fakes (candidates sharing >= 4 hits merged) / passing
    combinations = exp(-eta/2) within a few %, eta = window occupancy
    (other hits that could replace one of a candidate's hits), from 0.32
    at eta = 2.2 to 0.92 at eta = 0.17. eta_i ~ (n_i - 1) pi (c/2)
    (sigma_i^2 + V_i) / A_i agrees with the measured eta within 5-20%
    when the window is small compared with the plane.
