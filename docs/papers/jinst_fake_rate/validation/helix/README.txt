Helix validation of the fake-rate method (exact helix fit, tracks from the beam line)
====================================================================================

Counterpart of the straight-line tests in ../ (validate_line.py, fast_kofn.py).
Projective tower of N flat planes at Y = Y1*(1..N), sides W_i = (1000/1486)*Y_i,
sigma = 0.1 mm in both views, acceptance chi2 < one-sided 3-sigma cut AND
|kappa_fit| <= xi/Y_N (xi = 0.2). Towers of different w_min = W_1/sigma are exact
rescalings of each other.

Code
  helix_core.py         chi2 of the exact helix fit (bending: circle through the origin,
                        multi-start LM from exact_helix_fit_core.py; depth: Z = p + q r),
                        Tower class, depth-view prefilter (exact).
  helix_validation.py   calib <N> <trials> | extrap <N> <w list> <max> |
                        count <w_min> <hits> <events> [seed]
  kgeom.py              K FROM GEOMETRY (no calibration Monte Carlo):
                          K = V_d * Vol(M) / prod (W_i/sigma)^2,
                        Vol(M) = volume of the accepted track manifold in whitened hit
                        space = integral of sqrt(det J^T J) over track parameters whose
                        hits lie inside all planes (and |kappa| <= kappa_max).
                        Exact in the limit c -> 0 (tube volume around the manifold).
  kgeom_sets.py         K_geom for N=4 and for 6/6 and the six 5-plane subsets -> helix_kgeom.json
  kgeom_line.py         same for straight lines; checked against ../calib_sets.json
  ext4deep.py, ext6deep.py   direct P at small c' compared with K_geom c'^(d/2)
  pS_precise.py, pS_quick.py high-statistics P_S for the counting configurations
  refit.py              refits of the MC calibration with different ranges/orders
  make_helix_figure.py  -> validation_helix.png

Results (2026-10-09)
  1. Geometric K vs Monte Carlo calibration fits
     straight lines: K_geom / K_fit = 1.01-1.04 for all 7 sets (fits ~2% low)
     helix N=4: K_geom = 5.833e-3; order-2 fit over c'<=4 gave 5.22e-3 (-10%); refits
                over c'<=1 give 5.69-5.72e-3.  N=6: K_geom = 2.043e-6; fits 1.86e-6, and
                refits scatter by +-10%: the fitted K is not reliable for the helix.
  2. Extrapolation (direct P / K_geom c'^(d/2)) -> 1 as c' -> 0:
     N=4: 0.973, 0.986, 1.020, 1.013 (+-0.02-0.03) at c' = 0.070, 0.031, 0.017, 0.0077
          (P down to 3.5e-7); N=6: 1.00, 0.92, 0.86 (+-0.05) at c' = 0.70, 0.52, 0.40
          (P down to 4e-8), consistent with the 1+beta c' trend seen for N=4.
  3. Counting identity (every combination fitted, 6/6 + six 5/6 subsets):
     counted / (sum P_S prod n): 0.974+-0.049 (w8,n5), 1.045+-0.081 (w14,n8),
     0.929+-0.058 (w20,n10; two seeds); combined 0.97 +- 0.03.
  4. Distinct fakes (merge if >= 4 shared hits) / passing combinations vs eta:
     helix 0.226 (eta 2.43), 0.332 (1.86), 0.48 (1.32); exp(-eta/2) = 0.30, 0.39, 0.52.
     Helix lies below exp(-eta/2) by a fraction ~0.1*eta; converges as eta -> 0.
