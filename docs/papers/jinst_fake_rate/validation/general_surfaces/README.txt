General surfaces: the fake-probability formula for arbitrary detector layers
===========================================================================

For any surfaces on which hits measure two local coordinates (coaxial cylinders,
planes at any position and tilt, endcap disks, and mixtures), and any noise
density rho_i on each surface, the expected number of passing combinations on a
set S of layers is, to leading order,

    E_S = V_d c^(d/2) Int dtheta sqrt(det J^T J) prod_{i in S} rho_i(x_i(theta)) sigma_u,i sigma_v,i

theta = (phi0, kappa, z0, cot theta) of a helix from the beam line, integrated over
the tracks that cross the active area of every layer of S; J = derivatives of the
whitened local crossing coordinates. For uniform noise (rho_i = n_i/A_i) this is
P_S prod n_i, with P_S as computed by tube_P_fast.

  gsurf.py      surfaces (Cylinder, Plane), helix crossings with analytic derivatives,
                the integral (tube_P_fast) and a multi-start LM helix fit (fit_fast)
  geoms.py      the four test geometries (4 surfaces, sigma = 0.1 mm, size L in mm)
  run_val.py    formula vs direct Monte Carlo (random combinations, every one fitted)
  check_tower.py  cross-check: the flat projective tower of paper 1 (ratio 0.998 +- 0.007)
  make_fig.py   -> validation_general_surfaces.png/.pdf
  v_*.log       the results (JSON lines)

Results (12 October 2026): P_direct / P_formula against detector size; tends to 1
as the detector grows (finite-size corrections as in paper 1):
  cylinders             0.34, 0.75, 0.81, 0.95+-0.05, 1.17+-0.10  (L = 1, 2, 4, 8, 12)
  planes tilted 45 deg  0.72, 0.98+-0.05, 1.04+-0.05, 0.99+-0.08  (L = 1, 2, 4, 8)
  planes tilted 60 deg  0.69, 1.01+-0.05, 1.09+-0.05, 1.06+-0.09  (L = 1, 2, 4, 8)
  2 cylinders + 2 disks 0.38, 0.79, 0.87+-0.04, 0.91+-0.08         (L = 0.5, 1, 2, 4)
