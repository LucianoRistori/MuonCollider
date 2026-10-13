Data behind Sections 14-16 of the paper (v1.7), copied from runs in
~/Dropbox/Documents/MuonColliderSimulation/Analysis/runs/ (each holds the
settings files, log and plots of that run). All use the corrected BIB files
(ntu_bib_*_1bx.root, October 2026).

  section14_density_per_layer.csv    2026-10-08_155807_pdf  (./run_all, step4_cuts/)
  section15_time_scan.csv            2026-10-08_214609_scan_time
  section16_angle_scan.csv           2026-10-08_214636_scan_angle
  section16_angle_scan_no_pt_cut.csv 2026-10-08_214723_scan_angle (pT cut off in the OT)

E[#fakes] uses the fake probability scaled to the 50 um position resolution
of these settings, (0.05/0.1)^ndof (paper v1.7). The v1.6 scans
(2026-10-08_155550, _155612, _155659) used it unscaled, at 0.1 mm: their
E[#fakes] is exactly 256 times higher (512 for the B = 0 line); every other
column is identical. section14_density_per_layer.csv (densities only) is
unchanged.

Settings: position 0.05 mm, time 0.03 ns, angle 1 deg (the one not scanned),
seed 42; pT >= 2 GeV/c in the OT; z0 and time cuts of every subsystem set
to keep 98% of the signal hits in their N-1 distributions (Section 14: those
cuts at 1 deg, from the Section 16 scan). The settings files are in
../settings/ (see ../README.txt to rerun). Versions 1.5 used the older BIB
files: runs 2026-10-07_2057*.
