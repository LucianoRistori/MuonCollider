The 100 ps operating point (14 October 2026).
derive_cuts.py  z0 and time cuts of every subsystem at 98% signal containment (N-1), as
                resolution_scan.py does, for a given smearing; gave the cuts now in
                __cuts_config.txt (100 um, 100 ps, 1.5 deg), and reproduces the old 50 ps cuts.
hist_dump.py    histograms of t_corrected, z0 and curvature (IT/OT barrel, BIB and signal,
                all and N-1) for the figures of paper 2.
Then: ./run_all -> run 2026-10-14_162239_pdf; resolution scans with resolution_scan.py
(log 10 um-1 mm, log 10 ps-1 ns, list 0.1-10 deg); whole-barrel fakes with
../regions/regions_barrel_z0.py (REG_BIB, REG_DENS pointing to that run's BIB passing hits).
