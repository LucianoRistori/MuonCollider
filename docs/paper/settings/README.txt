Settings files of the runs behind Sections 14-16 (identical to the copies
in their run folders, Analysis/runs/2026-10-08_1555* and _155807_pdf):

  __cuts_config.txt, __smearing_config.txt, __input_files_config.txt
                                     both scans
  __scan_config_section15_time.txt   Section 15 (time scan)
  __scan_config_section16_angle.txt  Section 16 (angle scan)
  __cuts_config_no_pt_cut.txt        Section 16, the scan without the OT pT cut
  __cuts_config_section14.txt        Section 14, the ./run_all at the operating
                                     point (N-1 cuts at 1 deg from the angle scan)

To rerun one: save your own settings files from the Analysis folder first,
copy these into it (renaming the scan config to __scan_config.txt, and the
no-pT or Section 14 cuts file to __cuts_config.txt if wanted), run
./run_scan (or ./run_all for Section 14), then put your own files back.
