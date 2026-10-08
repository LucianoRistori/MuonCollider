Settings files of the scans of Sections 15-16 (identical to the copies in
their run folders, Analysis/runs/2026-10-07_2057*):

  __cuts_config.txt, __smearing_config.txt, __input_files_config.txt
                                     both scans
  __scan_config_section15_time.txt   Section 15 (time scan)
  __scan_config_section16_angle.txt  Section 16 (angle scan)
  __cuts_config_no_pt_cut.txt        Section 16, the scan without the OT pT cut

To rerun one: save your own settings files from the Analysis folder first,
copy these into it (renaming the scan config to __scan_config.txt, and the
no-pT cuts file to __cuts_config.txt if wanted), run ./run_scan, then put
your own files back.
