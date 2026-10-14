Paper 2 (JINST): hit-level BIB rejection in the IT+OT barrel of a muon collider and the
resolutions it requires; uses the fake-track lower limit of paper 1.

  main.tex, refs.bib      the source (main.pdf built from it)
  jinst-fallback.sty      stand-in for jinst.cls until the official class is added
  make_figures.py         makes fig_*.pdf/png from inputs/
  inputs/                 the numbers behind the paper, copied from the Analysis folder:
                          run 2026-10-14_162239_pdf (100 um, 100 ps, 1.5 deg; cuts at 98%
                          containment): summary, densities, efficiency, fake_rate.csv;
                          resolution scans (scan_*_results.csv); histograms (hists_op100.npz);
                          whole-barrel fakes from regions_barrel_z0.py (barrel_z0_*.json)
  versions/               frozen PDF of every draft
Build: pdflatex main && bibtex main && pdflatex main && pdflatex main
