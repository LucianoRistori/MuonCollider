The paper draft: "Estimating Combinatorial Fake-Track Rates in a
Multi-Plane Tracking Detector" (paper_draft.md, Pandoc Markdown).

Layout:
  paper_draft.md     the source (edit this)
  paper_draft.pdf    built from it (rebuild after every change, commit both)
  figures/           the figures it includes
  data/              scan results behind Sections 15-16 (see data/README.txt)
  settings/          the settings files of those scans, to rerun them
  make_figures.py    makes the Section 15-16 figures from data/

Rebuild the figures (only when their inputs change):
  cd docs/paper
  python3 make_figures.py                                    # Sections 15-16
  cd figures
  python3 ../../../fake_rate_framework/section13_density_resolution.py
  python3 ../../../fake_rate_framework/section14_real_tower.py

Rebuild the PDF:
  cd docs/paper
  pandoc paper_draft.md -o paper_draft.pdf --pdf-engine=xelatex \
         -V geometry:margin=1in -V colorlinks=true --resource-path=.

History: v1.4 (2026-10-06) is the first commit of this file, copied from
the archived "Track Fitting Simulation" Claude project; git log shows every
change since.
