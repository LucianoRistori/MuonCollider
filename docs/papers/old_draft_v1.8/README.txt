The original paper draft (v1.4-v1.8, 6-10 October 2026), superseded by the
two papers in docs/papers/jinst_fake_rate (paper 1, the method) and
docs/papers/muc_tracking (paper 2, the Muon Collider application).
Kept frozen as the record; it was in docs/paper/ until 13 October 2026, and the
text of paper_draft.md still refers to that old location.

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
  cd docs/papers/old_draft_v1.8
  python3 make_figures.py                                    # Sections 15-16
  cd figures
  python3 ../../../../fake_rate_framework/section13_density_resolution.py
  python3 ../../../../fake_rate_framework/section14_real_tower.py

Rebuild the PDF:
  cd docs/papers/old_draft_v1.8
  pandoc paper_draft.md -o paper_draft.pdf --pdf-engine=xelatex \
         -V geometry:margin=1in -V colorlinks=true --resource-path=.

History: v1.4 (2026-10-06) is the first commit of this file, copied from
the archived "Track Fitting Simulation" Claude project; git log shows every
change since.
