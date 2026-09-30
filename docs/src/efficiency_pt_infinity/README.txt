Source of docs/efficiency_pt_infinity.pdf

  "Signal efficiencies in the limit pT -> infinity: how the fit of
   eps_inf + c/pT^2 to the muons above 10 GeV/c works"

The note explains how the signal hit efficiency and the track-finding
efficiency quoted by a run (summary table, highlights PDF) are obtained, with
the numbers of run 2026-09-30_145059_pdf as the example.

To rebuild it, e.g. after changing the cuts:

  1. On the analysis machine (needs the analysis code and the data):

       python3 export_fitdata.py [analysis folder] [--label "run ..."]

     Reads the three settings files in the analysis folder (default
     ~/Dropbox/Documents/MuonColliderSimulation/Analysis) and writes
     fitdata.npz here: for every signal hit and every signal muon, the
     generated pT and whether it passes, exactly as apply_cuts.py and
     track_efficiency.py compute them in a run (it prints a check).
     --label names the settings in the note's subtitle.

  2. python3 ptinf_numbers.py     -> numbers.json (every number in the note)
  3. python3 figures.py           -> fig_tracks.pdf, fig_hits.pdf
  4. python3 make_note.py         -> efficiency_pt_infinity.tex
     latexmk -pdf efficiency_pt_infinity.tex   (or pdflatex, twice)
  5. Copy efficiency_pt_infinity.pdf to docs/.

Steps 2-4 need numpy, matplotlib and LaTeX (pdflatex with the mathpazo,
booktabs, tcolorbox, titlesec, lastpage and enumitem packages; figures.py
uses matplotlib's pgf backend, so the figures share the note's fonts).

Numbers, tables and figures follow the data. A few sentences name the cut
values of run 2026-09-30_145059_pdf (the 2 GeV/c pT cut in the OT only, the
0.07 ns time cut of the VXD barrel, the 0.03 ns time resolution); edit them
in make_note.py if those change.

fitdata.npz (about 1.3 MB) is not kept in git.
