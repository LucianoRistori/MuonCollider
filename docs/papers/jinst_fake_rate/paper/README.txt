Paper 1 (JINST): "An algorithm-independent lower limit on combinatorial fake
tracks in high-occupancy tracking detectors".

  main.tex            the paper (grey "Plan" boxes are removed as sections are written)
  refs.bib            references (entries marked % CHECK still need verifying)
  jinst-fallback.sty  stand-in for the official class, so the draft builds anywhere
  main.pdf            the current build

Figures are read directly from ../example/ and ../validation/ (made by the
scripts there). Build:

    pdflatex main && bibtex main && pdflatex main && pdflatex main
    rm -f main.aux main.log main.out main.toc main.bbl main.blg

For submission, put the official JINST class files (jinst.cls and JHEP.bst,
from the JINST author instructions at jinst.sissa.it) in this folder: main.tex
then uses them automatically instead of the stand-in.
