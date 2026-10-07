Source of docs/user_manual.pdf (how to use the package: folders, settings
files, ./run_all, ./run_scan, the archive). Keep it up to date when the
workflow changes. To rebuild:

    cd docs/src/user_manual
    pdflatex user_manual.tex && pdflatex user_manual.tex
    cp user_manual.pdf ../../user_manual.pdf
    cp user_manual.pdf ~/Dropbox/Documents/MuonColliderSimulation/Analysis/_START_HERE_user_manual.pdf
