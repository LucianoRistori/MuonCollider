Source of docs/user_manual.pdf (how to use the package: folders, settings
files, ./run_all, ./run_scan, the archive). Keep it up to date when the
workflow changes. To rebuild:

    cd docs/src/user_manual
    pdflatex user_manual.tex && pdflatex user_manual.tex
    cp user_manual.pdf ../../user_manual.pdf
    rm -f user_manual.aux user_manual.log user_manual.out user_manual.pdf

then commit docs/user_manual.pdf (the master copy). The copy in the
Analysis folder, _START_HERE_user_manual.pdf, is refreshed automatically
the next time ./run_all or ./run_scan is run - never edit it by hand.
