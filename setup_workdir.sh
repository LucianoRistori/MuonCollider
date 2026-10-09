#!/bin/bash
# =====================================================================
# setup_workdir.sh - create (or refresh) a working folder for the analysis.
#
#     /path/to/MuonCollider/setup_workdir.sh  <working folder>
#
# The working folder is where you edit the settings files and run
# ./run_all and ./run_scan. Its PARENT folder must hold the input files:
#
#     <parent>/Data/       the ROOT n-tuples (BIB and signal)
#     <parent>/Geometry/   the detector-geometry XML files
#     <parent>/<working folder>/   made by this script
#
# What it does:
#   1. creates the working folder if needed;
#   2. copies the four settings files from templates/ (never overwrites
#      an existing one - your settings are kept);
#   3. writes the launchers run_all and run_scan, pointing to THIS copy of
#      the code and to the Python in use now (python3 on your PATH, e.g.
#      that of an active virtual environment, or $PYTHON);
#   4. checks Python and its packages, and the input files named in
#      __input_files_config.txt (only warns: you can add them later).
# Run it again at any time (e.g. after moving the code or the virtual
# environment): it only rewrites the launchers.
# =====================================================================
set -u
say() { printf '%s\n' "$*"; }
stop() { printf '\nERROR: %s\n\n' "$*" >&2; exit 1; }
CODE_DIR="$(cd "$(dirname "$0")" && pwd -P)"

case "${1:-}" in ""|-h|--help) sed -n '3,/^# =====/{/^# =====/d;p;}' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;; esac
[ $# -eq 1 ] || stop "give exactly one argument: the working folder"

# --- Python
PY="${PYTHON:-$(command -v python3)}"
[ -n "$PY" ] || stop "python3 not found - install Python 3.9 or newer (see the user manual)"
case "$PY" in /*) ;; *) PY="$(command -v "$PY")" || stop "'$PY' not found" ;; esac
PY="$(cd "$(dirname "$PY")" && pwd)/$(basename "$PY")"     # absolute, venv link kept
"$PY" -c 'import sys; sys.exit(sys.version_info < (3, 9))' \
    || stop "$PY is Python $("$PY" -c 'import sys; print(sys.version.split()[0])'); 3.9 or newer is needed"

# --- the folder and the settings files
mkdir -p "$1" || stop "cannot create $1"
WORKDIR="$(cd "$1" && pwd -P)"
SIM_DIR="$(dirname "$WORKDIR")"
say "Working folder: $WORKDIR"
say "Code:           $CODE_DIR"
say "Python:         $PY ($("$PY" -c 'import sys; print(sys.version.split()[0])'))"
say ""
for f in "$CODE_DIR"/templates/__*_config.txt; do
    b="$(basename "$f")"
    if [ -f "$WORKDIR/$b" ]; then
        say "  kept     $b (already there)"
    else
        cp "$f" "$WORKDIR/$b" && chmod u+w "$WORKDIR/$b" && say "  created  $b (from templates/)"
    fi
done
for l in run_all run_scan; do
    sed -e "s#@CODE_DIR@#$CODE_DIR#g" -e "s#@PYTHON@#$PY#g" "$CODE_DIR/templates/$l.in" > "$WORKDIR/$l" \
        || stop "cannot write $WORKDIR/$l"
    chmod +x "$WORKDIR/$l"
    say "  wrote    $l"
done
mkdir -p "$WORKDIR/runs"

# --- checks (warnings only)
say ""
problems=0
missing=""
for pkg in numpy uproot awkward matplotlib; do
    "$PY" -c "import $pkg" >/dev/null 2>&1 || missing="$missing $pkg"
done
if [ -n "$missing" ]; then
    problems=1
    say "WARNING: Python packages missing:$missing"
    say "         install them with:  $PY -m pip install -r $CODE_DIR/requirements.txt"
else
    say "Python packages: OK"
fi
for d in Data Geometry; do
    [ -d "$SIM_DIR/$d" ] || { problems=1; say "WARNING: $SIM_DIR/$d/ does not exist yet - create it and put the input files there"; }
done
if [ -d "$SIM_DIR/Data" ] && [ -d "$SIM_DIR/Geometry" ] && [ -z "$missing" ]; then
    if "$PY" "$CODE_DIR/check_configs.py" --quiet "$WORKDIR/__cuts_config.txt" "$WORKDIR/__smearing_config.txt" \
           "$WORKDIR/__input_files_config.txt" --sim-dir "$SIM_DIR" >/dev/null 2>&1; then
        say "Settings and input files: OK"
    else
        problems=1
        say "WARNING: the settings or input files have problems; details with:"
        say "         $PY $CODE_DIR/check_configs.py $WORKDIR/__cuts_config.txt $WORKDIR/__smearing_config.txt $WORKDIR/__input_files_config.txt --sim-dir $SIM_DIR"
        say "         (typically: an input file named in __input_files_config.txt is not in Data/ or Geometry/ yet)"
    fi
fi
say ""
if [ "$problems" = 0 ]; then
    say "Ready. Next:  cd \"$WORKDIR\" && ./run_all"
else
    say "Fix the warnings above, then:  cd \"$WORKDIR\" && ./run_all"
fi
