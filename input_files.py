"""
Input files of the analysis, as listed in __input_files_config.txt (kept in
the working folder, next to __cuts_config.txt and __smearing_config.txt).

File names in the config are relative to the simulation folder
(MuonColliderSimulation, the folder above the working folder): ROOT files
in its Data/ subfolder, detector geometry XML files in its Geometry/
subfolder.

The combined BIB file (bib_combined) is DERIVED from bib_plus + bib_minus
(+ bib_ipp, if given). ensure_combined() checks that it really contains
exactly those files - the same number of entries, and the same run/event
numbers and hit count per entry, in that order - and (re)builds it with
merge_bib_files.merge() when it is missing or doesn't match (e.g. after
new plus/minus/ipp files were put in place).

Command line (used by run_all_steps.sh):
    python3 input_files.py shell   <__input_files_config.txt> <simulation folder>
        print shell assignments PLUS=... MINUS=... IPP=... COMBINED=...
        SIGNAL=... GEOM_MAIN=... GEOM_VERTEX=... GEOM_IT=... GEOM_OT=...
    python3 input_files.py prepare <__input_files_config.txt> <simulation folder>
        make sure the combined BIB file matches its parts, building it if needed
"""
import configparser
import shlex
import sys
from pathlib import Path

DATA_SUBDIR = "Data"
GEOMETRY_SUBDIR = "Geometry"
DATA_KEYS = ("bib_plus", "bib_minus", "bib_ipp", "bib_combined", "signal")
GEOMETRY_KEYS = ("main", "vertex", "inner_tracker", "outer_tracker")
MAY_BE_EMPTY = {"bib_ipp"}   # empty = combined file is plus + minus only


def input_paths(config_path, sim_dir):
    """{key: absolute Path, or None for an empty bib_ipp} for every entry
    of the [data] and [geometry] sections."""
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    with open(config_path) as f:
        cp.read_file(f)
    sim_dir = Path(sim_dir)
    paths = {}
    for key in DATA_KEYS:
        name = cp.get("data", key, fallback="").strip()
        paths[key] = (sim_dir / DATA_SUBDIR / name) if name else None
    for key in GEOMETRY_KEYS:
        name = cp.get("geometry", key, fallback="").strip()
        paths[key] = (sim_dir / GEOMETRY_SUBDIR / name) if name else None
    return paths


def combined_parts(paths):
    return [paths[k] for k in ("bib_plus", "bib_minus", "bib_ipp") if paths[k] is not None]


def _entry_signature(path):
    """(nRun, nEvt, n_hit) of every entry of the file's (latest) HTAtree."""
    import uproot
    from merge_bib_files import _latest_tree
    f = uproot.open(path)
    tree = f[_latest_tree(f)]
    cols = tree.arrays(["nRun", "nEvt", "n_hit"], library="np")
    return [(int(r), int(e), int(n)) for r, e, n in zip(cols["nRun"], cols["nEvt"], cols["n_hit"])]


def combined_matches(paths):
    """(True, "") if the combined file contains exactly its parts, in
    order; else (False, reason)."""
    combined = paths["bib_combined"]
    if not combined.exists():
        return False, "it does not exist yet"
    expected = []
    for part in combined_parts(paths):
        expected += _entry_signature(part)
    try:
        actual = _entry_signature(combined)
    except Exception as e:  # unreadable/corrupt file
        return False, f"it could not be read ({e})"
    if actual != expected:
        return False, (f"it does not match its parts (it has {len(actual)} entries with "
                       f"{sum(a[2] for a in actual):,} hits; the parts have {len(expected)} "
                       f"entries with {sum(e[2] for e in expected):,} hits)")
    return True, ""


def ensure_combined(paths):
    parts = combined_parts(paths)
    names = " + ".join(p.name for p in parts)
    roles = "plus + minus + ipp" if paths["bib_ipp"] is not None else "plus + minus"
    ok, why = combined_matches(paths)
    if ok:
        print(f"{paths['bib_combined'].name} matches its parts ({roles})")
        return
    print(f"Combined BIB file {paths['bib_combined'].name}: {why}.")
    print(f"Building it from {names} (this takes a minute or two)...")
    from merge_bib_files import merge
    merge(parts, paths["bib_combined"])
    ok, why = combined_matches(paths)
    if not ok:
        raise RuntimeError(f"the rebuilt combined file still does not match its parts: {why}")


def main(argv):
    if len(argv) != 3 or argv[0] not in ("shell", "prepare"):
        print(__doc__.split("Command line")[1])
        return 2
    mode, config_path, sim_dir = argv
    paths = input_paths(config_path, sim_dir)
    if mode == "shell":
        names = {"bib_plus": "PLUS", "bib_minus": "MINUS", "bib_ipp": "IPP",
                 "bib_combined": "COMBINED", "signal": "SIGNAL", "main": "GEOM_MAIN",
                 "vertex": "GEOM_VERTEX", "inner_tracker": "GEOM_IT", "outer_tracker": "GEOM_OT"}
        for key, var in names.items():
            value = str(paths[key]) if paths[key] is not None else ""
            print(f"{var}={shlex.quote(value)}")
        return 0
    ensure_combined(paths)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
