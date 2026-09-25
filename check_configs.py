"""
Check __cuts_config.txt, __smearing_config.txt and __input_files_config.txt
before a run, and print a short, readable summary of the settings they
contain.

Why this exists: Python's configparser (used to read these files) silently
falls back to a default whenever a section or setting it asks for isn't
there. So a typo - "sigma_u_m" instead of "sigma_u_mm", or "[angle-u]"
instead of "[angle_u]" - raises no error: that value just quietly stays at
its default (usually "no smearing", or "no cut"), and a whole run is
produced with different settings than intended. This script rejects,
before anything is run:
  - unknown or misspelled sections, settings and cut-table rows
    (suggesting the closest valid name), and missing ones,
  - values that aren't a number / true-false (or "off" in the cut table),
  - negative cuts or sigmas,
  - duplicated sections or settings,
  - input files that don't exist,
and warns about settings that are probably not what was meant (e.g. a
sigma set but its section left disabled, or a cut line that would fall
outside its zoomed plot).

The cuts file is read by cuts_table.py - the same code bib_common.load_cuts
uses, so the check and the analysis cannot disagree. The accepted
smearing settings mirror bib_common.load_smearing_config, and the input
file names input_files.py - keep them in sync if a setting is added there.

Usage:
    python3 check_configs.py <__cuts_config.txt> <__smearing_config.txt>
                             [<__input_files_config.txt> --sim-dir <folder>] [--brief]
Exit code 0 = OK (warnings allowed), 1 = at least one error.
--brief prints just a few summary lines (used for runs/<id>/summary.txt);
--quiet prints nothing unless there is a problem.
"""
import configparser
import difflib
import math
import sys
from pathlib import Path

import cuts_table as ct

NUMBER, WHOLE, TRUEFALSE, FILENAME = "number", "whole number", "true/false value", "file name"

# section -> {setting: (kind, required)}   (the cuts file: see cuts_table.py)
SMEAR_SPEC = {
    "position": {"sigma_u_mm": (NUMBER, True), "sigma_v_mm": (NUMBER, True),
                 "enabled": (TRUEFALSE, True)},
    "time": {"sigma_t_ns": (NUMBER, True), "enabled": (TRUEFALSE, True)},
    "angle_u": {"sigma_deg": (NUMBER, True), "enabled": (TRUEFALSE, True)},
    "angle_v": {"sigma_deg": (NUMBER, True), "enabled": (TRUEFALSE, True)},
    "general": {"seed": (WHOLE, True)},
}
# __input_files_config.txt - same keys as input_files.DATA_KEYS / GEOMETRY_KEYS
INPUTS_SPEC = {
    "data": {k: (FILENAME, True) for k in
             ("bib_plus", "bib_minus", "bib_ipp", "bib_combined", "signal")},
    "geometry": {k: (FILENAME, True) for k in
                 ("main", "vertex", "inner_tracker", "outer_tracker")},
}
MAY_BE_EMPTY = {"bib_ipp"}
MUST_NOT_BE_NEGATIVE = {"sigma_u_mm", "sigma_v_mm", "sigma_t_ns", "sigma_deg"}


def _suggest(name, candidates):
    close = difflib.get_close_matches(name, list(candidates), n=1, cutoff=0.6)
    return f" - did you mean '{close[0]}'?" if close else ""


def _one_line(exc):
    return " ".join(str(exc).split())


def read_and_check(path, spec, errors):
    """Parse `path` against `spec`; append problems to `errors`. Returns
    {(section, setting): value} for everything that parsed correctly."""
    label = path.split("/")[-1]
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    try:
        with open(path) as f:
            cp.read_file(f, source=label)
    except FileNotFoundError:
        errors.append(f"{label}: file not found: {path}")
        return {}
    except configparser.Error as e:
        errors.append(f"{label}: {_one_line(e)}")
        return {}

    for section in cp.sections():
        if section not in spec:
            errors.append(f"{label}: unknown section [{section}]{_suggest(section, spec)}")
            continue
        for key in cp[section]:
            if key not in spec[section]:
                errors.append(f"{label}: unknown setting '{key}' in [{section}]"
                              f"{_suggest(key, spec[section])}")

    values = {}
    for section, settings in spec.items():
        if not cp.has_section(section):
            errors.append(f"{label}: missing section [{section}] (to switch something off, "
                          f"set enabled = false instead of deleting its section)")
            continue
        for key, (kind, required) in settings.items():
            if key not in cp[section]:
                if required:
                    errors.append(f"{label}: missing setting '{key}' in [{section}]")
                continue
            raw = cp[section][key]
            if kind == FILENAME:
                name = raw.strip()
                if not name and key not in MAY_BE_EMPTY:
                    errors.append(f"{label}: [{section}] {key} is empty")
                    continue
                values[(section, key)] = name
                continue
            try:
                if kind == NUMBER:
                    val = cp.getfloat(section, key)
                    if not math.isfinite(val):
                        raise ValueError
                elif kind == WHOLE:
                    val = cp.getint(section, key)
                else:
                    val = cp.getboolean(section, key)
            except ValueError:
                errors.append(f"{label}: [{section}] {key} = {raw!r} is not a valid {kind}")
                continue
            if kind != TRUEFALSE and key in MUST_NOT_BE_NEGATIVE and val < 0:
                errors.append(f"{label}: [{section}] {key} = {raw} must not be negative")
                continue
            values[(section, key)] = val
    return values


def _g(x):
    return f"{x:g}"


def track_text(cuts):
    trk = f">= {cuts['track']['min_hits_found']} surviving hits"
    return trk + (", vertex-detector hits not counted" if cuts["track"]["exclude_vertex_hits"]
                  else ", all subsystems counted")


def describe(cuts, smear):
    """Human-readable lines, as (name, text): the cuts (a table, one row per
    subsystem), then the smearing. `cuts` is what cuts_table.read returns."""
    header, rows = ct.table_rows(cuts)
    widths = [max([len(header[i])] + [len(r[i + 1]) for r in rows]) for i in range(len(header))]

    def cells(r):
        return "   ".join(c.rjust(w) for c, w in zip(r, widths))
    cut_lines = [("hits kept if", "|t_corrected| <= time, |z0| <= z0 and pT >= pT, per subsystem:"),
                 ("", cells(header))]
    cut_lines += [(r[0], cells(r[1:])) for r in rows]
    cut_lines.append(("track found", track_text(cuts)))

    def sig(sec, key, unit):
        if not smear.get((sec, "enabled"), False):
            return "off"
        s = smear.get((sec, key), 0.0)
        return f"{_g(s)} {unit}" if s > 0 else "0 (off)"

    if smear.get(("position", "enabled"), False):
        pos = (f"sigma_u = {_g(smear.get(('position', 'sigma_u_mm'), 0.0))} mm, "
               f"sigma_v = {_g(smear.get(('position', 'sigma_v_mm'), 0.0))} mm")
    else:
        pos = "off"
    smear_lines = [
        ("position", pos),
        ("time", "sigma_t = " + sig("time", "sigma_t_ns", "ns") if sig("time", "sigma_t_ns", "ns") != "off" else "off"),
        ("angle", f"sigma_u = {sig('angle_u', 'sigma_deg', 'deg')}, "
                  f"sigma_v = {sig('angle_v', 'sigma_deg', 'deg')}"),
        ("random seed", str(smear.get(("general", "seed"), ""))),
    ]
    return cut_lines, smear_lines


def brief_lines(cuts, smear):
    """Three short lines summarizing the settings (for summary.txt, where
    the results table follows - when the cuts differ between subsystems,
    that table lists each subsystem's cuts)."""
    def on(sec):
        return smear.get((sec, "enabled"), False)

    def val(sec, key):
        s = smear.get((sec, key), 0.0)
        return _g(s) if s > 0 else "0"
    pos = (f"position u/v {val('position', 'sigma_u_mm')}/{val('position', 'sigma_v_mm')} mm"
           if on("position") else "position off")
    tim = f"time {val('time', 'sigma_t_ns')} ns" if on("time") else "time off"
    au = f"{val('angle_u', 'sigma_deg')}" if on("angle_u") else "off"
    av = f"{val('angle_v', 'sigma_deg')}" if on("angle_v") else "off"
    ang = "angle off" if au == av == "off" else f"angle u/v {au}/{av} deg"
    return [
        "Cuts:        " + (", ".join(ct.summary(cuts, c) for c in ct.COLUMNS)
                             if ct.is_uniform(cuts) else "set per subsystem - see the table below"),
        "Track found: " + track_text(cuts),
        "Resolutions: " + ", ".join([pos, tim, ang, f"seed {smear.get(('general', 'seed'), '')}"]),
    ]


def warnings_for(cuts, smear):
    warn = []

    def where(ids):
        ids = sorted(ids)
        return ("all subsystems" if len(ids) == len(ct.SYSTEM_IDS)
                else ", ".join(ct.NAMES[s] for s in ids))

    def zoom_setting(col):
        return (f"zoom_halfwidth in [{ct.CUT_NAMES[col]}]" if cuts["old_format"]
                else f"{ct.TITLES[col]} in [zoom]")
    for col in ("time", "z0"):
        vals, zoom = cuts["cuts"][col], cuts["zoom"][col]
        zero = [s for s, v in vals.items() if v == 0]
        wide = [s for s, v in vals.items() if v is not None and v > zoom]
        if zero:
            warn.append(f"__cuts_config.txt: {ct.TITLES[col]} = 0 rejects essentially every hit "
                        f"({where(zero)})")
        if wide:
            warn.append(f"__cuts_config.txt: the {ct.TITLES[col]} cut lies outside the zoomed N-1 "
                        f"plots' range (+/-{_g(zoom)} {ct.UNITS[col]}) in {where(wide)} - raise "
                        f"{zoom_setting(col)} to see the cut lines there")
    vals, zoom = cuts["cuts"]["pt"], cuts["zoom"]["pt"]
    zero = [s for s, v in vals.items() if v == 0]
    low = [s for s, v in vals.items() if v is not None and 0 < v < zoom]
    if zero:
        warn.append(f"__cuts_config.txt: pT = 0 keeps every hit, the same as off ({where(zero)})")
    if low:
        warn.append(f"__cuts_config.txt: the pT cut is below the zoomed N-1 plots' lower edge "
                    f"({_g(zoom)} GeV/c) in {where(low)} - lower {zoom_setting('pt')} to see the "
                    f"cut lines there")
    for sec, keys in (("position", ("sigma_u_mm", "sigma_v_mm")), ("time", ("sigma_t_ns",)),
                      ("angle_u", ("sigma_deg",)), ("angle_v", ("sigma_deg",))):
        if (sec, "enabled") in smear and not smear[(sec, "enabled")]:
            if any(smear.get((sec, k), 0.0) > 0 for k in keys):
                warn.append(f"__smearing_config.txt: [{sec}] has a nonzero sigma but enabled = false, "
                            f"so it will NOT be applied")
    return warn


def notes_for(cuts):
    if cuts["old_format"]:
        return ["__cuts_config.txt is in the old format (one section per cut), so each cut "
                "applies to every subsystem. For cuts per subsystem, use the table format - "
                "see templates/__cuts_config.txt."]
    return []


def _tilde(path):
    home = str(Path.home())
    s = str(path)
    return "~" + s[len(home):] if s == home or s.startswith(home + "/") else s


def check_input_files(inputs, sim_dir, errors):
    """Every file named in __input_files_config.txt must exist in Data/ or
    Geometry/ under the simulation folder - except the combined BIB file,
    which ./run_all builds when needed."""
    sim_dir = Path(sim_dir)
    for key, (kind, _) in INPUTS_SPEC["data"].items():
        name = inputs.get(("data", key), "")
        if not name or key == "bib_combined":
            continue
        if not (sim_dir / "Data" / name).is_file():
            errors.append(f"__input_files_config.txt: [data] {key}: file not found: "
                          f"{_tilde(sim_dir / 'Data' / name)}")
    for key in INPUTS_SPEC["geometry"]:
        name = inputs.get(("geometry", key), "")
        if name and not (sim_dir / "Geometry" / name).is_file():
            errors.append(f"__input_files_config.txt: [geometry] {key}: file not found: "
                          f"{_tilde(sim_dir / 'Geometry' / name)}")


def describe_inputs(inputs, sim_dir):
    ipp = inputs.get(("data", "bib_ipp"), "")
    made_of = "plus + minus + ipp" if ipp else "plus + minus"
    combined = inputs[("data", "bib_combined")]
    if sim_dir and not (Path(sim_dir) / "Data" / combined).is_file():
        made_of += ", will be built"
    geo = [inputs[("geometry", k)] for k in INPUTS_SPEC["geometry"]]
    lines = []
    if sim_dir:
        lines.append(("folder", f"{_tilde(sim_dir)}  (Data/ and Geometry/)"))
    lines += [
        ("BIB plus", inputs[("data", "bib_plus")]),
        ("BIB minus", inputs[("data", "bib_minus")]),
        ("BIB ipp", ipp or "(none)"),
        ("BIB combined", f"{combined}  (= {made_of})"),
        ("signal", inputs[("data", "signal")]),
        ("geometry", ", ".join(geo[:2]) + ","),
        ("", ", ".join(geo[2:])),
    ]
    brief = (f"BIB data:    {combined} ({made_of.split(',')[0]})\n"
             f"Signal:      {inputs[('data', 'signal')]}")
    return lines, brief


def main(argv):
    brief = "--brief" in argv
    quiet = "--quiet" in argv
    sim_dir = None
    args = []
    it = iter(argv)
    for a in it:
        if a == "--sim-dir":
            sim_dir = next(it, None)
        elif a not in ("--brief", "--quiet"):
            args.append(a)
    if len(args) not in (2, 3):
        print(__doc__.split("Usage:")[1].split("Exit code")[0].strip())
        return 2
    cuts_path, smear_path = args[0], args[1]
    inputs_path = args[2] if len(args) == 3 else None

    cuts, errors = ct.read(cuts_path)
    smear = read_and_check(smear_path, SMEAR_SPEC, errors)
    inputs = read_and_check(inputs_path, INPUTS_SPEC, errors) if inputs_path else None
    if inputs is not None and sim_dir:
        check_input_files(inputs, sim_dir, errors)
    if errors:
        print("Problems found in the settings files:")
        for e in errors:
            print(f"  - {e}")
        return 1

    if quiet:
        return 0
    cut_lines, smear_lines = describe(cuts, smear)
    if brief:
        if inputs is not None:
            print(describe_inputs(inputs, sim_dir)[1])
        for line in brief_lines(cuts, smear):
            print(line)
        return 0

    if inputs is not None:
        print(f"Input files  ({inputs_path.split('/')[-1]})")
        for n, t_ in describe_inputs(inputs, sim_dir)[0]:
            print(f"  {n:<12} {t_}")

    print(f"Cuts         ({cuts_path.split('/')[-1]})")
    for n, t in cut_lines:
        print(f"  {n:<12} {t}")
    print(f"Resolutions  ({smear_path.split('/')[-1]})")
    for n, t in smear_lines:
        print(f"  {n:<12} {t}")
    for w in warnings_for(cuts, smear):
        print(f"  WARNING: {w}")
    for w in notes_for(cuts):
        print(f"  NOTE: {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
