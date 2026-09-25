"""
Read and check __cuts_config.txt: the selection cuts, one row per detector
subsystem, plus the display range of the zoomed N-1 plots and the
track-finding settings.

    [cuts]
    #             time   z0     pT
    #             (ns)   (mm)   (GeV/c)
    vxd_barrel  = 0.3    15     5
    ...
    ot_endcap   = 0.5    40     off

    [zoom]
    time = 2.0
    z0   = 100
    pT   = 1.0

    [track]
    min_hits_found = 5
    exclude_vertex_hits = true

A hit in subsystem S is kept if |t_corrected| <= time(S), |z0| <= z0(S)
and pT >= pT(S). "off" switches that one cut off in that one subsystem.
(pT comes from the curvature measured at the hit, so that cut is applied
as |1/R| <= 0.3*B/pT(S) - see bib_common.apply_cuts.)

Files in the format used before the table - one section per cut,
[t_corrected_ns], [z_axis_intercept_mm] and [momentum_gev], each with a
halfwidth, an enabled flag and an optional zoom_halfwidth - are still
read (the settings copies in older runs/ folders are like that): each
cut then applies to all six subsystems.

bib_common.load_cuts / load_track_params, check_configs.py and
make_highlights_pdf.py all read the file through read() below. Only the
standard library is used, so the settings check can run before anything
heavy is loaded.
"""
import configparser
import difflib
import math
from collections import Counter

# rows of the [cuts] table -> system ids, as in bib_common.SYSTEM_NAMES
ROWS = {"vxd_barrel": 1, "vxd_endcap": 2, "it_barrel": 3,
        "it_endcap": 4, "ot_barrel": 5, "ot_endcap": 6}
SYSTEM_IDS = tuple(ROWS.values())
NAMES = {1: "VXD barrel", 2: "VXD endcap", 3: "IT barrel",
         4: "IT endcap", 5: "OT barrel", 6: "OT endcap"}   # = bib_common.SYSTEM_NAMES

# columns of the table, in order. Keys are lower case because configparser
# lower-cases setting names ("pT" in [zoom] arrives as "pt").
COLUMNS = ("time", "z0", "pt")
TITLES = {"time": "time", "z0": "z0", "pt": "pT"}
UNITS = {"time": "ns", "z0": "mm", "pt": "GeV/c"}
CUT_NAMES = {"time": "t_corrected_ns", "z0": "z_axis_intercept_mm", "pt": "momentum_gev"}
# half-range of each variable's zoomed N-1 plot when [zoom] doesn't give one
ZOOM_DEFAULTS = {"time": 2.0, "z0": 100.0, "pt": 1.0}

OLD_SECTIONS = {"t_corrected_ns": "time", "z_axis_intercept_mm": "z0", "momentum_gev": "pt"}


def _suggest(name, candidates, shown=None):
    close = difflib.get_close_matches(name, list(candidates), n=1, cutoff=0.6)
    if not close:
        return ""
    return f" - did you mean '{(shown or {}).get(close[0], close[0])}'?"


def _number(raw):
    """float(raw) if raw is a finite number, else None."""
    try:
        x = float(raw)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def read(path):
    """
    Parse and check a cuts file. Returns (settings, problems): `problems`
    is a list of messages, empty if the file is fine; `settings` is None
    if there are problems, otherwise
        {"cuts":  {"time": {system id: value, or None if off}, "z0": {...}, "pt": {...}},
         "zoom":  {"time": x, "z0": x, "pt": x},
         "track": {"min_hits_found": n, "exclude_vertex_hits": True/False},
         "old_format": True/False}
    """
    label = str(path).split("/")[-1]
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    try:
        with open(path) as f:
            cp.read_file(f, source=label)
    except FileNotFoundError:
        return None, [f"{label}: file not found: {path}"]
    except configparser.DuplicateOptionError as e:
        return None, [f"{label}: [{e.section}] {e.option} appears twice (line {e.lineno})"]
    except configparser.DuplicateSectionError as e:
        return None, [f"{label}: section [{e.section}] appears twice (line {e.lineno})"]
    except configparser.Error as e:
        return None, [f"{label}: {' '.join(str(e).split())}"]

    problems = []

    def problem(msg):
        problems.append(f"{label}: {msg}")

    if cp.has_section("cuts") or not any(s in OLD_SECTIONS for s in cp.sections()):
        settings = _read_table(cp, problem)
    else:
        settings = _read_old(cp, problem)
    settings["track"] = _read_track(cp, problem)
    return (None if problems else settings), problems


def load(path):
    """read(), raising ValueError with the list of problems if there are any."""
    settings, problems = read(path)
    if problems:
        raise ValueError("problems in the cuts settings file:\n  " + "\n  ".join(problems))
    return settings


def _read_table(cp, problem):
    for sec in cp.sections():
        if sec in OLD_SECTIONS:
            col = OLD_SECTIONS[sec]
            problem(f"[{sec}] is the old way of setting the {TITLES[col]} cut - it now goes in "
                    f"the {TITLES[col]} column of the [cuts] table, so delete [{sec}]")
        elif sec not in ("cuts", "zoom", "track"):
            problem(f"unknown section [{sec}]{_suggest(sec, ('cuts', 'zoom', 'track'))}")

    cuts = {c: {} for c in COLUMNS}
    if not cp.has_section("cuts"):
        problem("missing section [cuts] (the table of cuts, one row per subsystem)")
    else:
        explained = set()     # rows already reported as misspelled or indented
        for row, raw in cp["cuts"].items():
            if row not in ROWS:
                close = difflib.get_close_matches(row, list(ROWS), n=1, cutoff=0.6)
                explained.update(close)
                problem(f"[cuts] has an unknown row '{row}'{_suggest(row, ROWS)} "
                        f"(the rows are {', '.join(ROWS)})")
                continue
            if "\n" in raw:
                indented = [ln.split("=")[0].strip().lower() for ln in raw.split("\n")[1:]]
                explained.update(indented)
                problem(f"[cuts] {row}: the line after it starts with spaces, so it is read as "
                        f"part of the {row} row - remove the spaces at the start of that line")
                continue
            cells = raw.replace(",", " ").split()
            if len(cells) != len(COLUMNS):
                problem(f"[cuts] {row} = {raw.strip()!r}: expected 3 values - time, z0, pT "
                        f"(each a number, or off) - but found {len(cells)}")
                continue
            for col, cell in zip(COLUMNS, cells):
                if cell.lower() == "off":
                    cuts[col][ROWS[row]] = None
                    continue
                x = _number(cell)
                if x is None:
                    problem(f"[cuts] {row}: {TITLES[col]} = {cell!r} is not a number (or off)")
                elif x < 0:
                    problem(f"[cuts] {row}: {TITLES[col]} = {cell} must not be negative")
                else:
                    cuts[col][ROWS[row]] = x
        missing = [r for r in ROWS if r not in cp["cuts"] and r not in explained]
        if missing:
            problem(f"[cuts] has no row for {', '.join(missing)} - every subsystem needs "
                    f"one (write off for a cut that should not apply there)")

    zoom = dict(ZOOM_DEFAULTS)
    if cp.has_section("zoom"):
        for key, raw in cp["zoom"].items():
            if key not in ZOOM_DEFAULTS:
                problem(f"unknown setting '{key}' in [zoom]{_suggest(key, ZOOM_DEFAULTS, TITLES)} "
                        f"(the settings are time, z0 and pT)")
                continue
            x = _number(raw)
            if x is None or x <= 0:
                problem(f"[zoom] {TITLES[key]} = {raw!r} must be a number above zero")
            else:
                zoom[key] = x
    return {"cuts": cuts, "zoom": zoom, "old_format": False}


def _read_old(cp, problem):
    """The format used before the table: one section per cut, same cut everywhere."""
    known = list(OLD_SECTIONS) + ["track"]
    for sec in cp.sections():
        if sec not in known:
            problem(f"unknown section [{sec}]{_suggest(sec, known)}")
    cuts = {c: {s: None for s in SYSTEM_IDS} for c in COLUMNS}
    zoom = dict(ZOOM_DEFAULTS)
    for sec, col in OLD_SECTIONS.items():
        if not cp.has_section(sec):
            problem(f"missing section [{sec}]")
            continue
        s = cp[sec]
        for key in s:
            if key not in ("halfwidth", "enabled", "zoom_halfwidth"):
                problem(f"unknown setting '{key}' in [{sec}]"
                        f"{_suggest(key, ('halfwidth', 'enabled', 'zoom_halfwidth'))}")
        try:
            enabled = s.getboolean("enabled")
            if enabled is None:
                problem(f"missing setting 'enabled' in [{sec}]")
        except ValueError:
            problem(f"[{sec}] enabled = {s['enabled']!r} is not true or false")
            enabled = None
        x = _number(s.get("halfwidth"))
        if "halfwidth" not in s:
            problem(f"missing setting 'halfwidth' in [{sec}]")
        elif x is None or x < 0:
            problem(f"[{sec}] halfwidth = {s['halfwidth']!r} must be a number, not negative")
        elif enabled:
            cuts[col] = {i: x for i in SYSTEM_IDS}
        if "zoom_halfwidth" in s:
            z = _number(s["zoom_halfwidth"])
            if z is None or z <= 0:
                problem(f"[{sec}] zoom_halfwidth = {s['zoom_halfwidth']!r} must be a number "
                        f"above zero")
            else:
                zoom[col] = z
    return {"cuts": cuts, "zoom": zoom, "old_format": True}


def _read_track(cp, problem):
    track = {"min_hits_found": 5, "exclude_vertex_hits": False}
    if not cp.has_section("track"):
        problem("missing section [track] (min_hits_found and exclude_vertex_hits)")
        return track
    s = cp["track"]
    for key in s:
        if key not in track:
            problem(f"unknown setting '{key}' in [track]{_suggest(key, track)}")
    if "min_hits_found" not in s:
        problem("missing setting 'min_hits_found' in [track]")
    else:
        try:
            n = s.getint("min_hits_found")
            if n < 0:
                raise ValueError
            track["min_hits_found"] = n
        except ValueError:
            problem(f"[track] min_hits_found = {s['min_hits_found']!r} must be a whole number, "
                    f"not negative")
    if "exclude_vertex_hits" not in s:
        problem("missing setting 'exclude_vertex_hits' in [track]")
    else:
        try:
            track["exclude_vertex_hits"] = s.getboolean("exclude_vertex_hits")
        except ValueError:
            problem(f"[track] exclude_vertex_hits = {s['exclude_vertex_hits']!r} is not true "
                    f"or false")
    return track


# ---------------------------------------------------------------- describing
def fmt(x):
    """A table cell: the value, or 'off'."""
    return "off" if x is None else f"{x:g}"


def condition(col, x):
    """One cut as a condition on the hit, e.g. '|z0| <= 15 mm'."""
    if col == "time":
        return f"|t_corrected| <= {x:g} ns"
    if col == "z0":
        return f"|z0| <= {x:g} mm"
    return f"pT >= {x:g} GeV/c"


def values(settings, col):
    """The six subsystems' values of one cut (None = off), in subsystem order."""
    return [settings["cuts"][col][s] for s in SYSTEM_IDS]


def is_uniform(settings, col=None):
    """True if the cut (or, with col=None, every cut) is the same in all six subsystems."""
    return all(len(set(values(settings, c))) == 1 for c in ([col] if col else COLUMNS))


def _cuts_something(col, v):
    """False if a value doesn't cut anything: off, or pT >= 0."""
    return v is not None and not (col == "pt" and v == 0)


def is_on(settings, col):
    """True if the cut removes hits in at least one subsystem."""
    return any(_cuts_something(col, v) for v in values(settings, col))


def summary(settings, col):
    """
    One cut in words, compactly: the value most subsystems use, then the
    exceptions - e.g. '|z0| <= 15 mm (except IT endcap 30, OT endcap off)',
    or 'no pT cut' when it is off everywhere.
    """
    vals = values(settings, col)
    common = Counter(vals).most_common(1)[0][0]      # ties: the first one listed
    text = condition(col, common) if common is not None else f"no {TITLES[col]} cut"
    others = [f"{NAMES[s]} {fmt(v)}" for s, v in zip(SYSTEM_IDS, vals) if v != common]
    return text + (f" (except {', '.join(others)})" if others else "")


def value_range(settings, col, unit=True):
    """The span of a cut's values over the subsystems where it is on, e.g.
    '15 to 40 mm' or '15 mm' (None if it is off everywhere)."""
    on = sorted({v for v in values(settings, col) if _cuts_something(col, v)})
    if not on:
        return None
    text = f"{on[0]:g}" if len(on) == 1 else f"{on[0]:g} to {on[-1]:g}"
    return f"{text} {UNITS[col]}" if unit else text


def off_in(settings, col):
    """Names of the subsystems where a cut is off (or, for pT, set to 0)."""
    return [NAMES[s] for s, v in zip(SYSTEM_IDS, values(settings, col))
            if not _cuts_something(col, v)]


def table_rows(settings):
    """Header and rows of the cut table, as strings: subsystem, time, z0, pT."""
    header = [f"{TITLES[c]} ({UNITS[c]})" for c in COLUMNS]
    rows = [[NAMES[s]] + [fmt(settings["cuts"][c][s]) for c in COLUMNS] for s in SYSTEM_IDS]
    return header, rows
