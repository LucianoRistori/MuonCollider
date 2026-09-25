"""
Check cuts_config.txt and smearing_config.txt before a run, and print a
short, readable summary of the settings they contain.

Why this exists: Python's configparser (used by bib_common's loaders)
silently falls back to a default whenever a section or setting it asks
for isn't there. So a typo - "sigma_u_m" instead of "sigma_u_mm", or
"[angle-u]" instead of "[angle_u]" - raises no error: that value just
quietly stays at its default (usually "no smearing", or "no cut"), and a
whole run is produced with different settings than intended. This script
rejects, before anything is run:
  - unknown or misspelled sections/settings (suggesting the closest
    valid name), and missing ones,
  - values that aren't a number / true-false,
  - negative widths or sigmas,
  - duplicated sections or settings,
and warns about settings that are probably not what was meant (e.g. a
sigma set but its section left disabled).

The accepted sections/settings mirror what bib_common.load_cuts,
load_track_params and load_smearing_config read - keep them in sync if
a setting is added there.

Usage:
    python3 check_configs.py <cuts_config.txt> <smearing_config.txt> [--brief]
Exit code 0 = OK (warnings allowed), 1 = at least one error.
--brief prints just two summary lines (used for runs/<id>/summary.txt).
"""
import configparser
import difflib
import math
import sys

NUMBER, WHOLE, TRUEFALSE = "number", "whole number", "true/false value"

# section -> {setting: (kind, required)}
CUTS_SPEC = {
    "t_corrected_ns": {"halfwidth": (NUMBER, True), "enabled": (TRUEFALSE, True),
                       "zoom_halfwidth": (NUMBER, False), "center": (NUMBER, False)},
    "z_axis_intercept_mm": {"halfwidth": (NUMBER, True), "enabled": (TRUEFALSE, True),
                            "zoom_halfwidth": (NUMBER, False), "center": (NUMBER, False)},
    "momentum_gev": {"halfwidth": (NUMBER, True), "enabled": (TRUEFALSE, True),
                     "zoom_halfwidth": (NUMBER, False)},
    "track": {"min_hits_found": (WHOLE, True), "exclude_vertex_hits": (TRUEFALSE, True)},
}
SMEAR_SPEC = {
    "position": {"sigma_u_mm": (NUMBER, True), "sigma_v_mm": (NUMBER, True),
                 "enabled": (TRUEFALSE, True)},
    "time": {"sigma_t_ns": (NUMBER, True), "enabled": (TRUEFALSE, True)},
    "angle_u": {"sigma_deg": (NUMBER, True), "enabled": (TRUEFALSE, True)},
    "angle_v": {"sigma_deg": (NUMBER, True), "enabled": (TRUEFALSE, True)},
    "general": {"seed": (WHOLE, True)},
}
MUST_NOT_BE_NEGATIVE = {"halfwidth", "zoom_halfwidth", "sigma_u_mm", "sigma_v_mm",
                        "sigma_t_ns", "sigma_deg", "min_hits_found"}
# same defaults as bib_common.ZOOM_HALFWIDTH_DEFAULTS
ZOOM_DEFAULTS = {"t_corrected_ns": 2.0, "z_axis_intercept_mm": 100.0, "momentum_gev": 1.0}


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


def describe(cuts, smear):
    """Human-readable lines: list of (name, text) for cuts, then smearing."""
    def cut_line(sec, var, unit):
        if not cuts.get((sec, "enabled"), False):
            return "off"
        hw = cuts.get((sec, "halfwidth"))
        c = cuts.get((sec, "center"), 0.0)
        inner = var if c == 0 else f"{var} - ({_g(c)})"
        return f"|{inner}| <= {_g(hw)} {unit}"

    mom = (f"pT >= {_g(cuts.get(('momentum_gev', 'halfwidth')))} GeV/c"
           if cuts.get(("momentum_gev", "enabled"), False) else "off")
    trk = f">= {cuts.get(('track', 'min_hits_found'))} surviving hits"
    trk += (", vertex-detector hits not counted"
            if cuts.get(("track", "exclude_vertex_hits"), False) else ", all subsystems counted")
    cut_lines = [
        ("time", cut_line("t_corrected_ns", "t_corrected", "ns")),
        ("z-intercept", cut_line("z_axis_intercept_mm", "z0", "mm")),
        ("momentum", mom),
        ("track found", trk),
    ]

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


def warnings_for(cuts, smear):
    warn = []
    for sec, var, unit in (("t_corrected_ns", "time", "ns"), ("z_axis_intercept_mm", "z-intercept", "mm")):
        if cuts.get((sec, "enabled"), False):
            hw = cuts.get((sec, "halfwidth"), 0.0)
            zoom = cuts.get((sec, "zoom_halfwidth"), ZOOM_DEFAULTS[sec])
            if hw == 0:
                warn.append(f"cuts_config.txt: [{sec}] halfwidth = 0 rejects essentially every hit")
            elif hw > zoom:
                warn.append(f"cuts_config.txt: {var} cut at +/-{_g(hw)} {unit} lies outside the zoomed "
                            f"plots' range (+/-{_g(zoom)} {unit}) - raise zoom_halfwidth in "
                            f"[{sec}] to see the cut lines there")
    if cuts.get(("momentum_gev", "enabled"), False):
        hw = cuts.get(("momentum_gev", "halfwidth"), 0.0)
        zoom = cuts.get(("momentum_gev", "zoom_halfwidth"), ZOOM_DEFAULTS["momentum_gev"])
        if hw < zoom:
            warn.append(f"cuts_config.txt: momentum cut at {_g(hw)} GeV/c is below the zoomed plots' "
                        f"lower edge ({_g(zoom)} GeV/c) - lower zoom_halfwidth in [momentum_gev] "
                        f"to see the cut lines there")
    for sec, keys in (("position", ("sigma_u_mm", "sigma_v_mm")), ("time", ("sigma_t_ns",)),
                      ("angle_u", ("sigma_deg",)), ("angle_v", ("sigma_deg",))):
        if (sec, "enabled") in smear and not smear[(sec, "enabled")]:
            if any(smear.get((sec, k), 0.0) > 0 for k in keys):
                warn.append(f"smearing_config.txt: [{sec}] has a nonzero sigma but enabled = false, "
                            f"so it will NOT be applied")
    return warn


def main(argv):
    brief = "--brief" in argv
    paths = [a for a in argv if a != "--brief"]
    if len(paths) != 2:
        print(__doc__.split("Usage:")[1].split("Exit code")[0].strip())
        return 2
    cuts_path, smear_path = paths

    errors = []
    cuts = read_and_check(cuts_path, CUTS_SPEC, errors)
    smear = read_and_check(smear_path, SMEAR_SPEC, errors)
    if errors:
        print("Problems found in the settings files:")
        for e in errors:
            print(f"  - {e}")
        return 1

    cut_lines, smear_lines = describe(cuts, smear)
    if brief:
        print("Cuts:        " + "; ".join(f"{n} {t}" for n, t in cut_lines))
        print("Resolutions: " + "; ".join(f"{n} {t}" for n, t in smear_lines))
        return 0

    print(f"Cuts         ({cuts_path})")
    for n, t in cut_lines:
        print(f"  {n:<12} {t}")
    print(f"Resolutions  ({smear_path})")
    for n, t in smear_lines:
        print(f"  {n:<12} {t}")
    for w in warnings_for(cuts, smear):
        print(f"  WARNING: {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
