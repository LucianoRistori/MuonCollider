"""
Build a landscape (16:9) PDF presentation of a run's key plots: a title
page (input files, settings, BIB rejection factor and signal efficiency
for pT -> infinity, track-finding efficiency for pT -> infinity), a page
with the cuts and results per subsystem when the cuts differ between
subsystems, then one plot per page with a short caption underneath.
Captions take their numbers (cut values, resolutions, results) from the
run folder itself, so they are always right for that run.

./run_all makes this PDF at the end of every successful run. To make it
again for an earlier run:
Usage:
    python3 make_highlights_pdf.py <run folder>
writes <run folder>/_highlights/highlights_<run name>.pdf
"""
import csv
import re
import sys
import textwrap
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

import check_configs as cc
import cuts_table as ct
from bib_common import format_rejection_factor

PAGE = (13.333, 7.5)                     # inches, 16:9 landscape
INK, MUTED, RULE = "#1a1a1a", "#6b6b6b", "#c8c8c8"
B_FIELD_T = 5.0                          # as in incidence_angle_plots / add_time_of_flight
CUT_WORDS = {"time": "time cut", "z0": "$z_0$ cut", "pt": "$p_T$ cut"}


def _find(run_dir, *names):
    for n in names:
        if (run_dir / n).is_file():
            return run_dir / n
    return None


def load_run(run_dir):
    """Everything the pages need, read from the run folder."""
    cuts_file = _find(run_dir, "__cuts_config.txt", "cuts_config.txt")
    smear_file = _find(run_dir, "__smearing_config.txt", "smearing_config.txt")
    if cuts_file is None or smear_file is None:
        sys.exit(f"No settings files (__cuts_config.txt, __smearing_config.txt) in {run_dir}")
    cuts, errors = ct.read(str(cuts_file))
    smear = cc.read_and_check(str(smear_file), cc.SMEAR_SPEC, errors)
    if errors:
        sys.exit("Problems in the run's settings files:\n  " + "\n  ".join(errors))

    # input files: from the settings copy if the run has one, else from the log
    inputs = {"combined": "?", "made_of": "", "signal": "?"}
    p = _find(run_dir, "__input_files_config.txt", "input_files_config.txt")
    if p:
        v = cc.read_and_check(str(p), cc.INPUTS_SPEC, [])
        if ("data", "bib_combined") in v:
            inputs = {"combined": v[("data", "bib_combined")],
                      "made_of": "plus + minus + ipp" if v.get(("data", "bib_ipp")) else "plus + minus",
                      "signal": v[("data", "signal")]}
    else:
        log = (run_dir / "run_log.txt").read_text(errors="replace")
        for key, kind in (("combined", "BIB"), ("signal", "signal")):
            m = re.search(rf"Load(?:ing {kind} file:|ed {kind}) (\S+?)(?::\s|\s*$)", log, re.M)
            if m:
                inputs[key] = Path(m.group(1)).name

    def read_csv(path):
        with open(path) as f:
            return list(csv.DictReader(f))
    by_region = run_dir / "step4_cuts" / "summary_by_region.csv"
    limit = run_dir / "step4_track_efficiency" / "track_efficiency_limit.csv"
    if not (by_region.is_file() and limit.is_file()):
        sys.exit(f"Run {run_dir.name} was made before the results were quoted as a rejection "
                 f"factor and as efficiencies for pT -> infinity (it has no {by_region.name}); "
                 f"run ./run_all again to get a PDF.")
    # (region, BIB rejection factor as text, signal efficiency for pT -> inf, its uncertainty)
    table = [(r["system_name"],
              format_rejection_factor(int(r["bib_n_hits"]), int(r["bib_n_hits_after_cuts"])),
              100 * float(r["signal_efficiency_pt_inf"]),
              100 * float(r["signal_efficiency_pt_inf_unc"])) for r in read_csv(by_region)]
    lim = read_csv(limit)[0]
    track = {"eff": 100 * float(lim["efficiency_pt_inf"]),
             "unc": 100 * float(lim["efficiency_pt_inf_unc"]),
             "pt_min": float(lim["fit_pt_min_gev"]), "n_tracks": int(lim["n_tracks"])}
    signal_pt_min = sorted({float(r["fit_pt_min_gev"]) for r in read_csv(by_region)})

    return {"name": run_dir.name, "cuts": cuts, "smear": smear, "inputs": inputs,
            "table": table, "track": track, "signal_pt_min": signal_pt_min}


def nice(s):
    return s.replace("<=", "≤").replace(">=", "≥")


def wrap(text, width=140):
    """textwrap.fill, but counting only the characters that show on the
    page: mathtext markup ($, _{...}, ^) takes no room."""
    lines, line, n = [], [], 0
    for word in text.split():
        w = len(re.sub(r"[${}_^\\]", "", word))
        if line and n + 1 + w > width:
            lines.append(" ".join(line))
            line, n = [], 0
        n += w + (1 if line else 0)
        line.append(word)
    if line:
        lines.append(" ".join(line))
    return "\n".join(lines)


def run_date(name):
    try:
        return datetime.strptime(name[:17], "%Y-%m-%d_%H%M%S").strftime("%A %d %B %Y, %H:%M")
    except ValueError:
        return ""


def captions(run):
    c, s = run["cuts"], run["smear"]
    uniform = ct.is_uniform(c)
    on = lambda col: ct.is_on(c, col)
    zoom = c["zoom"]
    names = {"time": "time", "z0": "$z_0$", "pt": "momentum"}

    def others(skip):
        o = [names[k] for k in ct.COLUMNS if k != skip and on(k)]
        if not o:
            return "N-1: no other cut is on, so all hits are included."
        whose = "the" if uniform else "their own subsystem's"
        cuts = "cut" if len(o) == 1 else "cuts"
        return f"N-1: only hits passing {whose} {' and '.join(o)} {cuts} are included."

    def lines_note(col):
        if not on(col):
            return f"The {CUT_WORDS[col]} is switched off in this run."
        if ct.is_uniform(c, col):
            v = ct.values(c, col)[0]
            at = {"time": f"at ±{v:g} ns", "z0": f"at ±{v:g} mm", "pt": f"at {v:g} GeV/c"}[col]
            return f"Dashed red lines: the {CUT_WORDS[col]} {at}; hits between the lines are kept."
        off = ct.off_in(c, col)
        return (f"Dashed red lines: each subsystem's own {CUT_WORDS[col]} "
                f"({ct.value_range(c, col)}, given in the panel titles); hits between the lines "
                f"are kept." + (f" The cut is off in {', '.join(off)}." if off else ""))

    rej_all, eff_all = run["table"][-1][1], run["table"][-1][2]
    trk = run["track"]
    if uniform:
        v = {col: ct.values(c, col)[0] for col in ct.COLUMNS}
        cut_list = [f"|t - t$_{{TOF}}$| ≤ {v['time']:g} ns" if on("time") else None,
                    f"|$z_0$| ≤ {v['z0']:g} mm" if on("z0") else None,
                    f"$p_T$ ≥ {v['pt']:g} GeV/c" if on("pt") else None]
        cut_sentence = f"Cuts: {', '.join(x for x in cut_list if x) or 'none'}."
    else:
        cut_sentence = "Each subsystem has its own cuts (page 2)."
    made_of = f" ({run['inputs']['made_of']})" if run["inputs"]["made_of"] else ""
    trk_min = c["track"]["min_hits_found"]
    vtx = ", not counting vertex-detector hits" if c["track"]["exclude_vertex_hits"] else ""
    t_res = (f"the {s[('time', 'sigma_t_ns')]:g} ns time smearing"
             if s.get(("time", "enabled")) and s.get(("time", "sigma_t_ns"), 0) > 0 else "no time smearing")

    if not on("pt"):
        pt_tail = ""
    elif ct.is_uniform(c, "pt"):
        pt_tail = (f"; it falls off toward low $p_T$ because of the $p_T$ ≥ "
                   f"{ct.values(c, 'pt')[0]:g} GeV/c cut")
    else:
        pt_tail = (f"; it falls off toward low $p_T$ because of the subsystems' $p_T$ cuts "
                   f"({ct.value_range(c, 'pt')})")
    return [
        ("density_before_after_cuts.png",
         "BIB hit density before and after the cuts",
         f"Mean BIB hit density per subsystem (hits/mm$^2$, log scale) for one collision's worth of "
         f"background{made_of}, before the cuts (light blue) and after all cuts combined (dark blue). "
         f"{cut_sentence} Overall they reduce the BIB hits by a factor {rej_all} (rejection "
         f"factor), and keep {eff_all:.1f}% of the hits of a signal muon with $p_T$ → ∞."),
        ("track_efficiency_vs_pt.png",
         "Track-finding efficiency vs. transverse momentum",
         f"Fraction of signal muons counted as found – at least {trk_min} of their hits survive "
         f"all cuts{vtx} – vs. generated $p_T$ ({trk['n_tracks']:,} muons, one per event; the axis is "
         f"linear in 1/$p_T$ and cropped where the efficiency is below 1%). For $p_T$ → ∞ it is "
         f"{trk['eff']:.1f} ± {trk['unc']:.1f}% (red diamond), from a fit of ε$_∞$ + c/$p_T^2$ to the "
         f"muons above {trk['pt_min']:g} GeV/c (dashed){pt_tail}."),
        ("z_axis_intercept_per_subsystem_zoom_n1.png",
         "$z$-axis intercept $z_0$ (N-1)",
         f"Where each hit's direction, extrapolated in the $r$-$z$ plane, crosses the beam line, for "
         f"BIB (blue) and signal (green) in the six subsystems, zoomed to ±{zoom['z0']:g} mm, in hits "
         f"per collision per bin. {others('z0')} " + lines_note("z0")),
        ("inv_radius_per_subsystem_zoom_n1.png",
         "Transverse momentum from the hit curvature (N-1)",
         f"$p_T$ estimated from each hit alone: the circle through the beam line tangent to the hit's "
         f"direction (B = {B_FIELD_T:g} T). The axis is linear in the curvature 1/R, labelled in $p_T$ "
         f"(±∞ = straight track, at the center), for |$p_T$| ≥ {zoom['pt']:g} GeV/c; BIB (blue) vs. "
         f"signal (green), per collision. {others('pt')} " + lines_note("pt")),
        ("time_corrected_per_subsystem_zoom_n1.png",
         "Corrected hit time (N-1)",
         f"Hit time minus the time of flight expected for a particle coming from the interaction "
         f"point, for BIB (blue) and signal (green) in the six subsystems, zoomed to ±"
         f"{zoom['time']:g} ns, in hits per collision per bin; hit times include {t_res}. "
         f"{others('time')} " + lines_note("time")),
    ]


def page_frame(fig, title, run, page, n_pages, footer):
    if title:
        fig.text(0.04, 0.935, title, fontsize=19, weight="bold", color=INK, va="center")
        fig.text(0.96, 0.935, f"Run {run['name']}", fontsize=10, color=MUTED, ha="right", va="center")
        fig.add_artist(Line2D([0.04, 0.96], [0.892, 0.892], color=RULE, lw=0.8))
    fig.add_artist(Line2D([0.04, 0.96], [0.052, 0.052], color=RULE, lw=0.6))
    fig.text(0.04, 0.026, footer, fontsize=8, color=MUTED, va="center")
    fig.text(0.96, 0.026, f"{page} / {n_pages}", fontsize=8, color=MUTED, ha="right", va="center")


def draw_table(ax, header, cells, col_widths, fontsize, row_scale):
    """A table in the style of the title page: horizontal rules only, first
    column left-aligned, header and ALL row in bold. Header labels may
    have two lines ("a\nb"); the header row is then made taller."""
    ax.set_axis_off()
    tab = ax.table(cellText=cells, colLabels=header, colLoc="right", cellLoc="right",
                   loc="upper left", edges="horizontal", colWidths=col_widths)
    tab.auto_set_font_size(False)
    tab.set_fontsize(fontsize)
    tab.scale(1, row_scale)
    two_lines = any("\n" in h for h in header)
    for (r, col), cell in tab.get_celld().items():
        cell.set_edgecolor(RULE)
        if col == 0:
            cell.set_text_props(ha="left")
        if r == 0 or cells[r - 1][0] == "ALL":
            cell.set_text_props(weight="bold")
        if r == 0 and two_lines:
            cell.set_height(cell.get_height() * 1.75)
            cell.set_text_props(va="top" if col == 0 else "center")
    return tab


def main(argv):
    if len(argv) != 1:
        print(__doc__.split("Usage:")[1].strip())
        return 2
    run_dir = Path(argv[0]).resolve()
    hl = run_dir / "_highlights"
    run = load_run(run_dir)
    per_system = not ct.is_uniform(run["cuts"])
    brief = [nice(line) for line in cc.brief_lines(run["cuts"], run["smear"])]
    cuts_brief = ("cuts: set per subsystem (page 2)" if per_system
                  else brief[0].replace("Cuts:        ", "cuts: "))
    footer = (f"Muon Collider BIB study  ·  run {run['name']}  ·  {cuts_brief}  ·  "
              + brief[2].replace("Resolutions: ", "resolutions: "))
    pages = [p for p in captions(run) if (hl / p[0]).is_file()]
    missing = [p[0] for p in captions(run) if not (hl / p[0]).is_file()]
    n_pages = 1 + per_system + len(pages)
    out = hl / f"highlights_{run['name']}.pdf"

    with PdfPages(out) as pdf:
        # ---- title page
        fig = plt.figure(figsize=PAGE)
        fig.text(0.06, 0.84, "Muon Collider: beam-induced background suppression",
                 fontsize=26, weight="bold", color=INK)
        fig.text(0.06, 0.775, f"Highlights of run {run['name']}   ({run_date(run['name'])})",
                 fontsize=14, color=MUTED)
        fig.add_artist(Line2D([0.06, 0.94], [0.745, 0.745], color=RULE, lw=0.8))
        made_of = f"  ({run['inputs']['made_of']})" if run["inputs"]["made_of"] else ""
        settings = []
        for line in brief:
            label, value = [x.strip() for x in line.split(":", 1)]
            if label == "Cuts" and per_system:
                value = "set per subsystem - see page 2"
            value = value.replace("|z0|", "|$z_0$|").replace("pT", "$p_T$")
            if ", seed " in value:
                value, seed = value.rsplit(", seed ", 1)
                settings += [(label, value), ("Random seed", seed)]
            else:
                settings.append((label, value))
        blocks = [("Input files", [("BIB", run["inputs"]["combined"] + made_of),
                                   ("Signal", run["inputs"]["signal"])]),
                  ("Settings", settings)]
        y = 0.675
        for head, rows in blocks:
            fig.text(0.06, y, head, fontsize=14, weight="bold", color=INK)
            y -= 0.055
            for label, value in rows:
                fig.text(0.06, y, label, fontsize=11.5, color=MUTED)
                for part in textwrap.wrap(value, 62):
                    fig.text(0.165, y, part, fontsize=11.5, color=INK)
                    y -= 0.045
            y -= 0.04
        ax = fig.add_axes([0.6, 0.22, 0.34, 0.45])
        ax.set_title("BIB rejection and signal efficiency\n(all cuts combined)", fontsize=12,
                     weight="bold", color=INK, loc="left")
        tab = draw_table(ax, ["region", "BIB rejection\nfactor", "signal efficiency\n($p_T$ → ∞)"],
                         [[name, rej, f"{eff:.1f} ± {unc:.1f}%"] for name, rej, eff, unc in run["table"]],
                         [0.32, 0.3, 0.38], 11, 1.55)
        fig.canvas.draw()
        box = tab.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
        trk = run["track"]
        fig.text(0.6, box.y0 - 0.06, f"Track-finding efficiency for $p_T$ → ∞:  "
                 f"{trk['eff']:.1f} ± {trk['unc']:.1f}%", fontsize=11.5, weight="bold", color=INK)
        pt_min = run["signal_pt_min"]
        fig.text(0.6, box.y0 - 0.105,
                 "rejection factor = 1/(1 − fraction of BIB hits removed)\n"
                 "efficiencies for $p_T$ → ∞: fit of ε$_∞$ + c/$p_T^2$ to muons above "
                 + (f"{pt_min[0]:g} GeV/c" if len(pt_min) == 1 else "twice the $p_T$ cut"),
                 fontsize=8.5, color=MUTED, va="top", linespacing=1.5)
        page_frame(fig, "", run, 1, n_pages, footer)
        pdf.savefig(fig)
        plt.close(fig)

        # ---- cuts per subsystem, next to their results
        if per_system:
            fig = plt.figure(figsize=PAGE)
            page_frame(fig, "Cuts and results per subsystem", run, 2, n_pages, footer)
            cuts = run["cuts"]
            by_name = {ct.NAMES[s]: s for s in ct.SYSTEM_IDS}
            cells = []
            for name, rej, eff, unc in run["table"]:
                s = by_name.get(name)
                cut_cells = ([ct.fmt(cuts["cuts"][col][s]) for col in ct.COLUMNS] if s
                             else ["", "", ""])
                cells.append([name] + cut_cells + [rej, f"{eff:.1f} ± {unc:.1f}%"])
            ax = fig.add_axes([0.1, 0.27, 0.8, 0.56])
            draw_table(ax, ["region", "time\n(ns)", "$z_0$\n(mm)", "$p_T$\n(GeV/c)",
                            "BIB rejection\nfactor", "signal efficiency\n($p_T$ → ∞)"],
                       cells, [0.2, 0.12, 0.12, 0.14, 0.19, 0.23], 14, 2.1)
            caption = ("A hit is kept only if it passes all three cuts of its own subsystem: "
                       "|t - t$_{TOF}$| ≤ time, |$z_0$| ≤ $z_0$ and $p_T$ ≥ $p_T$ (off: that cut "
                       "is not applied there). BIB rejection factor: BIB hits before / after all "
                       "cuts, = 1/(1 − R) with R the fraction removed. Signal efficiency: fraction "
                       "of a signal muon's hits kept, in the limit $p_T$ → ∞ (fit of ε$_∞$ + "
                       "c/$p_T^2$ to muons well above the $p_T$ cut).")
            fig.text(0.04, 0.172, wrap(caption), fontsize=11.5, color=INK,
                     va="top", linespacing=1.4)
            pdf.savefig(fig)
            plt.close(fig)

        # ---- one plot per page
        for i, (png, title, caption) in enumerate(pages, start=2 + per_system):
            fig = plt.figure(figsize=PAGE)
            page_frame(fig, title, run, i, n_pages, footer)
            ax = fig.add_axes([0.03, 0.195, 0.94, 0.685])
            ax.imshow(plt.imread(hl / png), interpolation="none")
            ax.set_axis_off()
            fig.text(0.04, 0.172, wrap(caption), fontsize=11.5, color=INK,
                     va="top", linespacing=1.4)
            pdf.savefig(fig)
            plt.close(fig)

        info = pdf.infodict()
        info["Title"] = f"Muon Collider BIB study - highlights of run {run['name']}"
        info["Creator"] = "make_highlights_pdf.py"

    print(f"Wrote {out.name} ({n_pages} pages) to _highlights/")
    if missing:
        print(f"  (not found, so not included: {', '.join(missing)})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
