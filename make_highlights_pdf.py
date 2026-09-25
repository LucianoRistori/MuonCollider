"""
Build a landscape (16:9) PDF presentation of a run's key plots: a title
page (input files, settings, BIB rejection vs. signal efficiency), then
one plot per page with a short caption underneath. Captions take their
numbers (cut values, resolutions, results) from the run folder itself,
so they are always right for that run.

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

PAGE = (13.333, 7.5)                     # inches, 16:9 landscape
INK, MUTED, RULE = "#1a1a1a", "#6b6b6b", "#c8c8c8"
B_FIELD_T = 5.0                          # as in incidence_angle_plots / add_time_of_flight


def _find(run_dir, *names):
    for n in names:
        if (run_dir / n).is_file():
            return run_dir / n
    return None


def load_run(run_dir):
    """Everything the pages need, read from the run folder."""
    errors = []
    cuts = cc.read_and_check(str(_find(run_dir, "__cuts_config.txt", "cuts_config.txt")),
                             cc.CUTS_SPEC, errors)
    smear = cc.read_and_check(str(_find(run_dir, "__smearing_config.txt", "smearing_config.txt")),
                              cc.SMEAR_SPEC, errors)
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
    bib = {r["system_name"]: r for r in read_csv(run_dir / "step4_cuts" / "cutflow_bib.csv")}
    sig = {r["system_name"]: r for r in read_csv(run_dir / "step4_cuts" / "cutflow_signal.csv")}
    table = [(name, 100 * (1 - float(bib[name]["frac_after_all"])),
              100 * float(sig[name]["frac_after_all"])) for name in bib]

    eff = read_csv(run_dir / "step4_track_efficiency" / "track_efficiency_vs_pt.csv")
    high = [r for r in eff if float(r["pT_lo_gev"]) >= 10.0]
    n_high = sum(int(float(r["n_tracks_examined"])) for r in high)
    plateau = 100 * sum(int(float(r["n_tracks_found"])) for r in high) / n_high if n_high else float("nan")
    n_tracks = sum(int(float(r["n_tracks_examined"])) for r in eff)

    return {"name": run_dir.name, "cuts": cuts, "smear": smear, "inputs": inputs,
            "table": table, "plateau": plateau, "n_tracks": n_tracks}


def nice(s):
    return s.replace("<=", "≤").replace(">=", "≥")


def run_date(name):
    try:
        return datetime.strptime(name[:17], "%Y-%m-%d_%H%M%S").strftime("%A %d %B %Y, %H:%M")
    except ValueError:
        return ""


def captions(run):
    c, s = run["cuts"], run["smear"]
    on = lambda sec: c.get((sec, "enabled"), False)
    hw = lambda sec: c.get((sec, "halfwidth"))
    zoom = lambda sec: c.get((sec, "zoom_halfwidth"), cc.ZOOM_DEFAULTS[sec])
    names = {"t_corrected_ns": "time", "z_axis_intercept_mm": "$z_0$",
             "momentum_gev": "momentum"}

    def others(skip):
        o = [names[k] for k in ("t_corrected_ns", "z_axis_intercept_mm", "momentum_gev")
             if k != skip and on(k)]
        return (f"N-1: only hits passing the {' and '.join(o)} cuts are included."
                if o else "N-1: no other cut is enabled, so all hits are included.")

    rej_all, eff_all = run["table"][-1][1], run["table"][-1][2]
    cut_list = [f"|t - t$_{{TOF}}$| ≤ {hw('t_corrected_ns'):g} ns" if on("t_corrected_ns") else None,
                f"|$z_0$| ≤ {hw('z_axis_intercept_mm'):g} mm" if on("z_axis_intercept_mm") else None,
                f"$p_T$ ≥ {hw('momentum_gev'):g} GeV/c" if on("momentum_gev") else None]
    cut_text = ", ".join(x for x in cut_list if x) or "none"
    made_of = f" ({run['inputs']['made_of']})" if run["inputs"]["made_of"] else ""
    trk_min = c.get(("track", "min_hits_found"))
    vtx = ", not counting vertex-detector hits" if c.get(("track", "exclude_vertex_hits")) else ""
    t_res = (f"the {s[('time', 'sigma_t_ns')]:g} ns time smearing"
             if s.get(("time", "enabled")) and s.get(("time", "sigma_t_ns"), 0) > 0 else "no time smearing")

    def lines_note(sec, what):
        return (f"Dashed red lines: the {what}; hits between the lines are kept."
                if on(sec) else f"The {what.split(' at ')[0]} is switched off in this run.")

    pt_tail = (f", and drops to zero just below the $p_T$ ≥ {hw('momentum_gev'):g} GeV/c cut"
               if on("momentum_gev") else "")
    return [
        ("density_before_after_cuts.png",
         "BIB hit density before and after the cuts",
         f"Mean BIB hit density per subsystem (hits/mm$^2$, log scale) for one collision's worth of "
         f"background{made_of}, before the cuts (light blue) and after all cuts combined (dark blue). "
         f"Cuts: {cut_text}. Overall the cuts remove {rej_all:.2f}% of the BIB hits and keep "
         f"{eff_all:.2f}% of the signal hits."),
        ("track_efficiency_vs_pt.png",
         "Track-finding efficiency vs. transverse momentum",
         f"Fraction of signal muons counted as found \u2013 at least {trk_min} of their hits survive "
         f"all cuts{vtx} \u2013 vs. generated $p_T$ ({run['n_tracks']:,} muons, one per event; the axis is "
         f"linear in 1/$p_T$ and cropped where the efficiency is below 1%). Above 10 GeV/c the "
         f"efficiency averages {run['plateau']:.1f}%{pt_tail}."),
        ("z_axis_intercept_per_subsystem_zoom_n1.png",
         "$z$-axis intercept $z_0$ (N-1)",
         f"Where each hit's direction, extrapolated in the $r$-$z$ plane, crosses the beam line, for "
         f"BIB (blue) and signal (green) in the six subsystems, zoomed to ±"
         f"{zoom('z_axis_intercept_mm'):g} mm, in hits per collision per bin. "
         f"{others('z_axis_intercept_mm')} "
         + lines_note("z_axis_intercept_mm", f"$z_0$ cut at ±{hw('z_axis_intercept_mm'):g} mm")),
        ("inv_radius_per_subsystem_zoom_n1.png",
         "Transverse momentum from the hit curvature (N-1)",
         f"$p_T$ estimated from each hit alone: the circle through the beam line tangent to the hit's "
         f"direction (B = {B_FIELD_T:g} T). The axis is linear in the curvature 1/R, labelled in $p_T$ "
         f"(±∞ = straight track, at the center), for |$p_T$| ≥ "
         f"{zoom('momentum_gev'):g} GeV/c; BIB (blue) vs. signal (green), per collision. "
         f"{others('momentum_gev')} "
         + lines_note("momentum_gev", f"cut at $p_T$ ≥ {hw('momentum_gev'):g} GeV/c")),
        ("time_corrected_per_subsystem_zoom_n1.png",
         "Corrected hit time (N-1)",
         f"Hit time minus the time of flight expected for a particle coming from the interaction "
         f"point, for BIB (blue) and signal (green) in the six subsystems, zoomed to ±"
         f"{zoom('t_corrected_ns'):g} ns, in hits per collision per bin; hit times include {t_res}. "
         f"{others('t_corrected_ns')} "
         + lines_note("t_corrected_ns", f"time cut at ±{hw('t_corrected_ns'):g} ns")),
    ]


def page_frame(fig, title, run, page, n_pages, footer):
    if title:
        fig.text(0.04, 0.935, title, fontsize=19, weight="bold", color=INK, va="center")
        fig.text(0.96, 0.935, f"Run {run['name']}", fontsize=10, color=MUTED, ha="right", va="center")
        fig.add_artist(Line2D([0.04, 0.96], [0.892, 0.892], color=RULE, lw=0.8))
    fig.add_artist(Line2D([0.04, 0.96], [0.052, 0.052], color=RULE, lw=0.6))
    fig.text(0.04, 0.026, footer, fontsize=8, color=MUTED, va="center")
    fig.text(0.96, 0.026, f"{page} / {n_pages}", fontsize=8, color=MUTED, ha="right", va="center")


def main(argv):
    if len(argv) != 1:
        print(__doc__.split("Usage:")[1].strip())
        return 2
    run_dir = Path(argv[0]).resolve()
    hl = run_dir / "_highlights"
    run = load_run(run_dir)
    brief = [nice(line) for line in cc.brief_lines(run["cuts"], run["smear"])]
    footer = (f"Muon Collider BIB study  ·  run {run['name']}  ·  "
              + brief[0].replace("Cuts:        ", "cuts: ") + "  ·  "
              + brief[2].replace("Resolutions: ", "resolutions: "))
    pages = [p for p in captions(run) if (hl / p[0]).is_file()]
    missing = [p[0] for p in captions(run) if not (hl / p[0]).is_file()]
    n_pages = len(pages) + 1
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
        ax.set_axis_off()
        ax.set_title("BIB rejection vs. signal efficiency\n(all cuts combined)", fontsize=12,
                     weight="bold", color=INK, loc="left")
        cells = [[name, f"{rej:.2f}%", f"{eff:.2f}%"] for name, rej, eff in run["table"]]
        tab = ax.table(cellText=cells, colLabels=["region", "BIB rejection", "signal efficiency"],
                       colLoc="right", cellLoc="right", loc="upper left", edges="horizontal",
                       colWidths=[0.36, 0.3, 0.34])
        tab.auto_set_font_size(False)
        tab.set_fontsize(11)
        tab.scale(1, 1.55)
        for (r, col), cell in tab.get_celld().items():
            cell.set_edgecolor(RULE)
            if col == 0:
                cell.set_text_props(ha="left")
            if r == 0 or cells[r - 1][0] == "ALL":
                cell.set_text_props(weight="bold")
        page_frame(fig, "", run, 1, n_pages, footer)
        pdf.savefig(fig)
        plt.close(fig)

        # ---- one plot per page
        for i, (png, title, caption) in enumerate(pages, start=2):
            fig = plt.figure(figsize=PAGE)
            page_frame(fig, title, run, i, n_pages, footer)
            ax = fig.add_axes([0.03, 0.195, 0.94, 0.685])
            ax.imshow(plt.imread(hl / png), interpolation="none")
            ax.set_axis_off()
            fig.text(0.04, 0.172, textwrap.fill(caption, 146), fontsize=11.5, color=INK,
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
