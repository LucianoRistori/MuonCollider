"""
Figures of Sections 15-16 of the paper, from the scan results in data/
(see data/README.txt). Writes figures/section15_*.png and
figures/section16_*.png. The Section 13-14 figures are made by
fake_rate_framework/section13_density_resolution.py and
section14_real_tower.py, run from figures/.

    cd docs/paper && python3 make_figures.py
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
MODELS = (("efakes_exact_helix", "exact helix", "#1b9e77", "-"),
          ("efakes_conservative_quad", "conservative quad", "#d95f02", "--"),
          ("efakes_line_B0", "line (B=0)", "#7570b3", ":"))


def load(name):
    with open(HERE / "data" / name) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def efakes_plot(x, rows, xlabel, today, today_label, logx, out):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for key, label, c, ls in MODELS:
        ax.plot(x, [r[key] for r in rows], color=c, ls=ls, lw=2, marker="o", ms=4, label=label)
    ax.axhline(1.0, color="gray", lw=0.8)
    ax.axvline(today, color="black", lw=1.0, alpha=0.6)
    ax.text(0.03, 0.97, f"vertical line: {today_label}", transform=ax.transAxes, fontsize=8, va="top")
    ax.set_yscale("log")
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$E[\#\mathrm{fakes}]$ (6-layer IT+OT barrel tower)")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(HERE / "figures" / out, dpi=150)
    plt.close(fig)


def eff_plot(x, rows, xlabel, today, today_label, logx, out):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for key, label, c in (("eff_inf", "6-layer IT+OT barrel tower (as $E[\\#\\mathrm{fakes}]$)", "#1b9e77"),
                          ("eff_inf_tracker", "whole tracker (as in the analysis pipeline)", "#d95f02")):
        ax.errorbar(x, [100 * r[key] for r in rows], yerr=[100 * r[key + "_unc"] for r in rows],
                    color=c, lw=2, marker="o", ms=4, capsize=3, label=label)
    ax.axvline(today, color="black", lw=1.0, alpha=0.6)
    ax.text(0.03, 0.97, f"vertical line: {today_label}", transform=ax.transAxes, fontsize=8, va="top")
    if logx:
        ax.set_xscale("log")
    ax.set_ylim(95.5, 100.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"track-finding efficiency, $p_T\to\infty$ (%)")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(HERE / "figures" / out, dpi=150)
    plt.close(fig)


(HERE / "figures").mkdir(exist_ok=True)
t = load("section15_time_scan.csv")
xt = np.array([r["sigma_t_ns"] for r in t]) * 1000
lab_t = r"timing resolution $\sigma_t$ (ps), all subsystems"
efakes_plot(xt, t, lab_t, 30, "today's 30 ps", True, "section15_efakes_vs_sigma_t.png")
eff_plot(xt, t, lab_t, 30, "today's 30 ps", True, "section15_efficiency_vs_sigma_t.png")
a = load("section16_angle_scan.csv")
xa = np.array([r["sigma_angle_deg"] for r in a])
lab_a = r"pointing-angle resolution $\sigma_\mathrm{angle}$ (deg), all subsystems"
efakes_plot(xa, a, lab_a, 1.0, r"today's 1$^\circ$", False, "section16_efakes_vs_sigma_angle.png")
eff_plot(xa, a, lab_a, 1.0, r"today's 1$^\circ$", False, "section16_efficiency_vs_sigma_angle.png")
print("wrote figures/section15_*.png, figures/section16_*.png")
