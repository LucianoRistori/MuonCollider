"""
Step 3 of rebuilding docs/efficiency_pt_infinity.pdf (see README.txt): the two
figures, fig_tracks.pdf and fig_hits.pdf, from fitdata.npz and numbers.json.
Uses matplotlib's pgf backend with pdflatex, so the figures are set in the
same fonts as the note (Palatino, mathpazo).
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator, AutoMinorLocator, ScalarFormatter


class CapFormatter(ScalarFormatter):
    """Plain tick labels, none above 100 %."""
    def __call__(self, x, pos=None):
        return "" if x > 100.0001 else super().__call__(x, pos)

plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.rcfonts": False,
    "pgf.preamble": r"\usepackage[T1]{fontenc}\usepackage{mathpazo}",
    "font.family": "serif",
    "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 9, "legend.fontsize": 7.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": False, "ytick.right": True,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.4,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
INK, BLUE, BAND, SHADE, MUTED = "#0b0b0b", "#2a78d6", "#cfe0f6", "#eeeeee", "#52514e"
HERE = Path(__file__).resolve().parent
d = np.load(HERE / "fitdata.npz")
R = json.load(open(HERE / "numbers.json"))
NAMES = {1: "VXD barrel", 2: "VXD endcap", 3: "IT barrel", 4: "IT endcap", 5: "OT barrel", 6: "OT endcap"}
EDGES = [10, 15, 20, 30, 50, 100, 1e9]

def binned(pt, y, edges, groups=None):
    pt = np.asarray(pt, float); y = np.asarray(y, float); out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (pt > lo) & (pt <= hi); n = s.sum(); e = y[s].mean()
        if groups is None:
            err = np.sqrt(e * (1 - e) / n)
        else:
            g = np.unique(groups[s], return_inverse=True)[1]
            err = np.sqrt(np.sum((np.bincount(g, weights=y[s]) - e * np.bincount(g)) ** 2)) / n
        x = 1 / pt[s] ** 2
        out.append((x.mean(), 1 / hi ** 2, 1 / lo ** 2, e, err))
    return out

def top_pt_axis(ax, ticks=(10, 12.5, 15, 20, 30), label=True):
    sec = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
    pos = [1 / p ** 2 for p in ticks] + [0.0]
    sec.xaxis.set_major_locator(FixedLocator(pos))
    sec.xaxis.set_major_formatter(FuncFormatter(lambda v, _: r"$\infty$" if v == 0 else f"{1/np.sqrt(v):g}"))
    sec.tick_params(direction="in", length=2.5, width=0.6, labelsize=7)
    if label:
        sec.set_xlabel(r"generated $p_T$ (GeV/$c$)", fontsize=7.5, labelpad=3)
    return sec

def fit_panel(ax, pt, y, fit, groups=None, ylim=None, text_loc="lower right", show_legend=False, title=None,
              scale=100.0):
    xs = np.linspace(0, 0.0102, 200)
    cov = np.array(fit["cov"])
    line = fit["a"] + fit["c"] * xs
    sig = np.sqrt(cov[0, 0] + 2 * xs * cov[0, 1] + xs ** 2 * cov[1, 1])
    ax.fill_between(xs, (line - sig) * scale, (line + sig) * scale, color=BAND, lw=0, zorder=1,
                    label=r"fit $\pm1\sigma$")
    ax.plot(xs, line * scale, color=BLUE, lw=1.3, zorder=2, label=r"fit $\varepsilon_\infty + c/p_T^2$")
    b = binned(pt, y, EDGES, groups)
    xm = np.array([r[0] for r in b]); xlo = np.array([r[1] for r in b]); xhi = np.array([r[2] for r in b])
    e = np.array([r[3] for r in b]) * scale; err = np.array([r[4] for r in b]) * scale
    ax.errorbar(xm, e, yerr=err, xerr=[xm - xlo, xhi - xm], fmt="o", ms=3.2, mfc=INK, mec=INK, ecolor=INK,
                elinewidth=0.7, capsize=0, zorder=4, label="data, in bins of $p_T$")
    ax.errorbar([0], [fit["a"] * scale], yerr=[fit["sa"] * scale], fmt="D", ms=4.2, mfc="white", mec=BLUE,
                mew=1.2, ecolor=BLUE, elinewidth=1.2, capsize=2.5, zorder=5, clip_on=False,
                label=r"$\varepsilon_\infty$ (intercept)")
    ax.set_xlim(-0.0004, 0.0102)
    ax.xaxis.set_major_locator(MultipleLocator(0.002))
    ax.xaxis.set_minor_locator(MultipleLocator(0.001))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "0" if v == 0 else f"{v*1000:g}"))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if ylim is not None:
        ax.set_ylim(*ylim)
    # no tick labels above 100 %
    ax.yaxis.set_major_formatter(CapFormatter(useOffset=False))
    txt = (r"$\varepsilon_\infty = %.2f \pm %.2f\,\%%$" % (fit["a"] * 100, fit["sa"] * 100) + "\n"
           + r"$c = %+.2f \pm %.2f$ (GeV/$c$)$^2$" % (fit["c"], fit["sc"]))
    ha = "right" if "right" in text_loc else "left"
    va = "bottom" if "lower" in text_loc else "top"
    ax.text(0.97 if ha == "right" else 0.04, 0.05 if va == "bottom" else 0.95, txt, transform=ax.transAxes,
            ha=ha, va=va, fontsize=7.5, color=INK,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
    if title:
        ax.set_title(title, pad=17)
    if show_legend:
        h, l = ax.get_legend_handles_labels()
        order = [2, 1, 0, 3]
        ax.legend([h[i] for i in order], [l[i] for i in order], loc="upper right", frameon=True,
                  framealpha=0.9, edgecolor="none", borderpad=0.3, handlelength=1.6)

# ---------------- Figure 1: tracks ----------------
tpt, found = d["trk_pt"].astype(float), d["trk_found"].astype(float)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 2.75), gridspec_kw=dict(wspace=0.28))
# (a) efficiency vs 1/pT over the whole sample: the gun is flat in 1/pT, so every bin has the same statistics
u = 1 / tpt
edges = np.linspace(0, 2 / 3, 41)
idx = np.digitize(u, edges) - 1
n = np.bincount(idx, minlength=40)[:40]; k = np.bincount(idx, weights=found, minlength=40)[:40]
eff = k / n; err = np.sqrt(eff * (1 - eff) / n); uc = 0.5 * (edges[1:] + edges[:-1])
a1.axvspan(0, 0.1, color=SHADE, lw=0, zorder=0)
a1.text(0.05, 50, "fit\nrange", ha="center", va="center", fontsize=7.5, color=MUTED)
a1.axhline(R["tracks"]["a"] * 100, color=BLUE, lw=0.9, ls=(0, (4, 2)), zorder=2)
a1.errorbar(uc, eff * 100, yerr=err * 100, fmt="o", ms=2.6, mfc=INK, mec=INK, ecolor=INK, elinewidth=0.7,
            capsize=0, zorder=3)
a1.set_xlim(0, 2 / 3); a1.set_ylim(0, 104)
a1.xaxis.set_major_locator(MultipleLocator(0.1)); a1.xaxis.set_minor_locator(MultipleLocator(0.05))
a1.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
a1.yaxis.set_major_locator(MultipleLocator(20)); a1.yaxis.set_minor_locator(MultipleLocator(10))
a1.set_xlabel(r"$1/p_T$ (GeV/$c$)$^{-1}$")
a1.set_ylabel(r"track-finding efficiency (\%)")
a1.text(0.27, 58, r"turn-on of the" "\n" r"2 GeV/$c$ OT $p_T$ cut", fontsize=7.5, color=MUTED, ha="center",
        va="center")
a1.annotate("", xy=(0.49, 72), xytext=(0.395, 62), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.6,
                                                                    mutation_scale=7))
a1.text(0.99, 93.5, r"$\varepsilon_\infty$", color=BLUE, fontsize=8, ha="right", va="top", transform=a1.get_yaxis_transform())
sec = a1.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
sec.xaxis.set_major_locator(FixedLocator([0, 0.1, 1 / 5, 1 / 3, 1 / 2, 2 / 3]))
sec.xaxis.set_major_formatter(FuncFormatter(lambda v, _: r"$\infty$" if v == 0 else f"{1/v:.3g}"))
sec.tick_params(direction="in", length=2.5, width=0.6, labelsize=7)
sec.set_xlabel(r"generated $p_T$ (GeV/$c$)", fontsize=7.5, labelpad=3)
a1.set_title("(a) whole sample", pad=17)
# (b) the fit
fit_panel(a2, tpt, found, R["tracks"], ylim=(98.3, 100.0), show_legend=True)
a2.set_xlabel(r"$x = 1/p_T^2$ ($10^{-3}$ (GeV/$c$)$^{-2}$)")
a2.set_ylabel(r"track-finding efficiency (\%)")
top_pt_axis(a2)
a2.set_title(r"(b) the fit, $p_T > 10$ GeV/$c$", pad=17)
fig.savefig(HERE / "fig_tracks.pdf")

# ---------------- Figure 2: hits, six subsystems ----------------
hpt, hp, hs, hev = d["hit_pt"].astype(float), d["hit_passed"].astype(float), d["hit_system"], d["hit_ev"]
fig, axs = plt.subplots(2, 3, figsize=(6.5, 6.1), gridspec_kw=dict(wspace=0.36, hspace=0.55))
for ax, s in zip(axs.flat, range(1, 7)):
    sel = hs == s
    f = R["hits"][str(s)]
    b = binned(hpt[sel], hp[sel], EDGES, hev[sel])
    lo = min(min(r[3] - r[4] for r in b), f["a"] - f["sa"]) * 100
    hi = max(max(r[3] + r[4] for r in b), f["a"] + f["sa"]) * 100
    pad = max(0.12, 0.5 * (hi - lo))
    bottom = lo - pad * 1.1
    top = hi + pad * 0.9
    if top > 100.0:            # an efficiency axis stops just above 100 %
        top = max(100.0, hi) + 0.03 * (100.0 - bottom)
    ylim = (bottom, top)
    fit_panel(ax, hpt[sel], hp[sel], f, groups=hev[sel], ylim=ylim, title=NAMES[s])
    top_pt_axis(ax, ticks=(10, 15, 30), label=False)
    ax.set_xlabel(r"$1/p_T^2$ ($10^{-3}$ (GeV/$c$)$^{-2}$)", fontsize=7.5)
    if s in (1, 4):
        ax.set_ylabel(r"signal hit efficiency (\%)")
fig.savefig(HERE / "fig_hits.pdf")
print(f"-> {HERE / 'fig_tracks.pdf'}, {HERE / 'fig_hits.pdf'}")
