"""
make_example_figures.py - figures of the worked example (paper 1, Section 6),
from worked_example.json (python3 worked_example.py).

  example_density.pdf     E[#fakes] vs a common factor mu on all six hit densities
  example_resolution.pdf  E[#fakes] vs the hit position resolution (both coordinates),
                          at fixed densities: each m-layer term scales as sigma^(2(m-2))
Both show the at-least-5-of-6 sum, the six 5-of-6 terms and the 6-of-6 term.
PNG copies are written for quick viewing.
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = Path(__file__).resolve().parent
ex = json.load(open(HERE / "worked_example.json"))
SIG0 = ex["sigma"] * 1000.0                         # um

terms = []                                          # (label, m, E at mu = 1, sigma0)
for key, v in ex["sets"].items():
    S = [int(x) for x in key.split(",")]
    miss = sorted(set(range(6)) - set(S))
    label = "6 of 6" if not miss else f"5 of 6, without layer {miss[0] + 1}"
    terms.append((label, len(S), v["E"]))
five = sorted([t for t in terms if t[1] == 5], key=lambda t: -t[2])
six = [t for t in terms if t[1] == 6][0]

plt.rcParams.update({"font.size": 10, "axes.labelsize": 10.5, "legend.fontsize": 8,
                     "font.family": "serif", "mathtext.fontset": "dejavuserif"})
C_SUM, C_SIX = "#1b4f72", "#c0392b"
greys = plt.cm.Greys(np.linspace(0.75, 0.4, len(five)))


def crossing(x, y):
    i = np.nonzero(y >= 1)[0]
    if len(i) == 0 or i[0] == 0:
        return None
    j = i[0]
    lx = np.interp(0.0, np.log(y[j - 1:j + 1]), np.log(x[j - 1:j + 1]))
    return math.exp(lx)


def finish(ax, x0, x0_label, xlabel, out, ticks):
    ax.axhline(1.0, color="0.5", lw=0.8)
    ax.axvline(x0, color="0.3", lw=0.8, ls=":")
    ax.text(x0, 0.02, " " + x0_label, transform=ax.get_xaxis_transform(), fontsize=8, va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(ticks); ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"expected fake tracks per tower, $E$")
    ax.grid(True, which="major", alpha=0.3); ax.grid(True, which="minor", alpha=0.1)
    ax.legend(loc="upper left", frameon=False)
    fig = ax.figure
    fig.tight_layout()
    fig.savefig(HERE / f"{out}.pdf")
    fig.savefig(HERE / f"{out}.png", dpi=150)
    plt.close(fig)
    print(f"wrote {out}.pdf, {out}.png")


def plot(x, scale, x0, x0_label, xlabel, out, ylim, ticks):
    """scale(m): factor multiplying an m-layer term at the values x."""
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    tot = sum(E * scale(m) for _, m, E in terms)
    for (label, m, E), g in zip(five, greys):
        ax.plot(x, E * scale(m), color=g, lw=0.9)
    ax.plot([], [], color=greys[1], lw=0.9, label="each of the six 5-of-6 terms")
    ax.plot(x, six[2] * scale(6), color=C_SIX, lw=1.4, ls="--", label="6 of 6")
    ax.plot(x, tot, color=C_SUM, lw=2.2, label="at least 5 of 6 (sum)")
    xc = crossing(x, tot)
    if xc:
        ax.plot([xc], [1.0], "o", color=C_SUM, ms=6, mfc="white", mew=1.6, zorder=5)
        ax.annotate(f"one fake at {xc:.3g}" + ("$\\times$" if out == "example_density" else " $\\mu$m"),
                    (xc, 1.0), xytext=(-8, 8), textcoords="offset points", fontsize=8.5, color=C_SUM, ha="right")
    ax.set_ylim(*ylim)
    finish(ax, x0, x0_label, xlabel, out, ticks)


mu = np.geomspace(0.3, 30, 400)
plot(mu, lambda m: mu ** m, 1.0, "example densities",
     r"common factor $\mu$ on the hit density of all six layers", "example_density", (1e-9, 1e4), [0.3, 1, 3, 10, 30])

sig = np.geomspace(20, 500, 400)                     # um
plot(sig, lambda m: (sig / SIG0) ** (2 * (m - 2)), SIG0, f"{SIG0:g} $\\mu$m",
     r"hit position resolution $\sigma$ in both coordinates ($\mu$m)", "example_resolution", (1e-9, 1e4), [20, 50, 100, 200, 500])
