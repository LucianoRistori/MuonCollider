"""make_cut_stages.py - what each hit-level cut buys, at the 100 um / 100 ps / 1.5 deg operating point.
1) fig_cut_stages: mean BIB hit density in the four parts of the IT and OT (hits per mm^2 per bunch
   crossing) with no cuts, the time cut only, the z0 cut only, and the two combined (bar chart, as in
   the run highlights).
2) cut_stages.json: those densities, and the expected fake tracks per bunch crossing, per region, at
   each stage (from regions.py runs on the BIB hits passing each stage; track |z0| < 4.5 mm).
Inputs: bib_stages.npz (dump_stages.py: IT and OT BIB hits with the pass bits 1 = time, 2 = z0,
4 = pT), density_per_layer_before_after_cuts.csv (layer areas), stage_<s>.json (regions.py).
Usage: python3 make_cut_stages.py <folder with bib_stages.npz and the csv> <folder with the jsons>"""
import csv, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.size": 10,
                     "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
                     "ytick.right": True, "legend.frameon": False})
IN = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
JS = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
SYS = ((3, "IT barrel"), (4, "IT endcaps"), (5, "OT barrel"), (6, "OT endcaps"))
STAGES = (("none", "no cuts", lambda b: np.ones(len(b), bool)),
          ("time", "time cut only", lambda b: (b & 1) > 0),
          ("z0", "$z_0$ cut only", lambda b: (b & 2) > 0),
          ("tz0", "time and $z_0$ cuts", lambda b: (b & 3) == 3))
COLS = ("#c6dbef", "#6baed6", "#9e9ac8", "#08306b")

d = np.load(IN / "bib_stages.npz"); s, bits = d["system"], d["bits"]
area = {}
for r in csv.DictReader(open(IN / "density_per_layer_before_after_cuts.csv")):
    area[int(r["system"])] = area.get(int(r["system"]), 0.0) + float(r["area_mm2"])
dens = {k: [float(((s == sid) & f(bits)).sum() / area[sid]) for sid, _ in SYS] for k, _, f in STAGES}
dens["all"] = [float(((s == sid) & (bits == 7)).sum() / area[sid]) for sid, _ in SYS]

fig, ax = plt.subplots(figsize=(7.2, 3.6))
x = np.arange(len(SYS)); w = 0.2
for i, (k, lab, _) in enumerate(STAGES):
    bars = ax.bar(x + (i - 1.5) * w, dens[k], w, color=COLS[i], edgecolor="0.25", linewidth=0.4, label=lab)
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() * 1.15, f"{b.get_height():.2g}",
                ha="center", va="bottom", fontsize=6.5, rotation=90)
ax.set_yscale("log"); ax.set_ylim(1e-4, 3)
ax.set_xticks(x); ax.set_xticklabels([n for _, n in SYS]); ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.5, len(SYS) - 0.5)
ax.set_ylabel(r"mean BIB hit density (mm$^{-2}$)")
ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.13), fontsize=8.5, handlelength=1.2)
ax.yaxis.grid(True, which="major", color="0.88", lw=0.5); ax.set_axisbelow(True)
fig.savefig("fig_cut_stages.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig("fig_cut_stages.png", dpi=180, bbox_inches="tight", pad_inches=0.02)

out = {"density_mm2": {k: dict(zip([n for _, n in SYS], v)) for k, v in dens.items()}, "fakes": {}}
for k in ("none", "time", "z0", "tz0", "all"):
    j = json.load(open(JS / f"stage_{k}.json"))
    out["fakes"][k] = {r: [j[r]["E"], j[r]["E_err"]] for r in ("barrel", "intermediate", "endcap")}
json.dump(out, open("cut_stages.json", "w"), indent=1)
for k, v in dens.items():
    print(f"{k:5s}", " ".join(f"{x:.3g}" for x in v))
for k, v in out["fakes"].items():
    print(f"{k:5s}", " ".join(f"{e[0]:.3g}" for e in v.values()), f"{sum(e[0] for e in v.values()):.3g}")
