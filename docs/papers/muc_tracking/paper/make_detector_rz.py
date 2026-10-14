"""make_detector_rz.py - quarter view (|z|, r) of the MuSIC_v2 inner and outer trackers:
barrel layers and endcap disks at their true positions and extents, each coloured by the
BIB hit density after the cuts (the densities used in the fake-track integral, as functions
of |z| on the barrel layers and of r on the disks, both ends together); the vertex detector
in grey; three straight example tracks from the origin, one per region, with their hits.
Inputs: regions.py (geometry and densities; REG_BIB, REG_DENS point to the run's files) and
the extents of the VXD layers (extents.csv).  Output: fig_detector_rz.pdf/png"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

sys.path.insert(0, os.environ.get("REGIONS_DIR", "/home/claude/gen"))
import regions as R

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.size": 10,
                     "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})
EXT = os.environ.get("REG_EXT", f"{R.UP}/claude_scratch/regions/extents.csv")
surfs = R.load_geometry()
norm = LogNorm(1e-4, 5e-2); cmap = plt.get_cmap("viridis_r")

fig, ax = plt.subplots(figsize=(7.2, 4.9))
segs, vals = [], []
for S in surfs:
    e = S["edges"]
    for a, b, rho in zip(e[:-1], e[1:], S["rho"]):
        if S["kind"] == "barrel":
            segs.append([(a, S["R"]), (b, S["R"])])
        else:
            segs.append([(S["zd"], a), (S["zd"], b)])
        vals.append(max(rho, 1.01e-4))
lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=3.2, capstyle="butt", zorder=3)
lc.set_array(np.array(vals)); ax.add_collection(lc)

# vertex detector, in grey
for e in csv.DictReader(open(EXT)):
    s, side = int(e["system"]), int(e["side"])
    if s == 1:
        r = 0.5 * (float(e["rmin"]) + float(e["rmax"])); ax.plot([0, float(e["zmax"])], [r, r], color="0.6", lw=1.2)
    elif s == 2 and side == 1:
        z = 0.5 * (float(e["zmin"]) + float(e["zmax"])); ax.plot([z, z], [float(e["rmin"]), float(e["rmax"])], color="0.6", lw=1.2)

# labels of the subsystems
ax.text(400, 75, "VXD", color="0.45", fontsize=9, ha="left", va="center")
ax.text(60, 255, "IT barrel", fontsize=9, ha="left")
ax.text(60, 900, "OT barrel", fontsize=9, ha="left")
ax.text(1180, 150, "IT disks", fontsize=9, ha="center")
ax.text(1800, 1480, "OT disks", fontsize=9, ha="center")

# example tracks (straight, from the origin) and their hits
tracks = (("barrel", 65.0, "#c0392b"), ("intermediate", 40.0, "#7b3294"), ("endcap", 15.0, "#e08214"))
for name, th, col in tracks:
    c = 1 / np.tan(np.radians(th))
    hits = []
    for S in surfs:
        if S["kind"] == "barrel" and S["R"] * c <= S["Z"]:
            hits.append((S["R"] * c, S["R"]))
        if S["kind"] == "disk" and S["r0"] <= S["zd"] / c <= S["r1"]:
            hits.append((S["zd"], S["zd"] / c))
    hits = np.array(sorted(hits))
    rmax = min(1650, 2400 / c)
    ax.plot([0, rmax * c], [0, rmax], color=col, lw=1.1, zorder=5)
    ax.plot(hits[:, 0], hits[:, 1], "o", ms=4.5, mfc="white", mec=col, mew=1.3, zorder=6)
    lab = {"barrel": (790, 1560, "left"), "intermediate": (1980, 1560, "right"), "endcap": (2400, 690, "right")}[name]
    ax.text(lab[0], lab[1], f"{name} ($\\vartheta={th:.0f}^\\circ$)", color=col, fontsize=9,
            ha=lab[2], va="bottom")

ax.set_xlim(0, 2450); ax.set_ylim(0, 1650); ax.set_aspect("equal")
ax.set_xlabel(r"$|z|$ (mm)"); ax.set_ylabel(r"$r$ (mm)")
cb = fig.colorbar(lc, ax=ax, fraction=0.035, pad=0.02)
cb.set_label(r"BIB hit density after the cuts (mm$^{-2}$)")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"fig_detector_rz.{ext}", dpi=200, bbox_inches="tight", pad_inches=0.02)
print("done")
