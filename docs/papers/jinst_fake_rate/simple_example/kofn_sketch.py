"""kofn_sketch.py - schematic (one view, not to scale) for the revised section 2.5:
a candidate on planes 1-5, the replacement window on plane 3 with a second hit
in it (an alternative candidate sharing four hits), and the extension window on
plane 6. Output: kofn_sketch.pdf / .png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.size": 11})
Y = np.array([1, 2, 3, 4, 5, 6.0]); H = 0.25 + 0.25 * Y
RED, BLUE, GREY = "#c0392b", "#1f4e8c", "0.45"
rng = np.random.default_rng(5)
fig, ax = plt.subplots(figsize=(7.0, 3.3))
for y, h in zip(Y, H):
    ax.plot([y, y], [-h, h], color=BLUE, lw=1.6)
    ax.text(y, h + 0.12, f"${int(y)}$", ha="center", va="bottom", fontsize=12)
    m = rng.uniform(-h, h, int(4 + 2 * y))
    m = m[np.abs(m - (0.05 + 0.11 * y)) > 0.4]            # keep the windows free of other hits
    ax.plot(np.full_like(m, y), m, "o", ms=3, color=GREY, zorder=3)

a, b = 0.05, 0.11                                            # the candidate's line u = a + b Y
dev = {1: 0.03, 2: -0.04, 3: 0.05, 4: -0.03, 5: 0.02}
hits = {int(y): a + b * y + dev[int(y)] for y in Y[:5]}
yy = np.linspace(0.3, 6.3, 50)
ax.plot(yy, a + b * yy, color=RED, lw=1.5, zorder=4)
for y, u in hits.items():
    ax.plot(y, u, "o", ms=6, color=RED, zorder=5)

# replacement window on plane 3 (fit to planes 1, 2, 4, 5 predicts plane 3)
w3 = 0.22
ax.add_patch(plt.Rectangle((3 - 0.07, a + 3 * b - w3), 0.14, 2 * w3, color=RED, alpha=0.18, lw=0, zorder=2))
alt = a + 3 * b - 0.13
ax.plot(3, alt, "o", ms=6, mfc="white", mec=RED, mew=1.5, zorder=6)
ax.plot([2, 3, 4], [hits[2], alt, hits[4]], color=RED, lw=1.0, ls=(0, (3, 2)), zorder=4)
ax.annotate("replacement window\n(another hit here: a second\ncandidate sharing four hits)",
            xy=(2.93, a + 3 * b - w3 + 0.02), xytext=(0.35, -1.75), fontsize=9.5, ha="left",
            arrowprops=dict(arrowstyle="-", color="0.3", lw=0.7))

# extension window on plane 6 (fit to planes 1-5 predicts plane 6)
w6 = 0.34
ax.add_patch(plt.Rectangle((6 - 0.07, a + 6 * b - w6), 0.14, 2 * w6, color=BLUE, alpha=0.18, lw=0, zorder=2))
ax.annotate("extension window\n(a hit here: a candidate\non all six planes)", xy=(5.93, a + 6 * b + w6 - 0.02),
            xytext=(4.55, 1.95), fontsize=9.5, ha="left", arrowprops=dict(arrowstyle="-", color="0.3", lw=0.7))

ax.annotate("", xy=(6.55, 0), xytext=(0.2, 0), arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
ax.text(6.6, 0, "$Y$", va="center", fontsize=12)
ax.annotate("", xy=(0.2, 1.35), xytext=(0.2, -1.35), arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
ax.text(0.2, 1.42, "$u$", ha="center", va="bottom", fontsize=12)
ax.set_xlim(0, 7.0); ax.set_ylim(-1.9, 2.5); ax.axis("off")
fig.tight_layout(pad=0.1)
for ext in ("pdf", "png"):
    fig.savefig(f"kofn_sketch.{ext}", dpi=200, bbox_inches="tight", pad_inches=0.03)
