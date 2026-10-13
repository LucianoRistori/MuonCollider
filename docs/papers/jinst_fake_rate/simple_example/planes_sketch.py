"""planes_sketch.py - 3-D sketch of the simple case: six parallel square planes
(the projective tower of the numerical example, Y = 250..1500 mm, W = Y/1.5),
a straight track crossing them, a few noise hits, plane labels and axes.
Output: planes_sketch.pdf / .png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
plt.rcParams.update({"mathtext.fontset": "cm", "font.size": 11})

Y = np.array([250, 500, 750, 1000, 1250, 1500.0]); W = Y / 1.5
rng = np.random.default_rng(7)

fig = plt.figure(figsize=(7.0, 5.0))
ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
ax.set_proj_type("ortho")
ax.view_init(elev=16, azim=-58)

# drawing coordinates: x = Y (depth), y = u, z = v
for i, (y, w) in enumerate(zip(Y, W)):
    h = w / 2
    quad = [[(y, -h, -h), (y, h, -h), (y, h, h), (y, -h, h)]]
    ax.add_collection3d(Poly3DCollection(quad, facecolor=(0.55, 0.70, 0.90, 0.20),
                                         edgecolor=(0.15, 0.30, 0.55, 0.9), lw=0.9, zorder=i))
    nn = int(18 * (w / W[-1]) ** 2) + 4                               # noise hits, uniform
    uv = rng.uniform(-h, h, (nn, 2))
    ax.scatter(np.full(nn, y), uv[:, 0], uv[:, 1], s=4, c="0.35", depthshade=False, zorder=i + 0.5)
    ax.text(y, 0, h + 40, f"${i + 1}$", ha="center", va="bottom", fontsize=12, zorder=20)

# the track u = a_u + b_u Y, v = a_v + b_v Y
au, bu, av, bv = 0.0, 0.24, 0.0, 0.20
yy = np.array([0.0, 1650.0])
ax.plot(yy, au + bu * yy, av + bv * yy, color="#c0392b", lw=1.6, zorder=30)
ax.scatter(Y, au + bu * Y, av + bv * Y, s=16, c="#c0392b", depthshade=False, zorder=31)

# axes at the origin
L = 170
for d, lab, off in (((1, 0, 0), r"$Y$", (-20, -40, -55)), ((0, 1, 0), r"$u$", (0, 10, 60)),
                    ((0, 0, 1), r"$v$", (0, 0, 60))):
    ax.quiver(0, 0, 0, *(L * np.array(d)), color="k", lw=1.0, arrow_length_ratio=0.14, zorder=40)
    ax.text(*(L * np.array(d) + np.array(off)), lab, fontsize=14, ha="center", va="center", zorder=40)

# dimension marks: side W_6 and depth Y_6, below the tower on the near side
y6, h6 = Y[-1], W[-1] / 2
zf = -h6 - 80
ax.plot([y6, y6], [-h6, h6], [zf, zf], color="0.3", lw=0.8)
for s_ in (-h6, h6):
    ax.plot([y6, y6], [s_, s_], [zf - 20, zf + 20], color="0.3", lw=0.8)
ax.text(y6, 60, zf - 40, r"$W_6$", ha="center", va="top", fontsize=14)
ax.plot([0, y6], [-h6, -h6], [zf, zf], color="0.3", lw=0.8)
ax.plot([0, 0], [-h6, -h6], [zf - 20, zf + 20], color="0.3", lw=0.8)
for (y0, u0, v0), (y1, u1, v1) in (((0, 0, 0), (0, -h6, 0)), ((0, -h6, 0), (0, -h6, zf)),
                                   ((y6, -h6, -h6), (y6, -h6, zf))):
    ax.plot([y0, y1], [u0, u1], [v0, v1], color="0.65", lw=0.6, ls=(0, (2, 2)))
ax.text(y6 / 2, -h6, zf - 40, r"$Y_6$", ha="center", va="top", fontsize=14)

ax.set_xlim(0, 1650); ax.set_ylim(-560, 560); ax.set_zlim(-640, 560)
ax.set_box_aspect((1650, 1120, 1200), zoom=1.0)
ax.set_axis_off()
fig.subplots_adjust(0, 0, 1, 1)
# crop to the drawn content (3-D axes ignore bbox_inches="tight")
from matplotlib.transforms import Bbox
fig.canvas.draw()
img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
rows = np.nonzero((img < 250).any(axis=(1, 2)))[0]; cols = np.nonzero((img < 250).any(axis=(0, 2)))[0]
H, Wpx = img.shape[:2]; dpi = fig.dpi; m = 6
bb = Bbox([[(cols[0] - m) / dpi, (H - rows[-1] - m) / dpi], [(cols[-1] + m) / dpi, (H - rows[0] + m) / dpi]])
for ext in ("pdf", "png"):
    fig.savefig(f"planes_sketch.{ext}", dpi=200, bbox_inches=bb)
