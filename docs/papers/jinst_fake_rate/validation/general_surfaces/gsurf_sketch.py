"""gsurf_sketch.py - 3-D sketches of the four general-surface validation detectors
(geoms.py, lengths in units of L), in the style of planes_sketch.py: the four
surfaces, a few uniformly distributed noise hits on each, one helix crossing
all four (red) with its hits, surface labels, the axes and the beam line.
Output: gsurf_sketch_<name>.pdf / .png"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import gsurf as g
from geoms import GEOMS

plt.rcParams.update({"mathtext.fontset": "cm", "font.size": 11})
FACE = (0.55, 0.70, 0.90, 0.22); EDGE = (0.15, 0.30, 0.55, 0.9); RED = "#c0392b"
VIEWS = {"cylinders": (16, -42), "tilted45z": (48, -86), "tilted60_endcaplike": (32, -78), "mixed": (16, -38)}

FS = {"cylinders": 1.053, "tilted45z": 1.043, "tilted60_endcaplike": 1.064, "mixed": 0.86}   # same printed label size as figure 1
AXOFF = {"cylinders": {"y": (0.0, 0.25, -0.25)},
         "tilted45z": {"z": (-0.25, 0.0, 0.1), "y": (0.0, 0.15, 0.25)},
         "tilted60_endcaplike": {"z": (-0.15, 0.0, 0.1), "y": (0.0, 0.2, 0.25)}}


def view_vec(elev, azim):
    e, a = np.radians(elev), np.radians(azim)
    return np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


def outline(S, VIEW, n=60):
    """closed outline (n_pts, 3) and a mesh of quads for a surface"""
    if isinstance(S, g.Cylinder):
        ph = S.phi_c + np.linspace(-S.half_dphi, S.half_dphi, n)
        x, y = S.R * np.cos(ph), S.R * np.sin(ph)
        lo = np.stack([x, y, np.full(n, S.z_lo)], 1); hi = np.stack([x, y, np.full(n, S.z_hi)], 1)
        ring = np.concatenate([lo, hi[::-1], lo[:1]])
        quads = [[lo[i], lo[i + 1], hi[i + 1], hi[i]] for i in range(n - 1)]
        top = hi[n // 2]
        return ring, quads, top
    c = [(S.u_lo, S.v_lo), (S.u_hi, S.v_lo), (S.u_hi, S.v_hi), (S.u_lo, S.v_hi), (S.u_lo, S.v_lo)]
    pts = np.array([S.P0 + u * S.eu + v * S.ev for u, v in c])
    vmid = 0.5 * (S.v_lo + S.v_hi); umid = 0.5 * (S.u_lo + S.u_hi)
    # "top": the edge point highest in z (or largest u for horizontal planes)
    cand = np.array([S.P0 + umid * S.eu + S.v_hi * S.ev, S.P0 + umid * S.eu + S.v_lo * S.ev,
                     S.P0 + S.u_hi * S.eu + vmid * S.ev, S.P0 + S.u_lo * S.eu + vmid * S.ev])
    top = cand[np.argmax(cand[:, 2] + 1e-3 * cand @ VIEW)]
    return pts, [pts[:4].tolist()], top


def pick_track(surfs, box, kmax, rng, kfrac=0.6):
    lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
    th = lo + (hi - lo) * rng.random((400_000, 4))
    th[:, 1] = np.sign(th[:, 1]) * kfrac * kmax
    th[:, 2] = 0.0                                   # the track starts at the origin
    ins = np.all([S.cross(th)[3] for S in surfs], 0)
    th = th[ins]
    # most central: the crossings far from the edges, phi0 near 0
    score = np.abs(th[:, 0]) + np.abs(th[:, 2]) / (10 * abs(box[2][1]))
    return th[np.argmin(score)]


def helix_path(th, smax, n=300):
    s = np.linspace(0, smax, n)
    return g.helix_point(np.repeat(th[None], n, 0), s)


def draw(name, L=1.0, scale_axes=0.8, label_off=0.08, fname=None, extra=None):
    ELEV, AZIM = VIEWS[name]; VIEW = view_vec(ELEV, AZIM)
    surfs, kmax, box = GEOMS[name](L)
    rng = np.random.default_rng(3)
    th = pick_track(surfs, box, kmax, rng)
    fig = plt.figure(figsize=(5.2, 4.2))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    ax.set_proj_type("ortho"); ax.view_init(elev=ELEV, azim=AZIM)
    allpts = [np.zeros((1, 3))]
    depth = []
    for S in surfs:
        ring, quads, top = outline(S, VIEW)
        depth.append(np.mean(ring, 0) @ VIEW)
    order = np.argsort(np.argsort(depth))
    amax = max(S.area for S in surfs)
    for i, S in enumerate(surfs):
        ring, quads, top = outline(S, VIEW)
        z0 = 2 * order[i]
        ax.add_collection3d(Poly3DCollection(quads, facecolor=FACE, edgecolor="none", zorder=z0))
        ax.plot(*ring.T, color=EDGE, lw=0.9, zorder=z0 + 0.2)
        nn = int(16 * S.area / amax) + 4
        u, v = S.random_hits(rng, nn)
        p = g.local_to_global(S, u, v)
        ax.scatter(*p.T, s=4, c="0.35", depthshade=False, zorder=z0 + 0.5)
        lab = top + label_off * L * np.array([0, 0, 1.0])
        ax.text(*lab, f"${i + 1}$", ha="center", va="bottom", fontsize=12 * FS[name], zorder=50)
        allpts.append(ring)
    # the track and its hits
    hits = []
    for S in surfs:
        u, v, ok, inside = S.cross(th[None])
        hits.append(g.local_to_global(S, u, v)[0])
    hits = np.array(hits)
    smax = 1.12 * max(np.hypot(hits[:, 0], hits[:, 1]).max() * 1.0,
                      ((hits[:, 2] - th[2]) / th[3]).max() if abs(th[3]) > 1e-6 else 0)
    path = helix_path(th, smax)
    ax.plot(*path.T, color=RED, lw=1.6, zorder=60)
    ax.scatter(*hits.T, s=16, c=RED, depthshade=False, zorder=61)
    allpts.append(path)
    # beam line and axes at the origin
    P = np.concatenate(allpts)
    zlo, zhi = P[:, 2].min() - 0.3 * L, P[:, 2].max() + 0.3 * L
    ax.plot([0, 0], [0, 0], [zlo, zhi], color="0.5", lw=0.8, ls=(0, (4, 3)), zorder=1)
    A = scale_axes * L
    offs = {"x": (0.05, -0.18, -0.18), "y": (0.0, 0.1, 0.18), "z": (0.0, -0.12, 0.12)}
    offs.update(AXOFF.get(name, {}))
    for d, lab, off in (((1, 0, 0), r"$x$", offs["x"]), ((0, 1, 0), r"$y$", offs["y"]),
                        ((0, 0, 1), r"$z$", offs["z"])):
        ax.quiver(0, 0, 0, *(A * np.array(d)), color="k", lw=1.0, arrow_length_ratio=0.15, zorder=70)
        ax.text(*(A * np.array(d) + L * np.array(off)), lab, fontsize=14 * FS[name], ha="center", va="center", zorder=70)
    if extra:
        extra(ax, L)
    # equal scales
    P = np.concatenate([P, [[0, 0, zlo], [0, 0, zhi]]])
    c = 0.5 * (P.min(0) + P.max(0)); r = 0.5 * (P.max(0) - P.min(0)) + 0.1 * L
    ax.set_xlim(c[0] - r[0], c[0] + r[0]); ax.set_ylim(c[1] - r[1], c[1] + r[1]); ax.set_zlim(c[2] - r[2], c[2] + r[2])
    ax.set_box_aspect(tuple(2 * r), zoom=1.0)
    ax.set_axis_off(); fig.subplots_adjust(0, 0, 1, 1)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
    rows = np.nonzero((img < 250).any(axis=(1, 2)))[0]; cols = np.nonzero((img < 250).any(axis=(0, 2)))[0]
    H = img.shape[0]; dpi = fig.dpi; m = 6
    bb = Bbox([[(cols[0] - m) / dpi, (H - rows[-1] - m) / dpi], [(cols[-1] + m) / dpi, (H - rows[0] + m) / dpi]])
    for ext in ("pdf", "png"):
        fig.savefig(f"{fname or 'gsurf_sketch_' + name}.{ext}", dpi=200, bbox_inches=bb)
    plt.close(fig)
    print(name, "track", np.round(th, 5))


if __name__ == "__main__":
    for nm in GEOMS:
        draw(nm)
