"""make_towers_global_pv.py - figure 3 of paper 2: one "tower" per region, rendered with PyVista.
Each tower is the set of pieces of the real MuSIC_v2 layers crossed by a narrow bundle of
tracks from the interaction point (polar angle theta +- 4 deg, azimuth +- 6 deg): barrel
(theta = 65 deg: six barrel cylinders), intermediate (40 deg: barrel cylinders and OT disks),
endcap (15 deg: IT disks). Each piece carries a few uniformly distributed background hits
(illustrative numbers, not to scale) and one track crosses the tower with its hits.
Colours as in fig_detector_3d. All three towers in one scene, at their true positions, with global axes. Output: fig_towers.png."""
import os, sys
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, os.environ.get("REGIONS_DIR", "/home/claude/gen"))
import regions as R

pv.OFF_SCREEN = True
surfs = R.load_geometry()
BAR = {"IT": "#b3cdee", "OT": "#7d9fd0"}
DISK_IT = ["#e3c9a0", "#dcbe91", "#d5b383", "#cea875", "#c79d67", "#c09259", "#b9874b"]
DISK_OT = ["#c9b7d9", "#b9a4cc", "#a991bf", "#997eb2"]
EDGE = "#1d2f4f"
DTH, DPHI = np.radians(6.0), np.radians(9.0)
rng = np.random.default_rng(11)


def patch_barrel(Rr, Z, th):
    z0, z1 = Rr / np.tan(th + DTH), Rr / np.tan(th - DTH)
    z0, z1 = max(z0, 0.0), min(z1, Z)
    if z1 <= z0:
        return None
    P, Zg = np.meshgrid(np.linspace(-DPHI, DPHI, 30), np.linspace(z0, z1, 12), indexing="ij")
    mesh = pv.StructuredGrid(Rr * np.cos(P), Rr * np.sin(P), Zg).extract_surface(algorithm=None)
    area = Rr * 2 * DPHI * (z1 - z0)
    def sample(n):
        ph = rng.uniform(-DPHI, DPHI, n); z = rng.uniform(z0, z1, n)
        return np.stack([Rr * np.cos(ph), Rr * np.sin(ph), z], -1)
    return mesh, area, sample


def patch_disk(zd, r0, r1, th):
    a, b = zd * np.tan(th - DTH), zd * np.tan(th + DTH)
    a, b = max(a, r0), min(b, r1)
    if b <= a:
        return None
    P, Rg = np.meshgrid(np.linspace(-DPHI, DPHI, 30), np.linspace(a, b, 12), indexing="ij")
    mesh = pv.StructuredGrid(Rg * np.cos(P), Rg * np.sin(P), np.full_like(P, zd)).extract_surface(algorithm=None)
    area = DPHI * (b * b - a * a)
    def sample(n):
        ph = rng.uniform(-DPHI, DPHI, n); r = np.sqrt(rng.uniform(a * a, b * b, n))
        return np.stack([r * np.cos(ph), r * np.sin(ph), np.full(n, zd)], -1)
    return mesh, area, sample


PHI_C = float(os.environ.get("TOWER_PHI", "100"))       # azimuth of the towers (deg), as the tracks of fig 2
c, s_ = np.cos(np.radians(PHI_C)), np.sin(np.radians(PHI_C))
ROTM = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
pl = pv.Plotter(off_screen=True, window_size=(2600, 1800))
pl.set_background("white"); pl.enable_anti_aliasing("ssaa")
Lmax = 0
for name, th_deg, col_track in (("barrel", 65.0, "#c0392b"), ("intermediate", 40.0, "#7b3294"), ("endcap", 15.0, "#e08214")):
    th = np.radians(th_deg)
    it_d = ot_d = 0
    d = ROTM @ np.array([np.sin(th), 0.0, np.cos(th)])
    hits = []
    for S in surfs:
        if S["kind"] == "barrel":
            p = patch_barrel(S["R"], S["Z"], th)
            col = BAR[S["name"][:2]]
            hit = (S["R"], 0.0, S["R"] / np.tan(th)) if S["R"] / np.tan(th) <= S["Z"] else None
        else:
            if S["name"].startswith("IT"):
                col = DISK_IT[it_d]; it_d += 1
            else:
                col = DISK_OT[ot_d]; ot_d += 1
            p = patch_disk(S["zd"], S["r0"], S["r1"], th)
            rr = S["zd"] * np.tan(th)
            hit = (rr, 0.0, S["zd"]) if S["r0"] <= rr <= S["r1"] else None
        if p is None or hit is None:
            continue
        mesh, area, sample = p
        mesh = mesh.rotate_z(PHI_C, point=(0, 0, 0))
        pl.add_mesh(mesh, color=col, smooth_shading=True, ambient=0.4, diffuse=0.7, specular=0.1)
        pl.add_mesh(mesh.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                               non_manifold_edges=False), color=EDGE, line_width=2.5)
        nn = int(np.clip(area / 2.5e3, 6, 40))
        pts = sample(nn) @ ROTM.T
        pl.add_mesh(pv.PolyData(pts).glyph(geom=pv.Sphere(radius=5.5), scale=False, orient=False),
                    color="#555555", smooth_shading=True)
        hits.append(ROTM @ np.array(hit))
    hits = np.array(hits); L = np.linalg.norm(hits, axis=1).max() * 1.08; Lmax = max(Lmax, L)
    pl.add_mesh(pv.Line((0, 0, 0), tuple(d * L)).tube(radius=6), color=col_track, smooth_shading=True)
    pl.add_mesh(pv.PolyData(hits).glyph(geom=pv.Sphere(radius=17), scale=False, orient=False),
                color=col_track, smooth_shading=True)
    pl.add_point_labels([d * L * 1.04 + np.array([0, 40, 0])], [f"{name} ({th_deg:.0f}°)"], font_size=40,
                        text_color=col_track, font_family="times", shape=None, show_points=False, always_visible=True)

# beam line, interaction point and global axes
pl.add_mesh(pv.Line((0, 0, -400), (0, 0, 2500)), color="#555555", line_width=2)
pl.add_mesh(pv.Sphere(radius=14, center=(0, 0, 0)), color="black")
for dv, lab in (((1, 0, 0), "x"), ((0, 1, 0), "y"), ((0, 0, 1), "z")):
    v = np.array(dv, float)
    pl.add_mesh(pv.Arrow(start=(0, 0, 0), direction=v, scale=400, tip_length=0.2, tip_radius=0.06, shaft_radius=0.018),
                color="black")
    pl.add_point_labels([v * 490], [lab], font_size=44, text_color="black", font_family="times", italic=True,
                        shape=None, show_points=False, always_visible=True)

V = np.array([-0.80, 0.45, -0.40]); V /= np.linalg.norm(V)       # as fig_detector_3d (a)
focal = np.array([0.0, 750.0, 1000.0])
pl.camera.focal_point = focal
pl.camera.position = focal + V * 9000
pl.camera.up = (0, 1, 0)
pl.camera.view_angle = 26
img = pl.screenshot(return_img=True)
m = (img[..., :3] < 250).any(2); r = np.nonzero(m.any(1))[0]; cc = np.nonzero(m.any(0))[0]
Image.fromarray(img).crop((cc[0] - 30, r[0] - 30, cc[-1] + 30, r[-1] + 30)).save(os.environ.get("OUT_PNG", "fig_towers.png"))
print("done")
