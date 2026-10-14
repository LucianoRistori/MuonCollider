"""make_detector_3d_pv.py - 3-D view of the MuSIC_v2 inner and outer trackers, rendered with
PyVista (VTK): one half (z >= 0) with the quarter facing the viewer cut away, opaque layers with
real lighting and dark outlines, and three helices from the interaction point, one per region.
Usage: python3 make_detector_3d_pv.py [a|b]
  a: z to the right, y up, x pointing away from the viewer
  b: z to the left,  y up, x pointing towards the viewer
Output: fig_detector_3d_<a|b>.png"""
import os, sys
import numpy as np
import pyvista as pv

sys.path.insert(0, os.environ.get("REGIONS_DIR", "/home/claude/gen"))
import regions as R

MODE = sys.argv[1] if len(sys.argv) > 1 else "a"
pv.OFF_SCREEN = True
surfs = R.load_geometry()
B_T = 5.0
if MODE == "a":
    V = np.array([-0.80, 0.45, -0.40])
else:
    V = np.array([0.80, 0.45, -0.40])
V = V / np.linalg.norm(V)
PHI_CAM = np.degrees(np.arctan2(V[1], V[0]))
ROT = float(os.environ.get("CUT_ROT", "15"))        # detector (cut) rotated by +ROT deg in phi: clockwise as seen from the camera
TRK = float(os.environ.get("TRK_GAP", "-2"))         # tracks start TRK deg outside the upper edge of the cut, on drawn layers
CUT_C = PHI_CAM + ROT
PHI_LO, PHI_HI = np.radians(CUT_C + 45), np.radians(CUT_C + 315)     # the quarter around CUT_C is cut away
IT_COL, OT_COL = "#b3cdee", "#7d9fd0"
DISK_IT = ["#e3c9a0", "#dcbe91", "#d5b383", "#cea875", "#c79d67", "#c09259", "#b9874b"]
DISK_OT = ["#c9b7d9", "#b9a4cc", "#a991bf", "#997eb2"]
EDGE = "#1d2f4f"

pl = pv.Plotter(off_screen=True, window_size=(2200, 1500))
pl.set_background("white")
pl.enable_anti_aliasing("ssaa")


def sector_surface(r_or_R, phis, z_or_zs, kind):
    if kind == "cyl":
        Rr, zs = r_or_R, z_or_zs
        P, Z = np.meshgrid(phis, zs, indexing="ij")
        X, Y = Rr * np.cos(P), Rr * np.sin(P)
    else:
        rs, zd = r_or_R, z_or_zs
        P, Rg = np.meshgrid(phis, rs, indexing="ij")
        X, Y, Z = Rg * np.cos(P), Rg * np.sin(P), np.full_like(P, zd)
    return pv.StructuredGrid(X, Y, Z).extract_surface()


phis = np.linspace(PHI_LO, PHI_HI, 220)
it_d = ot_d = 0
for S in surfs:
    if S["kind"] == "barrel":
        mesh = sector_surface(S["R"], phis, np.linspace(0, S["Z"], 40), "cyl")
        col = IT_COL if S["name"].startswith("IT") else OT_COL
    else:
        mesh = sector_surface(np.linspace(S["r0"], S["r1"], 30), phis, S["zd"], "disk")
        if S["name"].startswith("IT"):
            col = DISK_IT[it_d]; it_d += 1
        else:
            col = DISK_OT[ot_d]; ot_d += 1
    pl.add_mesh(mesh, color=col, smooth_shading=True, specular=0.1, ambient=0.4, diffuse=0.7)
    pl.add_mesh(mesh.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                           non_manifold_edges=False), color=EDGE, line_width=2.2)


def helix_pts(phi0, cot, pt, s):
    k = 0.2998 * B_T / pt / 1000.0
    x = (np.sin(phi0 + k * s) - np.sin(phi0)) / k
    y = (np.cos(phi0) - np.cos(phi0 + k * s)) / k
    return np.stack([x, y, cot * s], -1)


def crossings(phi0, cot, pt):
    k = 0.2998 * B_T / pt / 1000.0; out = []
    for S in surfs:
        if S["kind"] == "barrel":
            s = 2 / k * np.arcsin(k * S["R"] / 2)
            if cot * s <= S["Z"]:
                out.append(s)
        else:
            s = S["zd"] / cot
            r = np.hypot(*helix_pts(phi0, cot, pt, np.array([s]))[0, :2])
            if S["r0"] <= r <= S["r1"]:
                out.append(s)
    return np.array(sorted(out))


phi0 = np.radians(CUT_C - 45 - TRK)                 # just past the upper edge of the cut (mode a)
Q_PT = -4.0                                          # negative charge: the track bends away from the cut
labels, lab_pts = [], []
for name, th, col in (("barrel", 65.0, "#c0392b"), ("intermediate", 40.0, "#7b3294"), ("endcap", 15.0, "#e08214")):
    cot = 1 / np.tan(np.radians(th))
    sh = crossings(phi0, cot, Q_PT)
    s = np.linspace(0, sh.max() * 1.1, 300)
    pts = helix_pts(phi0, cot, Q_PT, s)
    pl.add_mesh(pv.Spline(pts, 300).tube(radius=9), color=col, smooth_shading=True)
    for q in helix_pts(phi0, cot, Q_PT, sh):
        pl.add_mesh(pv.Sphere(radius=24, center=q), color=col, smooth_shading=True)
    labels.append(name); lab_pts.append(pts[-1] + np.array([0, 60, 40]))
    for sc, tc, mg in (("black", "black", 8), ("white", col, 6)):     # black rectangle under a white one: a thin border
        pl.add_point_labels([pts[-1] + np.array([0, 70, 40])], [name], font_size=34, text_color=tc, font_family="times",
                            shape="rect", fill_shape=True, shape_color=sc, shape_opacity=1.0, margin=mg,
                            show_points=False, always_visible=True, bold=False)

# beam line and axes
pl.add_mesh(pv.Line((0, 0, -300), (0, 0, 2500)), color="#555555", line_width=2)
for d, lab in (((1, 0, 0), "x"), ((0, 1, 0), "y"), ((0, 0, 1), "z")):
    v = np.array(d, float)
    pl.add_mesh(pv.Arrow(start=(0, 0, 0), direction=v, scale=450, tip_length=0.18, tip_radius=0.05, shaft_radius=0.015),
                color="black")
    pl.add_point_labels([v * 540], [lab], font_size=40, text_color="black", font_family="times", italic=True,
                        shape=None, show_points=False, always_visible=True)

focal = np.array([0.0, 0.0, 1150.0])
pl.camera.focal_point = focal
pl.camera.position = focal + V * 9000
pl.camera.up = (0, 1, 0)
pl.camera.view_angle = 30
img = pl.screenshot(return_img=True)
from PIL import Image
im = Image.fromarray(img)
import numpy as _np
m = (_np.asarray(im)[..., :3] < 250).any(2)
r = _np.nonzero(m.any(1))[0]; c = _np.nonzero(m.any(0))[0]
im.crop((c[0] - 20, r[0] - 20, c[-1] + 20, r[-1] + 20)).save(os.environ.get("OUT_PNG", f"fig_detector_3d_{MODE}.png"))
print("done", MODE, round(PHI_CAM, 1))
