"""Test geometries (4 surfaces each, sigma = 0.1 mm), all lengths proportional to L (mm)."""
import numpy as np
import gsurf as g
S = 0.1


def cylinders(L):
    surfs = [g.Cylinder(R, 0.0, 0.12, -1.5 * L, 1.5 * L, S) for R in L * np.array([1, 2, 3, 4.0])]
    kmax = 0.2 / (4 * L)
    box = [(-0.6, 0.6), (-kmax, kmax), (-6 * L, 6 * L), (-3.0, 3.0)]
    return surfs, kmax, box


def tilted(L, alpha_deg=45.0, axis="z"):
    a = np.radians(alpha_deg); surfs = []
    for i, x in enumerate(L * np.array([1, 2, 3, 4.0])):
        sgn = 1 if i % 2 == 0 else -1          # alternate the tilt direction
        if axis == "z":    # normal rotated in the transverse plane
            n = [np.cos(sgn * a), np.sin(sgn * a), 0.0]; eu = [-np.sin(sgn * a), np.cos(sgn * a), 0.0]
        else:              # normal rotated towards the beam axis (endcap-like)
            n = [np.cos(sgn * a), 0.0, np.sin(sgn * a)]; eu = [0.0, 1.0, 0.0]
        surfs.append(g.Plane([x, 0, 0], n, eu, -0.25 * x, 0.25 * x, -0.4 * x, 0.4 * x, S))
    kmax = 0.2 / (4 * L)
    box = [(-1.0, 1.0), (-kmax, kmax), (-8 * L, 8 * L), (-3.0, 3.0)]
    return surfs, kmax, box


def mixed(L):
    """two barrel cylinders, then two endcap disks (planes z = const)"""
    surfs = [g.Cylinder(R, 0.0, 0.25, 0.0, 3.0 * L, S) for R in (L, 2 * L)]
    for zd in (4 * L, 5 * L):
        surfs.append(g.Plane([0, 0, zd], [0, 0, 1], [1, 0, 0], 1.0 * L, 3.5 * L, -0.6 * L, 0.6 * L, S))
    kmax = 0.2 / (3 * L)
    box = [(-1.0, 1.0), (-kmax, kmax), (-6 * L, 6 * L), (0.0, 8.0)]
    return surfs, kmax, box


GEOMS = {"cylinders": cylinders, "tilted45z": lambda L: tilted(L, 45, "z"),
         "tilted60_endcaplike": lambda L: tilted(L, 60, "y"), "mixed": mixed}
