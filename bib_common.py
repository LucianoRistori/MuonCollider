"""
Common utilities for analyzing Muon Collider Beam-Induced Background (BIB)
n-tuples (HTAtree).

Detector cellID encoding of hit_id0 (32-bit unsigned int), as provided:
    "system:5,side:-2,layer:6,module:11,sensor:8"

  system:
    1 = VXD barrel   (vertex detector, barrel)
    2 = VXD endcap   (vertex detector, endcap disks)
    3 = IT barrel    (inner tracker, barrel)
    4 = IT endcap    (inner tracker, endcap disks)
    5 = OT barrel    (outer tracker, barrel)
    6 = OT endcap    (outer tracker, endcap disks)

  side:
    0   = barrel (no side distinction)
    1/3 = endcap +z / -z

  layer: layer/disk index within the subsystem (0-based)

Units (confirmed): x,y,z,u,v in mm; t in ns; px,py,pz in GeV.
"""

import hashlib
import os
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot

SYSTEM_NAMES = {
    1: "VXD barrel",
    2: "VXD endcap",
    3: "IT barrel",
    4: "IT endcap",
    5: "OT barrel",
    6: "OT endcap",
}

C_MM_PER_NS = 299.792458  # exact: c = 299,792,458 m/s
MUON_MASS_GEV = 0.1056583755  # PDG muon mass, GeV/c^2


def decode_id0(id0):
    """Unpack hit_id0 (uint32 array) into system, side, layer, module, sensor."""
    id0 = np.asarray(id0).astype(np.uint32)
    system = (id0 & 0x1F).astype(np.int32)
    side = ((id0 >> 5) & 0x3).astype(np.int32)
    layer = ((id0 >> 7) & 0x3F).astype(np.int32)
    module = ((id0 >> 13) & 0x7FF).astype(np.int32)
    sensor = ((id0 >> 24) & 0xFF).astype(np.int32)
    return system, side, layer, module, sensor


HIT_BRANCHES = [
    "hit_index", "hit_id0", "hit_x", "hit_y", "hit_z", "hit_t",
    "hit_u", "hit_v", "hit_du", "hit_dv",
    "hit_ux", "hit_uy", "hit_uz", "hit_vx", "hit_vy", "hit_vz",
    "hit_px", "hit_py", "hit_pz",
]
PRIMARY_BRANCHES = [
    "nRun", "nEvt", "part_pdg", "part_px", "part_py", "part_pz", "part_e",
    "part_vx", "part_vy", "part_vz", "part_t0", "part_q", "n_hit",
]


def load_hits(root_file, tree_name=None, entry=None):
    """
    Load hit-level branches from the HTAtree, flattened across events.

    entry=None (the default) loads and concatenates hits from ALL entries
    (events) in the tree - the natural behavior for a multi-event file:
    the merged plus+minus BIB file (see merge_bib_files.py, 2 entries), a
    signal/particle-gun sample (one entry per generated particle, e.g.
    the muon-gun file, 100k entries), or any other multi-event ntuple.
    For a single-entry file (the original plus/minus BIB files) this is
    identical to the old default of entry=0, so no existing single-event
    caller needs to change. Pass an integer `entry` to load just one
    specific event instead (old single-entry behavior).

    Returns a dict of numpy arrays (one row per hit, across however many
    events were loaded), including the decoded 'system', 'side', 'layer',
    'r' (transverse radius) fields, plus per-event primary-particle info
    under '_primary' (a dict of scalars if entry was given, or a dict of
    arrays - one value per event - if entry=None) and the number of
    events actually loaded under '_n_events'.
    """
    f = uproot.open(root_file)
    if tree_name is None:
        # pick the highest ROOT "cycle" of HTAtree, e.g. HTAtree;2
        candidates = sorted(
            [k for k in f.keys() if k.split(";")[0] == "HTAtree"],
            key=lambda k: int(k.split(";")[1]),
        )
        tree_name = candidates[-1]
    tree = f[tree_name]

    if entry is not None:
        arrs = tree.arrays(library="np")
        hits = {b[4:]: arrs[b][entry] for b in HIT_BRANCHES}
        primary = {k[5:] if k.startswith("part_") else k: arrs[k][entry]
                   for k in PRIMARY_BRANCHES}
        n_events = 1
        hits["event_id"] = np.zeros(len(hits["x"]), dtype=np.int64)
    else:
        import awkward as ak
        ak_hits = tree.arrays(HIT_BRANCHES, library="ak")
        hits = {b[4:]: ak.to_numpy(ak.flatten(ak_hits[b], axis=1)) for b in HIT_BRANCHES}
        np_primary = tree.arrays(PRIMARY_BRANCHES, library="np")
        primary = {k[5:] if k.startswith("part_") else k: np_primary[k]
                   for k in PRIMARY_BRANCHES}
        n_events = tree.num_entries
        # Per-hit event index (which event/track each hit came from) - needed
        # to group hits back into per-track counts, e.g. for a tracking-
        # efficiency calculation (see bib_common.apply_cuts / track_efficiency.py).
        # ak.flatten preserves event order (event 0's hits, then event 1's,
        # ...), so repeating each event index by its own hit count lines up
        # exactly with the flattened hit arrays above.
        n_hits_per_event = ak.to_numpy(ak.num(ak_hits[HIT_BRANCHES[0]], axis=1))
        hits["event_id"] = np.repeat(np.arange(n_events, dtype=np.int64), n_hits_per_event)

    system, side, layer, module, sensor = decode_id0(hits["id0"])
    hits["system"] = system
    hits["side"] = side
    hits["layer"] = layer
    hits["module"] = module
    hits["sensor"] = sensor
    hits["r"] = np.sqrt(hits["x"] ** 2 + hits["y"] ** 2)

    hits["_primary"] = primary
    hits["_tree_name"] = tree_name
    hits["_n_events"] = n_events
    return hits


def add_incidence_angles(hits):
    """
    Compute the angle of incidence of each hit's track (from its global
    momentum hit_px/py/pz) relative to the local sensor-plane normal,
    decomposed into a "longitudinal" and a "transverse" component.

    Local sensor normal: n_hat = normalize(u x v), where u, v are the
    module's own local axis unit vectors (hit_ux/uy/uz, hit_vx/vy/vz),
    oriented to point outward (away from the beam axis / interaction
    point: n_hat . hit_position > 0).

    Local cylindrical frame at the hit position (x, y, z), r = sqrt(x^2+y^2):
        rho_hat = (x, y, 0) / r    - radially outward, in the X-Y plane
        phi_hat = (-y, x, 0) / r   - azimuthal, in the X-Y plane
        z_hat   = (0, 0, 1)

    Two incidence-angle components, both measured relative to n_hat:
      - Longitudinal (theta_long_deg): signed angle, taken *within the
        meridian plane* (the plane containing the Z axis and the hit
        position, spanned by rho_hat and z_hat), between the track
        direction and the normal:
            theta_long = atan2(p_rho*n_z - p_z*n_rho, p_rho*n_rho + p_z*n_z)
        Well-defined unless the normal is purely azimuthal (n_rho, n_z
        both ~0), which never happens for a realistic tracker module.
      - Transverse (theta_trans_deg): signed angle representing the
        track's azimuthal (X-Y plane) tilt relative to the normal:
            theta_trans = atan2(p . phi_hat, p_hat . n_hat)

    Real tracker layers are built with (near-exact) rotational symmetry
    about the Z axis, so the sensor normal lies very nearly *in* the
    meridian plane itself (n . phi_hat =~ 0); under that condition
    tan(theta_long)^2 + tan(theta_trans)^2 =~ tan(theta_full)^2, where
    theta_full = arccos(p_hat . n_hat) is the plain 3D incidence angle
    (also stored, as theta_full_deg, as a sanity cross-check). The
    per-hit out-of-meridian-plane tilt of the normal itself is stored
    as n_tilt_deg = arcsin(|n . phi_hat|) - it should be small; where
    it isn't, the two components no longer cleanly add up to the full
    angle for that hit.

    Adds to `hits` (all angles in degrees): n_x/n_y/n_z (outward unit
    normal), theta_long_deg, theta_trans_deg, theta_full_deg,
    n_tilt_deg, z_axis_intercept_mm (the meridian-plane extrapolation
    of the track back to the Z axis) and psi_transverse_deg /
    inv_radius_per_mm (the transverse-plane circle-through-the-origin
    curvature) - see the comments just above where each is computed.
    Mutates and returns `hits`.
    """
    x, y, z = hits["x"], hits["y"], hits["z"]
    r = hits["r"]
    ux, uy, uz = hits["ux"], hits["uy"], hits["uz"]
    vx, vy, vz = hits["vx"], hits["vy"], hits["vz"]
    px, py, pz = hits["px"], hits["py"], hits["pz"]

    # local sensor normal = u x v, normalized, oriented outward
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    nnorm = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
    nx, ny, nz = nx / nnorm, ny / nnorm, nz / nnorm
    inward = (nx * x + ny * y + nz * z) < 0
    nx = np.where(inward, -nx, nx)
    ny = np.where(inward, -ny, ny)
    nz = np.where(inward, -nz, nz)

    # local cylindrical frame at the hit position
    rho_x, rho_y = x / r, y / r
    phi_x, phi_y = -y / r, x / r

    p_norm = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    p_rho = (px * rho_x + py * rho_y) / p_norm
    p_phi = (px * phi_x + py * phi_y) / p_norm
    p_z = pz / p_norm

    n_rho = nx * rho_x + ny * rho_y
    n_phi = nx * phi_x + ny * phi_y
    n_z = nz

    p_dot_n = p_rho * n_rho + p_phi * n_phi + p_z * n_z  # = cos(theta_full)

    theta_long = np.degrees(np.arctan2(
        p_rho * n_z - p_z * n_rho, p_rho * n_rho + p_z * n_z
    ))
    theta_trans = np.degrees(np.arctan2(p_phi, p_dot_n))
    theta_full = np.degrees(np.arccos(np.clip(p_dot_n, -1.0, 1.0)))
    n_tilt = np.degrees(np.arcsin(np.clip(np.abs(n_phi), 0.0, 1.0)))

    # Companion to theta_long: extrapolate the track's meridian-plane
    # (rho, z) trajectory - a straight line through (r, z) with local
    # slope (p_rho, p_z) - back to where it crosses the Z axis (rho=0):
    #     z_axis_intercept = z - r * (p_z / p_rho)
    # This is the "virtual" z origin of the hit's local track direction,
    # ignoring curvature and the transverse (phi) component entirely -
    # a purely kinematic quantity from the hit's own (r, z) position and
    # its meridian-plane momentum direction (no dependence on the sensor
    # normal). Undefined (+/-inf) for tracks running along z within the
    # meridian plane (p_rho ~ 0, i.e. theta_long ~ +/-90 deg relative to
    # rho_hat) - those are left as NaN rather than actual infinities.
    with np.errstate(divide="ignore", invalid="ignore"):
        z_axis_intercept = np.where(
            np.abs(p_rho) > 1e-9, z - r * (p_z / p_rho), np.nan
        )

    # Companion to theta_trans: the curvature 1/R of the circle, in the
    # transverse (X-Y) plane, that passes through both the origin (the
    # beam axis) and the hit, and is tangent there to the track's
    # transverse momentum direction - i.e. the circular arc a particle
    # produced at the origin would have to follow to reach this hit
    # going in this transverse direction. This is the standard
    # single-hit "sagitta" curvature estimate used to infer a track's
    # radius of curvature (and hence, given B, its transverse momentum -
    # not computed here, but R[m] and B[T] would give
    # pT[GeV/c] = 0.3 * B * R).
    #
    # Let psi = angle between the transverse momentum direction and the
    # radial direction rho_hat at the hit (rho_hat is exactly the
    # direction of the chord from the origin to the hit, so psi is the
    # standard tangent-chord angle). The tangent-chord relation for a
    # circle gives chord length r = 2*R*sin(psi), so:
    #     1/R = 2*sin(psi) / r
    # Bounded (|1/R| <= 2/r), unlike z_axis_intercept, since sin is
    # bounded and every hit has r > 0 - no exclusion needed. Sign
    # follows the sign of psi (i.e. which way the track curves).
    psi = np.arctan2(p_phi, p_rho)
    inv_radius_per_mm = 2.0 * np.sin(psi) / r

    hits["n_x"], hits["n_y"], hits["n_z"] = nx, ny, nz
    hits["theta_long_deg"] = theta_long
    hits["theta_trans_deg"] = theta_trans
    hits["theta_full_deg"] = theta_full
    hits["n_tilt_deg"] = n_tilt
    hits["z_axis_intercept_mm"] = z_axis_intercept
    hits["psi_transverse_deg"] = np.degrees(psi)
    hits["inv_radius_per_mm"] = inv_radius_per_mm
    return hits


def add_time_of_flight(hits, B_FIELD_T=5.0, mass_gev=MUON_MASS_GEV):
    """
    For each hit, compute the time a muon would take to travel from the
    origin (the IP, at the transverse-plane point rho=0) to the hit
    position, following the actual curved (helical) trajectory implied
    by the hit's own transverse curvature - and subtract that expected
    time-of-flight from the measured hit time. Requires
    add_incidence_angles(hits) to have been called first (uses
    psi_transverse_deg and inv_radius_per_mm).

    For a genuine signal muon produced at the IP at t=0 and travelling
    at close to c, t_corrected should cluster tightly around 0,
    regardless of how far out the hit is - unlike the raw hit time,
    which grows with distance from the IP. BIB hits, which don't
    necessarily originate at the IP at t=0 (secondaries from
    downstream interactions, slower non-relativistic particles, etc.),
    are not expected to peak at zero the same way - so t_corrected is a
    potential BIB-vs-signal discriminant in its own right.

    Geometry (two separate pieces):

    1. Transverse arc length from the origin to the hit, along the
       actual curved (circular, in the transverse plane) trajectory -
       NOT the straight-line chord. Reuses the same circle already
       used for inv_radius_per_mm: it passes through the origin and
       the hit, tangent there to the transverse momentum direction,
       and subtends angle 2*psi at its center (psi = psi_transverse_deg,
       the tangent-chord angle already computed by
       add_incidence_angles - see its docstring). Chord length is r
       (the hit's transverse radius) = 2*R*sin(psi), so:
           s_transverse = R * (2*psi) = r * psi / sin(psi)
       This form doesn't need R (=1/inv_radius_per_mm) explicitly and
       is numerically stable as psi -> 0 (straight track), where it
       correctly limits to s_transverse -> r (psi/sin(psi) -> 1).

    2. Total (3D) path length along the helix: a helical trajectory in
       a solenoidal (z-directed) field has constant pitch, since B only
       bends the transverse motion - pz and the transverse momentum
       magnitude are each separately conserved along the true path
       (ignoring small energy loss/scattering). So the ratio of
       transverse path length to total path length is fixed along the
       whole trajectory, equal to sin(theta_polar) = pT/p - using the
       hit's own (raw, measured) local momentum direction, which is
       the actual known pitch angle, independent of how pT's magnitude
       is estimated below:
           s_total = s_transverse / sin(theta_polar)

    Momentum/speed estimate: per the request, the transverse momentum
    magnitude used for the speed (beta) calculation comes from the
    curvature (inv_radius_per_mm -> pT = 0.3*B*R, same as the pT axis
    on inv_radius_per_subsystem.png), not directly from the raw hit
    momentum - this is what a real curvature-based track-pT
    measurement would give. Combined with the same (raw-direction)
    polar angle as above to get a total momentum estimate:
        pT_curv = 0.3 * B_FIELD_T / |inv_radius_per_m|
        p_estimate = pT_curv / sin(theta_polar)
        beta = p_estimate / sqrt(p_estimate^2 + mass_gev^2)
    A fixed particle mass (muon, by default) is assumed for every hit,
    regardless of what actually produced it - that's the point: this
    tests which hits are consistent with an on-time muon from the IP,
    whatever they actually are.

    Expected time-of-flight (ns, since c is in mm/ns and lengths in mm):
        tof_expected_ns = s_total / (beta * C_MM_PER_NS)
        t_corrected_ns = hit_t - tof_expected_ns

    Adds to `hits`: path_length_mm (s_total), p_estimate_gev, beta_estimate,
    tof_expected_ns, t_corrected_ns. Mutates and returns `hits`.
    """
    if "psi_transverse_deg" not in hits or "inv_radius_per_mm" not in hits:
        raise ValueError("add_time_of_flight requires add_incidence_angles(hits) first")

    r = hits["r"]
    px, py, pz = hits["px"], hits["py"], hits["pz"]
    t = hits["t"]
    psi = np.radians(hits["psi_transverse_deg"])
    inv_radius_per_mm = hits["inv_radius_per_mm"]

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(np.abs(psi) > 1e-9, psi / np.sin(psi), 1.0)
    s_transverse = r * ratio

    pT_raw = np.sqrt(px ** 2 + py ** 2)
    p_raw = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        sin_theta_polar = np.where(p_raw > 0, pT_raw / p_raw, np.nan)
        s_total = np.where(sin_theta_polar > 1e-12, s_transverse / sin_theta_polar, np.nan)

    GEV_PER_INV_M = 0.3 * B_FIELD_T
    inv_radius_per_m = np.abs(inv_radius_per_mm) * 1000.0
    with np.errstate(divide="ignore", invalid="ignore"):
        pT_curv = np.where(inv_radius_per_m > 0, GEV_PER_INV_M / inv_radius_per_m, np.inf)
        p_estimate = np.where(sin_theta_polar > 1e-12, pT_curv / sin_theta_polar, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        beta = p_estimate / np.sqrt(p_estimate ** 2 + mass_gev ** 2)
        tof_expected_ns = s_total / (beta * C_MM_PER_NS)

    # Signed version of the curvature-derived transverse momentum, sign
    # following the sign of the curvature itself (same convention as the
    # pT axis relabeling on inv_radius_per_subsystem.png: "sign of pT
    # follows the sign of the curvature") - this is the "momentum
    # estimated from the transverse angle" quantity used for the
    # momentum selection cut (see cuts_config.txt / load_cuts /
    # apply_cuts below).
    with np.errstate(invalid="ignore"):
        pT_curv_signed = np.sign(inv_radius_per_mm) * pT_curv

    hits["path_length_mm"] = s_total
    hits["p_estimate_gev"] = p_estimate
    hits["beta_estimate"] = beta
    hits["tof_expected_ns"] = tof_expected_ns
    hits["t_corrected_ns"] = t - tof_expected_ns
    hits["pT_curv_gev"] = pT_curv_signed
    return hits


# ---------------------------------------------------------------------------
# Selection cuts: an editable text config (cuts_config.txt) defines a
# symmetric acceptance window around zero for each of three quantities that
# should be consistent with a genuine signal muon from the IP:
#   t_corrected_ns      - time-of-flight-corrected hit time (see
#                          add_time_of_flight); expected ~0 for signal.
#                          Simple window: accept |t_corrected_ns| <= halfwidth.
#   z_axis_intercept_mm - meridian-plane track extrapolated back to the Z
#                          axis (see add_incidence_angles); expected ~0
#                          (the IP) for signal. Simple window: accept
#                          |z_axis_intercept_mm| <= halfwidth.
#   momentum_gev         - transverse momentum estimated from the hit's own
#                          curvature. This one is NOT a symmetric window on
#                          pT itself: pT = 0.3*B/|1/R| diverges as 1/R -> 0,
#                          so a perfectly straight track (R = Infinity, the
#                          best-reconstructed, highest-momentum case) has
#                          pT = +/-Infinity depending which side of R=0 you
#                          approach from - a genuine discontinuity that no
#                          finite window in pT-space can straddle. Instead
#                          the cut is applied directly on the (bounded,
#                          continuous) curvature inv_radius_per_mm: accept
#                          |1/R| <= (0.3*B_FIELD_T/1000) / halfwidth_gev,
#                          i.e. "reconstructed |pT| >= halfwidth_gev, with
#                          R = Infinity always accepted since 1/R = 0 sits
#                          at the exact center of that window".
# Each hit must pass every *enabled* cut to be accepted. A cut whose value
# is NaN (e.g. an undefined z_axis_intercept) fails that cut, same as it's
# already excluded/reported separately elsewhere in this codebase.
# ---------------------------------------------------------------------------

CUT_NAMES = ("t_corrected_ns", "z_axis_intercept_mm", "momentum_gev")

# Default half-range of each cut variable's *zoomed* histogram (same units
# as that cut's halfwidth), used when cuts_config.txt doesn't set its own
# zoom_halfwidth - see load_cuts(). These match the values the zoomed
# plots used before zoom_halfwidth became configurable (Z0_ZOOM_RANGE_MM,
# PT_ZOOM_RANGE_GEV, TC_ZOOM_RANGE_NS in incidence_angle_plots.py /
# time_of_flight_plots.py), kept here too as the fallback for scripts
# (currently n1_cut_plots.py) that size their zoom window from the cuts
# config so it can track a cut's own halfwidth as that's changed.
ZOOM_HALFWIDTH_DEFAULTS = {
    "t_corrected_ns": 2.0,          # ns
    "z_axis_intercept_mm": 100.0,   # mm
    "momentum_gev": 1.0,            # GeV/c
}

# system IDs (see decode_id0 / SYSTEM_NAMES above) that make up the vertex
# detector (VXD barrel + endcap) - used by the [track] section's
# exclude_vertex_hits option, see load_track_params().
VERTEX_SYSTEM_IDS = {1, 2}


def load_cuts(path):
    """
    Load a selection-cuts config file (see cuts_config.txt for the
    canonical example/format). Uses the standard library configparser, so
    the file is plain INI: one section per cut (t_corrected_ns,
    z_axis_intercept_mm, momentum_gev), each with a `halfwidth` and an
    `enabled` key (t_corrected_ns / z_axis_intercept_mm also take an
    optional `center`, default 0.0; momentum_gev has no center - see the
    module-level comment above for why it's inherently centered on
    1/R = 0). Each cut also takes an optional `zoom_halfwidth` - the
    half-range (same units as `halfwidth`) of that variable's zoomed
    histogram, used by n1_cut_plots.py; defaults to
    ZOOM_HALFWIDTH_DEFAULTS[name] if not set. Lines starting with '#' or
    ';' are comments.

    Returns a dict {cut_name: {...}}. A cut with no section in the file is
    treated as enabled=False (no cut applied), with its zoom_halfwidth
    still defaulted from ZOOM_HALFWIDTH_DEFAULTS.
    """
    import configparser

    cp = configparser.ConfigParser()
    read_ok = cp.read(path)
    if not read_ok:
        raise FileNotFoundError(f"cuts config not found: {path}")

    cuts = {}
    for name in CUT_NAMES:
        if cp.has_section(name):
            sec = cp[name]
            cuts[name] = {
                "center": sec.getfloat("center", 0.0),
                "halfwidth": sec.getfloat("halfwidth", fallback=float("inf")),
                "enabled": sec.getboolean("enabled", fallback=True),
                "zoom_halfwidth": sec.getfloat(
                    "zoom_halfwidth", fallback=ZOOM_HALFWIDTH_DEFAULTS[name]),
            }
        else:
            cuts[name] = {"center": 0.0, "halfwidth": float("inf"), "enabled": False,
                           "zoom_halfwidth": ZOOM_HALFWIDTH_DEFAULTS[name]}
    return cuts


def load_track_params(path):
    """
    Load track-level parameters from the same cuts config file (see
    cuts_config.txt, [track] section) used by load_cuts():
      - min_hits_found: the minimum number of hits surviving the
        selection cuts (see apply_cuts()) for a simulated track (all the
        hits from one event/one generated particle, grouped by
        hits["event_id"]) to be counted as "found" - see
        track_efficiency.py. Default 5.
      - exclude_vertex_hits: if true, hits in the vertex detector
        (VXD barrel/endcap, system in VERTEX_SYSTEM_IDS) are excluded
        from that per-track surviving-hit count, even if they pass the
        selection cuts - i.e. "found" then means >=min_hits_found
        surviving hits *outside* the vertex detector. Default false
        (vertex hits count same as any other).
    Returns {"min_hits_found": int, "exclude_vertex_hits": bool}.
    """
    import configparser

    cp = configparser.ConfigParser()
    read_ok = cp.read(path)
    if not read_ok:
        raise FileNotFoundError(f"cuts config not found: {path}")
    min_hits_found = cp.getint("track", "min_hits_found", fallback=5) \
        if cp.has_section("track") else 5
    exclude_vertex_hits = cp.getboolean("track", "exclude_vertex_hits", fallback=False) \
        if cp.has_section("track") else False
    return {"min_hits_found": min_hits_found, "exclude_vertex_hits": exclude_vertex_hits}


def apply_cuts(hits, cuts, B_FIELD_T=5.0):
    """
    Apply a cuts dict (from load_cuts()) to `hits`. `hits` must already
    have had add_incidence_angles() and add_time_of_flight() called on it.
    B_FIELD_T must match whatever was used for add_time_of_flight (and the
    pT axis on the angle plots), since the momentum_gev cut converts its
    GeV/c halfwidth to a curvature threshold using it.

    Returns (combined_mask, per_cut_masks) - see CUT_NAMES / the
    module-level comment above for what each cut does. per_cut_masks maps
    cut name -> its own accept mask on its own (before combining), useful
    for a cutflow table. A disabled cut contributes an all-True mask (it's
    simply skipped); NaN values in a cut's underlying quantity fail that
    cut even when it's enabled.
    """
    n = len(hits["t_corrected_ns"])
    combined = np.ones(n, dtype=bool)
    per_cut = {}

    for name in ("t_corrected_ns", "z_axis_intercept_mm"):
        spec = cuts[name]
        if not spec["enabled"]:
            per_cut[name] = np.ones(n, dtype=bool)
            continue
        v = hits[name]
        lo, hi = spec["center"] - spec["halfwidth"], spec["center"] + spec["halfwidth"]
        with np.errstate(invalid="ignore"):
            mask = (v >= lo) & (v <= hi)
        per_cut[name] = mask
        combined &= mask

    spec = cuts["momentum_gev"]
    if not spec["enabled"]:
        per_cut["momentum_gev"] = np.ones(n, dtype=bool)
    else:
        GEV_PER_INV_M = 0.3 * B_FIELD_T
        GEV_PER_INV_MM = GEV_PER_INV_M / 1000.0  # inv_radius_per_mm is in 1/mm
        inv_radius_threshold = GEV_PER_INV_MM / spec["halfwidth"]
        v = hits["inv_radius_per_mm"]
        with np.errstate(invalid="ignore"):
            mask = np.abs(v) <= inv_radius_threshold
        per_cut["momentum_gev"] = mask
        combined &= mask

    return combined, per_cut


def mask_hits(hits, mask):
    """
    Return a new hits dict with every per-hit array field filtered by
    `mask` (a boolean array the same length as the hits), leaving
    per-file metadata fields (_primary, _tree_name, _n_events) unchanged.
    Useful for re-running region_table/add_peak_density/etc. on only the
    hits that survive a selection (see apply_cuts()).
    """
    n = len(mask)
    out = {}
    for k, v in hits.items():
        if isinstance(v, np.ndarray) and v.shape[:1] == (n,):
            out[k] = v[mask]
        else:
            out[k] = v
    return out


def region_table(hits):
    """
    Build a per-(system, side, layer) summary table: hit counts, basic
    geometric extent, an estimated sensitive area, and hit density
    (hits/mm^2). Returns a list of dict rows, sorted by (system, side,
    layer).

    Area estimate: within a (system, side, layer) region, sensor local
    coordinates (u, v) are pooled across all hit modules to get a single
    "typical module footprint" (u_max-u_min) x (v_max-v_min) - this is
    robust because module design is uniform within a layer/disk, so
    pooling statistics across modules fills in the footprint even where
    a single module has too few hits to trace its own edges. That
    footprint is multiplied by the number of distinct modules that
    registered >=1 hit in this file to get the region's sensitive area.

    Caveat: modules that never got hit in this event are invisible to
    us, so n_modules (and hence area and density) is a lower bound on
    the true installed area, more so for low-occupancy regions. Ratios
    between cuts/samples computed with the same recipe are still valid
    since the normalization is applied consistently.
    """
    system, side, layer, r, z, u, v, module = (
        hits["system"], hits["side"], hits["layer"], hits["r"], hits["z"],
        hits["u"], hits["v"], hits["module"],
    )
    combo = np.stack([system, side, layer], axis=1)
    uniq, inv, counts = np.unique(combo, axis=0, return_inverse=True, return_counts=True)
    order = np.lexsort((uniq[:, 2], uniq[:, 1], uniq[:, 0]))

    rows = []
    for i in order:
        mask = inv == i
        sysv, sidev, layv = uniq[i]
        n_modules = int(np.unique(module[mask]).size)
        u_lo, u_hi = float(u[mask].min()), float(u[mask].max())
        v_lo, v_hi = float(v[mask].min()), float(v[mask].max())
        module_area = (u_hi - u_lo) * (v_hi - v_lo)
        area = module_area * n_modules
        n_hits = int(counts[i])
        rows.append({
            "system": int(sysv),
            "system_name": SYSTEM_NAMES.get(int(sysv), f"unknown({sysv})"),
            "side": int(sidev),
            "layer": int(layv),
            "n_hits": n_hits,
            "n_modules_hit": n_modules,
            "module_area_mm2": module_area,
            "area_mm2": area,
            "density_hits_per_mm2": n_hits / area if area > 0 else float("nan"),
            "mean_r": float(r[mask].mean()),
            "min_r": float(r[mask].min()),
            "max_r": float(r[mask].max()),
            "mean_z": float(z[mask].mean()),
            "min_z": float(z[mask].min()),
            "max_z": float(z[mask].max()),
        })
    return rows


def add_peak_density(hits, rows, bin_size_mm=None, percentile=99.0,
                      target_hits_per_bin=20.0, min_bin_mm=1.0, max_bin_mm=200.0):
    """
    Add local ("peak") hit density to each region_table() row, alongside
    its already-computed mean density.

    Unlike the mean density (which uses the *local* sensor u,v footprint,
    since it needs the module's own area), the peak density needs real
    spatial position, so it bins the region's actual *global* hit
    positions onto a physical surface:
      - barrel regions (side==0): unrolled (arc_length, z), where
        arc_length = mean_r * phi - a flat approximation, valid because
        each barrel layer is a thin shell at ~constant r.
      - endcap regions (side in {1,3}): (x, y) directly - the disk is
        already flat, so this is exact.
    Each surface is binned into bin_size_mm x bin_size_mm cells; local
    density = hits in a cell / cell area.

    With only one simulated event, a single stray hit in a tiny bin is
    indistinguishable from a real hot spot, so the reported "peak" is
    the `percentile`-th percentile of nonzero bin densities (robust to
    single-bin noise), not the strict maximum - though the strict max is
    also recorded for reference. bin_size_mm therefore isn't a detail:
    it *is* part of the definition of "peak density" here, and is
    recorded on every row so the number stays interpretable.

    Bin size: mean hit density varies by 5+ orders of magnitude across
    regions (dense inner VXD layers vs. sparse outer OT endcap disks),
    so a single fixed bin size is either too fine (peak collapses to
    small-integer bin counts dominated by single-hit Poisson noise in
    sparse regions - the bin-occupancy quantization problem) or too
    coarse (throws away real structure in dense regions). By default
    (bin_size_mm=None) the bin size is therefore chosen per region,
    from that region's own mean density, to target roughly
    `target_hits_per_bin` hits per occupied bin on average:
        bin_size = sqrt(target_hits_per_bin / mean_density)
    clipped to [min_bin_mm, max_bin_mm]. Pass an explicit bin_size_mm
    to force the old fixed-bin-size behavior instead (e.g. for a direct
    side-by-side comparison across regions at one physical scale).

    Mutates and returns `rows`.
    """
    system, side, layer = hits["system"], hits["side"], hits["layer"]
    x, y, z = hits["x"], hits["y"], hits["z"]
    phi = np.arctan2(y, x)

    for row in rows:
        mask = (
            (system == row["system"]) & (side == row["side"]) & (layer == row["layer"])
        )
        if row["side"] == 0:  # barrel: unroll to (arc length, z)
            coord1 = row["mean_r"] * phi[mask]
            coord2 = z[mask]
        else:  # endcap: flat disk, use (x, y) directly
            coord1 = x[mask]
            coord2 = y[mask]

        if bin_size_mm is not None:
            row_bin_mm = bin_size_mm
        else:
            mean_density = row.get("density_hits_per_mm2", float("nan"))
            if not mean_density or not np.isfinite(mean_density) or mean_density <= 0:
                row_bin_mm = max_bin_mm
            else:
                row_bin_mm = float(np.clip(
                    np.sqrt(target_hits_per_bin / mean_density), min_bin_mm, max_bin_mm
                ))
        bin_area = row_bin_mm ** 2

        b1 = np.floor(coord1 / row_bin_mm).astype(np.int64)
        b2 = np.floor(coord2 / row_bin_mm).astype(np.int64)
        _, counts = np.unique(np.stack([b1, b2], axis=1), axis=0, return_counts=True)
        bin_densities = counts / bin_area

        row["peak_bin_size_mm"] = row_bin_mm
        row["peak_n_bins_occupied"] = int(counts.size)
        row["peak_mean_hits_per_bin"] = float(counts.mean()) if counts.size else 0.0
        row[f"peak_density_p{percentile:.0f}_hits_per_mm2"] = float(
            np.percentile(bin_densities, percentile)
        )
        row["peak_density_max_hits_per_mm2"] = float(bin_densities.max())
    return rows


def _display_path(path):
    """
    Rewrite a filesystem path for human-facing display (logs, print
    statements) using the BIB_PATH_REMAP environment variable, if set.

    This exists only because scripts are sometimes run through a
    remote-device bridge, where the real path on the user's machine
    (e.g. /Users/luciano/Dropbox/Documents/MuonColliderSimulation) is
    mounted at a different, session-internal path
    (e.g. /sessions/<id>/mnt/MuonColliderSimulation) - so os.path
    calls inside the script see the internal path, not the one the
    user actually has on disk. When invoking through such a bridge,
    set BIB_PATH_REMAP="internal1=>real1;internal2=>real2" so logs
    show the real, human-meaningful path instead. Run directly on the
    user's own machine (a normal terminal), this is a no-op.
    """
    remap = os.environ.get("BIB_PATH_REMAP", "")
    path = str(path)
    pairs = []
    for pair in remap.split(";"):
        if "=>" not in pair:
            continue
        old, new = pair.split("=>", 1)
        old, new = old.strip(), new.strip()
        if old:
            pairs.append((old, new))
    # longest prefix first, so e.g. ".../MuonCollider" doesn't shadow
    # ".../MuonColliderSimulation"
    for old, new in sorted(pairs, key=lambda p: -len(p[0])):
        if path.startswith(old):
            return new + path[len(old):]
    return path


def prepare_output_dir(outdir):
    """
    Ensure `outdir` exists and is empty, without ever silently
    overwriting a previous run's results: if `outdir` already exists
    (e.g. from a previous run with the same name), it is moved whole
    into an "OldResults" folder next to it, tagged with the move
    timestamp, before a fresh empty `outdir` is created.

    e.g. .../results/step1_basic_plots/  (existing)
      -> .../results/OldResults/step1_basic_plots_20260921-150200/
      then a new, empty .../results/step1_basic_plots/ is created.
    """
    import shutil

    outdir = Path(outdir)
    if outdir.exists() and any(outdir.iterdir()):
        old_results = outdir.parent / "OldResults"
        old_results.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        dest = old_results / f"{outdir.name}_{stamp}"
        shutil.move(str(outdir), str(dest))
        print(f"Archived previous results: {_display_path(outdir)} -> {_display_path(dest)}")
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def _sha256_of_file(path, chunk_size=8 * 1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def write_logfile(log_path, root_file, hits, sys_rows, script_name, outdir,
                   hash_input=True, area_source="hit-inferred"):
    """
    Write a plain-text run log documenting how a set of results was
    produced: when, from which input file (with size/checksum for
    provenance), what the file contained, and a summary of what came
    out (hit counts and densities per subsystem).

    area_source: "geometry" if mean-density areas came from the true
    detector geometry (geometry.build_area_lookup), else "hit-inferred".
    """
    n_hit = len(hits["x"])
    st = os.stat(root_file)
    input_size_mb = st.st_size / (1024 ** 2)
    input_mtime = datetime.fromtimestamp(st.st_mtime, tz=timezone.utc)
    sha256 = _sha256_of_file(root_file) if hash_input else None

    lines = []
    lines.append("=" * 70)
    lines.append("Muon Collider BIB analysis - run log")
    lines.append("=" * 70)
    lines.append(f"Run timestamp (UTC)   : {datetime.now(timezone.utc).isoformat()}")
    lines.append(f"Script                : {script_name}")
    lines.append(f"Output directory      : {_display_path(outdir)}")
    lines.append("")
    lines.append("-- Input file --")
    lines.append(f"Path                  : {_display_path(os.path.abspath(root_file))}")
    lines.append(f"Size                  : {input_size_mb:.1f} MB")
    lines.append(f"Last modified (UTC)   : {input_mtime.isoformat()}")
    if sha256:
        lines.append(f"SHA-256               : {sha256}")
    lines.append(f"ROOT tree read        : {hits['_tree_name']}")
    lines.append("")
    lines.append("-- Hits read out --")
    lines.append(f"Total hits (n_hit)    : {n_hit:,}")
    lines.append(f"Detector regions      : {len(sys_rows)} subsystems, "
                 f"see region_summary.csv for full (system,side,layer) breakdown")
    lines.append("")
    peak_p_key = next(
        (k for k in (sys_rows[0] if sys_rows else {}) if k.startswith("peak_density_p")),
        None,
    )
    peak_percentile_str = None
    if peak_p_key:
        m = re.match(r"peak_density_p(\d+(?:\.\d+)?)_", peak_p_key)
        peak_percentile_str = m.group(1) if m else "?"
    area_note = ("true geometry" if area_source == "geometry"
                 else "hit-inferred, approximate")
    if peak_p_key:
        lines.append("-- Hit density summary (hits/mm^2), per subsystem --")
        lines.append(f"   (mean = total hits / total sensitive area [{area_note}]; "
                      f"peak = hottest {peak_percentile_str}th-percentile bin found in "
                      f"that subsystem; bin size is chosen per (system,side,layer) "
                      f"region to target ~20 hits/bin - see peak_bin_size_mm and "
                      f"peak_mean_hits_per_bin in region_summary.csv - so it is not "
                      f"one fixed size across a subsystem, only for the hottest region)")
        for r in sys_rows:
            lines.append(f"  {r['system_name']:12s}  n_hits={r['n_hits']:>9,}  "
                         f"area={r['area_mm2']:>12,.1f} mm^2  "
                         f"mean={r['density_hits_per_mm2']:.4g}  "
                         f"peak={r[peak_p_key]:.4g}  "
                         f"peak_max={r['peak_density_max_hits_per_mm2']:.4g}  "
                         f"(bin={r.get('peak_bin_size_mm', 0):.3g}mm)")
    else:
        lines.append("-- Hit density summary (hits/mm^2), per subsystem --")
        for r in sys_rows:
            lines.append(f"  {r['system_name']:12s}  n_hits={r['n_hits']:>9,}  "
                         f"area={r['area_mm2']:>12,.1f} mm^2  "
                         f"density={r['density_hits_per_mm2']:.4g}")
    lines.append("")
    lines.append("-- Software environment --")
    lines.append(f"Python                : {sys.version.split()[0]}")
    lines.append(f"uproot                : {uproot.__version__}")
    lines.append(f"numpy                 : {np.__version__}")
    lines.append(f"Platform              : {platform.platform()}")
    lines.append("=" * 70)

    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return log_path


def subsystem_density_table(rows):
    """
    Aggregate region_table() rows up to one row per subsystem (system),
    summing hits and area across all (side, layer) regions, then taking
    the ratio - NOT the mean of per-region densities.

    If per-region peak density fields (from add_peak_density()) are
    present, the subsystem-level peak is taken as the max across its
    regions' peak values - the hottest spot found anywhere in that
    subsystem, not an average of hot spots.
    """
    peak_p_key = next(
        (k for k in (rows[0] if rows else {}) if k.startswith("peak_density_p")), None
    )
    by_sys = {}
    for row in rows:
        s = row["system"]
        d = by_sys.setdefault(s, {
            "system": s, "system_name": row["system_name"],
            "n_hits": 0, "area_mm2": 0.0,
        })
        d["n_hits"] += row["n_hits"]
        d["area_mm2"] += row["area_mm2"]
        if peak_p_key:
            if row[peak_p_key] >= d.get(peak_p_key, -1.0):
                # bin size used by the region that set the hottest-spot
                # value, since bin size now varies per region (adaptive
                # sizing - see add_peak_density) and a single subsystem
                # no longer has one bin size to report.
                d["peak_bin_size_mm"] = row.get("peak_bin_size_mm")
            d[peak_p_key] = max(d.get(peak_p_key, 0.0), row[peak_p_key])
            d["peak_density_max_hits_per_mm2"] = max(
                d.get("peak_density_max_hits_per_mm2", 0.0),
                row["peak_density_max_hits_per_mm2"],
            )
    out = []
    for s in sorted(by_sys):
        d = by_sys[s]
        d["density_hits_per_mm2"] = (
            d["n_hits"] / d["area_mm2"] if d["area_mm2"] > 0 else float("nan")
        )
        out.append(d)
    return out


# ---------------------------------------------------------------------------
# Measurement-error ("smearing") config, for studying the effect of assumed
# detector resolution on the analysis. The current input ROOT files carry
# zero measurement error (hit_du/hit_dv are exactly 0.0 for every hit, and
# hit_u/hit_v are exactly self-consistent with hit_x/y/z - verified via
# check_smearing.py), unlike earlier samples which had resolution baked in
# at generation time. Rather than regenerating ROOT files per hypothesis,
# smearing here is applied on load, controlled by an editable config (see
# smearing_config.txt), the same way cuts_config.txt controls the selection
# cuts - so different error sizes can be scanned with no regeneration.
SMEARING_SECTIONS = ["position", "time", "angle_long", "angle_trans"]


def load_smearing_config(path):
    """
    Load a measurement-smearing config file (see smearing_config.txt for
    the canonical example/format). Plain INI (configparser), one section
    per smeared quantity:
      - [position]: sigma_uv_mm - Gaussian sigma (mm) applied
        independently to the hit's local hit_u and hit_v (NOT global
        hit_x/hit_y - see apply_position_time_smearing() docstring for
        why local u,v is the physically correct quantity to smear).
      - [time]: sigma_t_ns - Gaussian sigma (ns) applied to raw hit_t.
      - [angle_long] / [angle_trans]: sigma_deg - Gaussian sigma
        (degrees) applied directly to the already-computed
        theta_long_deg / theta_trans_deg (bib_common.add_incidence_angles
        output) - representing angular-reconstruction uncertainty on top
        of position/time smearing, not a re-derivation from smeared
        positions.
      - [general]: seed (int, default 42) - RNG seed, so re-running with
        the same sigmas reproduces the same smeared sample.
    Each of the four quantity sections also takes `enabled` (default
    false). A missing section, or enabled=false, or sigma=0, all mean
    "leave this quantity unsmeared" - so the default config (all
    disabled) reproduces the current, exact input files unchanged.

    Returns a dict {section_name: {"sigma": float, "enabled": bool}} for
    each of SMEARING_SECTIONS, plus {"seed": int}.
    """
    import configparser

    cp = configparser.ConfigParser()
    read_ok = cp.read(path)
    if not read_ok:
        raise FileNotFoundError(f"smearing config not found: {path}")

    sigma_key = {
        "position": "sigma_uv_mm",
        "time": "sigma_t_ns",
        "angle_long": "sigma_deg",
        "angle_trans": "sigma_deg",
    }
    cfg = {}
    for name in SMEARING_SECTIONS:
        if cp.has_section(name):
            sec = cp[name]
            cfg[name] = {
                "sigma": sec.getfloat(sigma_key[name], fallback=0.0),
                "enabled": sec.getboolean("enabled", fallback=False),
            }
        else:
            cfg[name] = {"sigma": 0.0, "enabled": False}
    cfg["seed"] = cp.getint("general", "seed", fallback=42) if cp.has_section("general") else 42
    return cfg


def smearing_rng(smear_cfg):
    """Build the numpy Generator to use for a run, from a loaded smearing config's seed."""
    return np.random.default_rng(smear_cfg["seed"])


def apply_position_time_smearing(hits, smear_cfg, rng):
    """
    Apply Gaussian position and/or time smearing to `hits` in place, per
    the [position]/[time] sections of `smear_cfg` (see
    load_smearing_config()). No-op for a disabled/zero-sigma quantity.
    Call this right after load_hits(), before add_incidence_angles() /
    add_time_of_flight() - both are computed from hit_x/y/z, hit_t and so
    automatically reflect the smearing applied here.

    Position smearing acts on the hit's LOCAL sensor-plane coordinates
    hit_u/hit_v, not the global hit_x/hit_y/hit_z - this is the correct
    physical quantity: a real tracker module measures position within its
    own local readout plane (see the README's "hit_u/hit_v" and
    "hit_du/hit_dv" description), so resolution belongs there, not on the
    global frame (which would ignore module tilt/orientation, wrong for
    endcap disks and tilted barrel ladders alike - see
    add_incidence_angles()'s n_tilt_deg).

    After drawing delta_u, delta_v ~ N(0, sigma_uv_mm) independently per
    hit and adding them to hits["u"]/hits["v"], hits["x"]/["y"]/["z"] are
    updated by the SAME delta, expressed in the global frame via the
    module's own local axes (hits["ux"/"uy"/"uz"], ["vx"/"vy"/"vz"]):
        new_global = old_global + delta_u * u_axis + delta_v * v_axis
    This keeps hit_u/hit_v and hit_x/y/z mutually consistent (the same
    self-consistency check that confirmed the input files are currently
    unsmeared - see check_smearing.py - continues to pass after this
    smearing is applied), and only moves the hit within its sensor's
    plane (out-of-plane position is unaffected, as for a real planar
    sensor). hits["du"]/["dv"] are set to the sigma actually applied, so
    the assumed resolution is visible in the (smeared) hits dict too.
    hits["r"] = sqrt(x^2+y^2) is recomputed to stay consistent.

    Time smearing simply adds N(0, sigma_t_ns) to hits["t"].
    """
    n = len(hits["x"])

    # Memory note: this runs on samples up to ~17M hits in a
    # memory-constrained environment, so every step below is written to
    # avoid allocating extra full-length temporaries where an in-place
    # (+=) update or a reused scratch buffer will do, rather than the
    # more readable `hits["x"] = hits["x"] + a*b + c*d` form (which
    # would allocate 3-4 extra full arrays at once).
    pos_cfg = smear_cfg["position"]
    if pos_cfg["enabled"] and pos_cfg["sigma"] > 0:
        sigma = pos_cfg["sigma"]
        du = rng.normal(0.0, sigma, size=n)
        dv = rng.normal(0.0, sigma, size=n)
        hits["u"] += du
        hits["v"] += dv

        scratch = du * hits["ux"]
        scratch += dv * hits["vx"]
        hits["x"] += scratch

        scratch = du * hits["uy"]
        scratch += dv * hits["vy"]
        hits["y"] += scratch

        scratch = du * hits["uz"]
        scratch += dv * hits["vz"]
        hits["z"] += scratch
        del du, dv, scratch

        hits["r"] = np.hypot(hits["x"], hits["y"], out=hits["r"])
        hits["du"].fill(sigma)
        hits["dv"].fill(sigma)

    time_cfg = smear_cfg["time"]
    if time_cfg["enabled"] and time_cfg["sigma"] > 0:
        hits["t"] += rng.normal(0.0, time_cfg["sigma"], size=n)

    import gc
    gc.collect()
    return hits


def apply_angle_smearing(hits, smear_cfg, rng):
    """
    Apply Gaussian smearing directly to the already-computed
    theta_long_deg / theta_trans_deg (per the [angle_long]/[angle_trans]
    sections of `smear_cfg` - see load_smearing_config()), in place.
    Must be called AFTER add_incidence_angles(). No-op for a disabled/
    zero-sigma quantity.

    This represents angular-reconstruction uncertainty as an independent
    error source on top of position/time smearing (e.g. multiple
    scattering or track-fit angular resolution), rather than deriving it
    from smeared hit positions - it does NOT touch z_axis_intercept_mm,
    inv_radius_per_mm or psi_transverse_deg, which are computed
    independently from the raw momentum direction, not from
    theta_long_deg/theta_trans_deg.
    """
    n = len(hits["theta_long_deg"])

    long_cfg = smear_cfg["angle_long"]
    if long_cfg["enabled"] and long_cfg["sigma"] > 0:
        hits["theta_long_deg"] += rng.normal(0.0, long_cfg["sigma"], size=n)

    trans_cfg = smear_cfg["angle_trans"]
    if trans_cfg["enabled"] and trans_cfg["sigma"] > 0:
        hits["theta_trans_deg"] += rng.normal(0.0, trans_cfg["sigma"], size=n)

    return hits
