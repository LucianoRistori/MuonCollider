# CORRECTION (10 October 2026): the depth-view fit in this framework had r and z
# swapped (it fitted r = p + q z, with residuals in r instead of z), so the fake
# probabilities computed or used here are about 450 times too low for the helix.
# Kept as the record of the earlier paper draft; the corrected values come from
# fake_rate_framework/geometric_K.py and are used by fake_rate.py.
"""
angle_resolution_scan.py

Scans the assumed pointing-angle (incidence-angle reconstruction)
resolution sigma_angle from 0 deg (no smearing) to 5 deg, 11 points
equally spaced on a LINEAR scale (0.5 deg steps), applied equally to
both local angle components (angle_u, angle_v) and to ALL systems
(not just IT/OT barrel -- see note below). Position resolution and
timing resolution are held fixed at today's values throughout. The
per-hit momentum_gev (pT) selection cut is REINTRODUCED here, held
fixed at today's threshold VALUE (see REVISION2 note below) -- a prior
version of this script removed it entirely, then put the ENDCAP-cut
bug (see REVISION note) back in play; this version has both the
endcap fix AND the pT cut, so its numbers are directly comparable to
the no-pT-cut run to isolate the pT cut's own effect.

Unlike the time-resolution scan, BOTH the z0 and the time cut must be
re-derived at every grid point here: apply_angle_smearing() recomputes
z_axis_intercept_mm, psi_transverse_deg and inv_radius_per_mm (and
hence t_corrected_ns, via add_time_of_flight) directly from the
smeared measured direction -- see bib_common.apply_angle_smearing's
own docstring. So both selection cuts depend on sigma_angle here,
unlike the time-resolution scan where z0 was independent of timing and
could be derived once and held fixed. Still uses the same 98%-per-hit
signal containment convention and an N-1 hierarchy: z0 first (N-1 =
no other cuts, see below), then time (N-1 = the (now point-dependent)
z0 cut).

REVISION (per Luciano, after seeing the first version of this scan):
the per-hit momentum_gev selection cut was REMOVED ENTIRELY in an
earlier revision of this script -- see the git/project history. That
revision showed the efficiency-vs-sigma_angle curve dropping to ~67%
at sigma=5 deg, essentially UNCHANGED from the version with the pT cut
held fixed (99.5%->67% either way) -- proving the pT cut was NOT the
cause of that drop.

REVISION2 (the actual bug, per Luciano: "I do not understand... the
same [flat efficiency] must happen for this scan in angle. So
something does not check"): the z0/time cuts were only being
re-derived at 98% containment for the IT/OT BARREL systems, while
angle smearing (and hence z0/time degradation) is applied tracker-wide
-- so IT/OT ENDCAP hits kept failing their frozen sigma=0 cuts as
sigma_angle grew, and since min_hits_found counts across all four
non-vertex systems (IT/OT barrel+endcap), this alone produced the
entire spurious efficiency drop. Fixed by re-deriving z0/time cuts for
every counted system. With that fix alone (pT cut still removed),
efficiency came back flat (~98-99%) across the whole scan, as it must
by construction of the 98%-containment convention -- exactly like the
time-resolution scan.

With the endcap bug now fixed, the per-hit momentum_gev selection cut
(the real analysis's own single-hit curvature-based pT pre-filter,
cuts_config's third column) is REINTRODUCED in THIS version, held
fixed at today's threshold VALUE (not re-derived -- it's a physical pT
selection, not a containment-derived cut) and used as the N-1
population for the z0 cut, matching the very first version's
convention. Goal: isolate the pT cut's own effect on efficiency and
E[#fakes] now that the endcap bug is no longer confounding the
comparison. The pointing angle is still used exactly as before for the
z0 cut (and, unavoidably, for the time-of-flight correction underlying
the time cut, since both come from the same smeared direction). Track
efficiency continues to use MC-truth pT for the x-axis extrapolation.

A single load+baseline-smear pass is done once for signal and BIB;
apply_angle_smearing() and add_time_of_flight() are then re-run at
each scan point directly on the same hits dict (both are idempotent
functions of the never-mutated truth momentum/geometry fields, not of
whatever was left over from the previous iteration) -- add_incidence_angles()
is also re-run first at each point, purely for robustness/order-
independence (so a sigma=0 point always means "truth baseline", not
"whatever the previous iteration left behind").

NOTE on scope: "pointing angle resolution" is applied uniformly across
ALL systems here (not restricted to IT/OT barrel, unlike the timing
scan) -- angular reconstruction precision is treated as a tracker-wide
quantity rather than per-subsystem hardware. Position and timing
resolution stay at today's per-system baseline throughout. If a
per-subsystem angle-resolution scan is wanted instead, this is the
place to restrict the smearing the same way the timing scan restricted
sigma_t to IT_BARREL/OT_BARREL.

REVISION3 (per Luciano, for consistency with the fake-rate side of this
calculation): the E[#fakes] side of this script has ALWAYS been
defined purely for the 6-layer IT+OT barrel tower (the P_TRUE constants
come from real_tower_calibration.py's own barrel-only MC), but the
efficiency side had always used the FULL tracker (4 systems, 17 layers
total: IT barrel=3, IT endcap=7, OT barrel=3, OT endcap=4) with
min_hits_found=5 -- i.e. the efficiency and the fake rate reported
side by side never actually described the same geometric object. Fixed
here by restricting BOTH the muon population and the "found" condition
to the barrel tower:
  - population: only muons with ZERO truth hits in IT_ENDCAP/OT_ENDCAP
    are kept for the efficiency denominator (a purely geometric
    selection on raw/truth hit system IDs, made before any smearing or
    cuts, so it is not circular with reconstruction efficiency). Only
    33.7% of the full muon-gun sample survives this cut (checked
    directly: 66.3% of events have >=1 endcap hit) -- still ~34k of the
    100k generated events, plenty for the pT->inf extrapolation fit.
  - found condition: now counts hits ONLY in SCAN_SYSTEMS (IT_BARREL,
    OT_BARREL -- 6 layers total), not the other two (endcap) systems.
  - min_hits_found: kept at 5 (now "5 of 6" barrel layers, a much
    tighter requirement than the old "5 of 17" full-tracker one) --
    per Luciano, start here and compare against 4/6 if needed.
z0/time cut re-derivation (REVISION2 above) is correspondingly scoped
back down to SCAN_SYSTEMS only, since endcap systems no longer enter
the found-mask at all.

Writes angle_resolution_scan_results.csv next to this script.
"""
import sys, os, time, csv
import numpy as np

CODE = os.path.expanduser("~/mnt/code/MuonCollider")
ANALYSIS = os.path.expanduser("~/mnt/Analysis")
DATA = os.path.expanduser("~/mnt/Data")
sys.path.insert(0, CODE)
import bib_common as bc

B_FIELD_T = 5.0
IT_BARREL, OT_BARREL = 3, 5
SCAN_SYSTEMS = (IT_BARREL, OT_BARREL)  # the fake-rate tower's systems (density/efakes computed here)
IT_ENDCAP, OT_ENDCAP = 4, 6
CONTAINMENT = 0.98
MIN_HITS_FOUND_TOWER = 5  # per Luciano: start with 5 of 6 barrel layers (try 4/6 next if needed)

CUTS_PATH = f"{ANALYSIS}/__cuts_config.txt"
SMEAR_PATH = f"{ANALYSIS}/__smearing_config.txt"
SIGNAL_FILE = f"{DATA}/ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root"
BIB_FILE = f"{DATA}/ntu_bib_ipp_3evt.root"

# Real per-layer areas (mm^2), from step4_cuts/density_per_layer_before_after_cuts.csv
AREA_MM2 = {
    (IT_BARREL, 0): 1043723.52, (IT_BARREL, 1): 2203416.32, (IT_BARREL, 2): 5167881.04,
    (OT_BARREL, 0): 14003290.56, (OT_BARREL, 1): 19482839.04, (OT_BARREL, 2): 28615419.84,
}
# Synthetic-tower plane areas W_i^2 (mm^2), Section 14's 1 m^2-outermost projective tower
Y_TRUE = {(IT_BARREL, 0): 164.0, (IT_BARREL, 1): 354.0, (IT_BARREL, 2): 554.0,
          (OT_BARREL, 0): 819.0, (OT_BARREL, 1): 1153.0, (OT_BARREL, 2): 1486.0}
W2_SYNTH = {k: (1000.0 * y / 1486.0) ** 2 for k, y in Y_TRUE.items()}

# Already-calibrated P_true(cut) for this real-radii tower (real_tower_calibration.py)
P_TRUE = {
    "line_B0": 3.6895e-28,
    "conservative_quad": 6.2186e-26,
    "exact_helix": 7.3054e-26,
}

t0 = time.time()
print("Loading cuts/smearing config...")
cuts_today = bc.load_cuts(CUTS_PATH)
track_params = bc.load_track_params(CUTS_PATH)
min_hits_found = track_params["min_hits_found"]
exclude_vertex_hits = track_params["exclude_vertex_hits"]
smear_today = bc.load_smearing_config(SMEAR_PATH)

rng_geom = bc.smearing_rng(smear_today)  # seed 42, shared position/time smearing (today's convention)


def load_and_baseline_smear(root_file, muon_hits_only, restrict_systems=None):
    hits = bc.load_hits(root_file, muon_hits_only=muon_hits_only)
    # Restrict to the systems actually needed BEFORE any derived-field
    # computation -- this machine has only ~3.8GB RAM and no swap; see
    # time_resolution_scan.py's header comment for the full story (a
    # silent OOM-kill was traced to exactly this ordering mistake).
    if restrict_systems is not None:
        mask = np.isin(hits["system"], restrict_systems)
        hits = bc.mask_hits(hits, mask)
    # Position + time smearing only, at today's baseline values -- held
    # fixed for the whole scan (angle smearing is applied/re-applied
    # separately, per scan point, below).
    smear_pos_time_only = {
        "position": dict(smear_today["position"]),
        "time": dict(smear_today["time"]),
        "angle_u": {"sigma": 0.0, "enabled": False},
        "angle_v": {"sigma": 0.0, "enabled": False},
        "seed": smear_today["seed"],
    }
    bc.apply_position_time_smearing(hits, smear_pos_time_only, rng_geom)
    return hits


print("Loading signal ntuple (all systems, for track efficiency)...")
sig = load_and_baseline_smear(SIGNAL_FILE, muon_hits_only=True, restrict_systems=None)
n_events = sig["_n_events"]
print(f"  {len(sig['x']):,} signal hits, {n_events:,} events, {time.time()-t0:.0f}s elapsed")

print("Loading BIB ntuple (IT/OT barrel only, for density)...")
bib = load_and_baseline_smear(BIB_FILE, muon_hits_only=False, restrict_systems=SCAN_SYSTEMS)
print(f"  {len(bib['x']):,} BIB hits (IT+OT barrel), {time.time()-t0:.0f}s elapsed")

# pT-only pass mask is recomputed fresh each scan point below (inv_radius_per_mm
# -- and hence momentum_gev -- changes with sigma_angle), so no fixed version here.

# Per REVISION3: efficiency is now scoped to the same 6-layer IT+OT barrel
# tower as the fake rate (SCAN_SYSTEMS), not the full 4-system tracker.
counted_systems = list(SCAN_SYSTEMS)

# Barrel-confined muon population (geometric, truth-hit-based, computed once
# -- system/layer assignment is fixed truth geometry, unaffected by any of
# the smearing applied in the scan loop below): keep only events with ZERO
# truth hits in IT_ENDCAP/OT_ENDCAP.
has_endcap_hit = np.isin(sig["system"], (IT_ENDCAP, OT_ENDCAP))
events_with_endcap = np.unique(sig["event_id"][has_endcap_hit])
barrel_confined = np.ones(n_events, dtype=bool)
barrel_confined[events_with_endcap] = False
print(f"  barrel-confined muons (no endcap hit at all): {barrel_confined.sum():,} / {n_events:,} "
      f"({100.0*barrel_confined.sum()/n_events:.1f}%)")

px_evt, py_evt = sig["_primary"]["px"], sig["_primary"]["py"]
pT_gen = np.sqrt(np.asarray(px_evt, dtype=float) ** 2 + np.asarray(py_evt, dtype=float) ** 2)
pT_gen_barrel = pT_gen[barrel_confined]

sigma_angle_grid = np.arange(0.0, 5.0001, 0.5)  # deg, 11 points
rows = []
for i, sigma_angle in enumerate(sigma_angle_grid):
    rng_scan_sig = np.random.default_rng(3000 + i)
    rng_scan_bib = np.random.default_rng(4000 + i)
    angle_cfg = {
        "angle_u": {"sigma": float(sigma_angle), "enabled": True},
        "angle_v": {"sigma": float(sigma_angle), "enabled": True},
    }

    for hits, rng_scan in ((sig, rng_scan_sig), (bib, rng_scan_bib)):
        cfg = {"angle_u": angle_cfg["angle_u"], "angle_v": angle_cfg["angle_v"]}
        bc.add_incidence_angles(hits)              # reset to truth baseline (order-independence)
        bc.apply_angle_smearing(hits, cfg, rng_scan)  # recompute z0/psi/inv_radius from smeared direction
        bc.add_time_of_flight(hits, B_FIELD_T=B_FIELD_T)  # recompute t_corrected_ns from the new psi/inv_radius

    # pT cut (momentum_gev), REINTRODUCED per Luciano: held fixed at today's
    # threshold VALUE (not re-derived -- it's a physical pT selection, not a
    # containment-derived cut), applied to both signal and BIB, at every
    # sigma_angle. The per-hit momentum_gev value itself still changes with
    # sigma_angle (inv_radius_per_mm depends on the smeared direction), so
    # the mask it produces is point-dependent even though the threshold is
    # not. Used as the N-1 population for the z0 cut below, matching the
    # original (first-version) convention.
    cuts_pt_only = {name: dict(cuts_today[name]) for name in bc.CUT_NAMES}
    for name in ("t_corrected_ns", "z_axis_intercept_mm"):
        cuts_pt_only[name]["enabled"] = False
    pt_only_mask_sig, _ = bc.apply_cuts(sig, cuts_pt_only, B_FIELD_T=B_FIELD_T)
    pt_only_mask_bib, _ = bc.apply_cuts(bib, cuts_pt_only, B_FIELD_T=B_FIELD_T)

    # z0 cut (98% containment, N-1 = pT cut only), re-derived at this
    # sigma_angle.
    #
    # IMPORTANT: re-derived for EVERY counted (non-vertex) system, not just
    # the IT/OT-barrel pair used for the fake-rate density calculation
    # below. Angle smearing here is applied tracker-wide (see header NOTE
    # on scope), so IT/OT *endcap* hits' z0/time distributions degrade with
    # sigma_angle exactly like the barrel ones do. min_hits_found counts
    # across all of counted_systems (IT barrel+endcap, OT barrel+endcap),
    # so leaving the endcap cuts frozen at their sigma=0 threshold -- as an
    # earlier version of this script did -- makes endcap hits fail an
    # increasingly-too-tight cut as sigma_angle grows, producing a spurious
    # efficiency drop that has nothing to do with real track-finding: it
    # would violate the by-construction 98%-containment guarantee exactly
    # the way it should not (per Luciano, after seeing that version's
    # numbers -- the same guarantee that kept the time-resolution scan's
    # efficiency curve flat applies here once every relevant system's cuts
    # are actually re-derived).
    z0_limit = {}
    for sysid in counted_systems:
        sel = (sig["system"] == sysid) & pt_only_mask_sig
        vals = sig["z_axis_intercept_mm"][sel]
        vals = vals[np.isfinite(vals)]
        z0_limit[sysid] = float(np.percentile(np.abs(vals), CONTAINMENT * 100))

    sig_z0_pass = np.zeros(len(sig["system"]), dtype=bool)
    for sysid in counted_systems:
        m = sig["system"] == sysid
        sig_z0_pass |= m & (np.abs(sig["z_axis_intercept_mm"]) <= z0_limit[sysid])
    # systems not in counted_systems (e.g. vertex) are irrelevant to the
    # N-1 population for the time cut below and to the found-mask, so they
    # are simply left False here -- they never enter combined_mask's count.

    # time cut (98% containment, N-1 = pT cut AND the (point-dependent) z0 cut)
    time_limit = {}
    for sysid in counted_systems:
        sel = (sig["system"] == sysid) & pt_only_mask_sig & sig_z0_pass
        vals = sig["t_corrected_ns"][sel]
        vals = vals[np.isfinite(vals)]
        time_limit[sysid] = float(np.percentile(np.abs(vals), CONTAINMENT * 100))

    # --- track efficiency (pT -> inf), full-detector combined cuts ---
    # momentum_gev REINTRODUCED here, at today's fixed threshold values
    # (cuts_today's own per_system/enabled state, unmodified).
    cuts_now = {name: {"per_system": dict(cuts_today[name]["per_system"]),
                        "enabled": cuts_today[name]["enabled"],
                        "zoom_per_system": cuts_today[name]["zoom_per_system"]}
                for name in bc.CUT_NAMES}
    for sysid in counted_systems:
        cuts_now["z_axis_intercept_mm"]["per_system"][sysid] = z0_limit[sysid]
        cuts_now["t_corrected_ns"]["per_system"][sysid] = time_limit[sysid]

    # --- track efficiency, scoped to the 6-layer IT+OT barrel tower (REVISION3) ---
    # found-mask restricted to SCAN_SYSTEMS (IT_BARREL, OT_BARREL) only --
    # endcap hits no longer count toward "found", consistent with the
    # barrel-confined population used for the denominator below.
    combined_mask, _ = bc.apply_cuts(sig, cuts_now, B_FIELD_T=B_FIELD_T)
    count_mask = combined_mask & np.isin(sig["system"], SCAN_SYSTEMS)
    n_hits_survive = np.bincount(sig["event_id"][count_mask], minlength=n_events)
    found_all = n_hits_survive >= MIN_HITS_FOUND_TOWER  # "5 of 6" barrel layers
    # efficiency denominator/numerator both restricted to barrel-confined muons
    found_barrel = found_all[barrel_confined]
    fit_pt_min = bc.pt_inf_fit_min(cuts_now, counted_systems)
    eff_inf, eff_inf_unc, n_fit = bc.efficiency_at_infinite_pt(pT_gen_barrel, found_barrel, fit_pt_min)

    # --- BIB density per layer, at this sigma_angle ---
    bib_z0_pass = (np.abs(bib["z_axis_intercept_mm"]) <= np.where(
        bib["system"] == IT_BARREL, z0_limit[IT_BARREL], z0_limit[OT_BARREL]))
    bib_time_pass = (np.abs(bib["t_corrected_ns"]) <= np.where(
        bib["system"] == IT_BARREL, time_limit[IT_BARREL], time_limit[OT_BARREL]))
    # momentum_gev REINTRODUCED -- pT pre-filter applied to BIB again too
    bib_pass = bib_z0_pass & bib_time_pass & pt_only_mask_bib

    n_i = {}
    for sysid in SCAN_SYSTEMS:
        for layer in (0, 1, 2):
            sel = bib_pass & (bib["system"] == sysid) & (bib["layer"] == layer)
            n_after = int(np.sum(sel))
            density = n_after / AREA_MM2[(sysid, layer)]
            n_i[(sysid, layer)] = density * W2_SYNTH[(sysid, layer)]

    prod_n = 1.0
    for v in n_i.values():
        prod_n *= v

    row = dict(sigma_angle_deg=float(sigma_angle),
               z0_limit_IT=z0_limit[IT_BARREL], z0_limit_OT=z0_limit[OT_BARREL],
               time_limit_IT=time_limit[IT_BARREL], time_limit_OT=time_limit[OT_BARREL],
               prod_n=prod_n, eff_inf=eff_inf, eff_inf_unc=eff_inf_unc, n_fit=n_fit)
    for name, p in P_TRUE.items():
        row[f"efakes_{name}"] = p * prod_n
    rows.append(row)
    print(f"sigma_angle={sigma_angle:4.1f} deg  z0 IT/OT={z0_limit[IT_BARREL]:.2f}/{z0_limit[OT_BARREL]:.2f} mm  "
          f"time_cut IT/OT={time_limit[IT_BARREL]:.4f}/{time_limit[OT_BARREL]:.4f} ns  "
          f"prod_n={prod_n:.3e}  eff_inf={eff_inf*100:.2f}+-{eff_inf_unc*100:.2f}%  "
          f"E[fakes] exact={row['efakes_exact_helix']:.3e}   ({time.time()-t0:.0f}s elapsed)")

out_csv = f"{ANALYSIS}/claude_scratch/angle_resolution_scan_results.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("Wrote", out_csv)
