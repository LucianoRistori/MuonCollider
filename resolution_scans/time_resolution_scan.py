# CORRECTION (10 October 2026): the depth-view fit in this framework had r and z
# swapped (it fitted r = p + q z, with residuals in r instead of z), so the fake
# probabilities computed or used here are about 450 times too low for the helix.
# Kept as the record of the earlier paper draft; the corrected values come from
# fake_rate_framework/geometric_K.py and are used by fake_rate.py.
"""
time_resolution_scan.py

Scans assumed IT/OT-barrel timing resolution sigma_t from 10 ps to 1 ns
(10 log-spaced points), re-deriving the z0 and time selection cuts at
each point via 98%-per-hit signal containment (|value| <= limit, matching
apply_cuts' own convention), holding position/angle resolution and the
pT cut fixed at today's __smearing_config.txt / __cuts_config.txt
values. Reuses bib_common.py's real functions throughout (load_hits,
add_incidence_angles, add_time_of_flight, apply_position_time_smearing,
apply_angle_smearing, apply_cuts, efficiency_at_infinite_pt) rather than
reimplementing any of the geometry/TOF/efficiency logic.

For each sigma_t: computes (a) the resulting BIB hit density in each of
the 6 IT-barrel/OT-barrel layers -> E[#fakes] via the track-fitting-sim
project's already-calibrated P_true(cut) (3 models), and (b) the
track-finding efficiency for pT -> infinity, with the scanned sigma_t
applied only to IT/OT barrel hits (other systems keep today's 0.03 ns).

REVISION (per Luciano, for consistency with the angle-resolution scan's
REVISION3 fix -- see angle_resolution_scan.py's header for the full
story): the efficiency side of this script originally used the FULL
tracker (4 systems, 17 layers: IT barrel=3, IT endcap=7, OT barrel=3,
OT endcap=4) with min_hits_found=5, while the fake-rate side (the
P_TRUE constants from real_tower_calibration.py) has always been
defined purely for the 6-layer IT+OT barrel tower -- i.e. the
efficiency and fake-rate numbers reported side by side never described
the same geometric object. Fixed identically to the angle scan:
  - population: only muons with ZERO truth hits in IT_ENDCAP/OT_ENDCAP
    are kept for the efficiency denominator (geometric, truth-hit-based,
    computed once -- not circular with reconstruction efficiency).
  - found condition: counts hits ONLY in SCAN_SYSTEMS (IT_BARREL,
    OT_BARREL -- 6 layers total), not the two endcap systems.
  - min_hits_found: MIN_HITS_FOUND_TOWER = 5 ("5 of 6" barrel layers),
    matching the angle scan's choice.
Note this scan was never hit by the angle scan's REVISION2 endcap-cut-
freezing bug, because sigma_t here is only ever applied to IT/OT barrel
hits (unlike the angle scan's tracker-wide smearing) -- endcap z0/time
distributions never moved, so their frozen cuts were still valid. This
REVISION is purely about matching the efficiency population/found-scope
to the barrel tower, for consistency across Sections 15 and 16.

Writes time_resolution_scan_results.csv next to this script.
"""
import sys, os, time
import numpy as np

CODE = os.path.expanduser("~/mnt/code/MuonCollider")
ANALYSIS = os.path.expanduser("~/mnt/Analysis")
DATA = os.path.expanduser("~/mnt/Data")
sys.path.insert(0, CODE)
import bib_common as bc

B_FIELD_T = 5.0
IT_BARREL, OT_BARREL = 3, 5
SCAN_SYSTEMS = (IT_BARREL, OT_BARREL)
IT_ENDCAP, OT_ENDCAP = 4, 6
CONTAINMENT = 0.98
BASELINE_SIGMA_T_NS = 0.03  # today's value, used for every system NOT being scanned
MIN_HITS_FOUND_TOWER = 5  # "5 of 6" barrel layers, matching the angle scan

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

# position+angle smearing config (today's values), time OFF (handled by hand below)
smear_geom_only = {
    "position": dict(smear_today["position"]),
    "time": {"sigma": 0.0, "enabled": False},
    "angle_u": dict(smear_today["angle_u"]),
    "angle_v": dict(smear_today["angle_v"]),
    "seed": smear_today["seed"],
}

rng_geom = bc.smearing_rng(smear_today)  # seed 42, shared for position+angle (today's convention)


def prepare(root_file, muon_hits_only, restrict_systems=None):
    hits = bc.load_hits(root_file, muon_hits_only=muon_hits_only)
    # Restrict to the systems we actually need BEFORE computing the extra
    # per-hit derived fields below -- this machine has only ~3.8GB RAM and
    # no swap, and materializing incidence-angle/smearing/TOF arrays for
    # all ~17M BIB hits (vs. ~7.7M after restricting to IT+OT barrel)
    # pushes peak memory past the limit and gets the process killed with
    # no error message. Filtering early keeps peak RSS around 2.8GB.
    if restrict_systems is not None:
        mask = np.isin(hits["system"], restrict_systems)
        hits = bc.mask_hits(hits, mask)
    bc.apply_position_time_smearing(hits, smear_geom_only, rng_geom)
    bc.add_incidence_angles(hits)
    bc.apply_angle_smearing(hits, smear_geom_only, rng_geom)
    bc.add_time_of_flight(hits, B_FIELD_T=B_FIELD_T)  # hit_t is still raw/unsmeared-in-time here
    hits["t_corrected_ns_base"] = hits["t_corrected_ns"].copy()  # pure geometric residual, no time smear yet
    return hits


print("Loading signal ntuple (all systems, for track efficiency)...")
sig = prepare(SIGNAL_FILE, muon_hits_only=True, restrict_systems=None)
n_events = sig["_n_events"]
print(f"  {len(sig['x']):,} signal hits, {n_events:,} events, {time.time()-t0:.0f}s elapsed")

print("Loading BIB ntuple (IT/OT barrel only, for density)...")
bib = prepare(BIB_FILE, muon_hits_only=False, restrict_systems=SCAN_SYSTEMS)
print(f"  {len(bib['x']):,} BIB hits (IT+OT barrel), {time.time()-t0:.0f}s elapsed")

# Fixed baseline time residual (today's 0.03 ns) for every hit NOT in the scanned systems
rng_base_sig = np.random.default_rng(142)
sig_is_scan = np.isin(sig["system"], SCAN_SYSTEMS)
sig_t_base_smear = sig["t_corrected_ns_base"] + rng_base_sig.normal(0.0, BASELINE_SIGMA_T_NS, len(sig["x"]))

# pT-only pass mask (momentum cut fixed at today's values, z0/time cuts disabled)
cuts_pt_only = {name: dict(cuts_today[name]) for name in bc.CUT_NAMES}
for name in ("t_corrected_ns", "z_axis_intercept_mm"):
    cuts_pt_only[name]["enabled"] = False
pt_only_mask_sig, _ = bc.apply_cuts(sig, cuts_pt_only, B_FIELD_T=B_FIELD_T)

# z0 cut (98% containment, N-1 = pT cut only): FIXED across the whole sigma_t scan,
# since z0 does not depend on time resolution at all.
z0_limit = {}
for sysid in SCAN_SYSTEMS:
    sel = (sig["system"] == sysid) & pt_only_mask_sig
    vals = sig["z_axis_intercept_mm"][sel]
    vals = vals[np.isfinite(vals)]
    z0_limit[sysid] = float(np.percentile(np.abs(vals), CONTAINMENT * 100))
print("z0 98%-containment limits (mm):", z0_limit)

sig_z0_pass = (np.abs(sig["z_axis_intercept_mm"]) <= np.where(
    sig["system"] == IT_BARREL, z0_limit[IT_BARREL],
    np.where(sig["system"] == OT_BARREL, z0_limit[OT_BARREL], np.inf)))
bib_z0_pass = (np.abs(bib["z_axis_intercept_mm"]) <= np.where(
    bib["system"] == IT_BARREL, z0_limit[IT_BARREL], z0_limit[OT_BARREL]))

# Per REVISION: efficiency is now scoped to the same 6-layer IT+OT barrel
# tower as the fake rate (SCAN_SYSTEMS), not the full 4-system tracker.
counted_systems = list(SCAN_SYSTEMS)

# Barrel-confined muon population (geometric, truth-hit-based, computed
# once -- system/layer assignment is fixed truth geometry, unaffected by
# any smearing): keep only events with ZERO truth hits in IT_ENDCAP/OT_ENDCAP.
has_endcap_hit = np.isin(sig["system"], (IT_ENDCAP, OT_ENDCAP))
events_with_endcap = np.unique(sig["event_id"][has_endcap_hit])
barrel_confined = np.ones(n_events, dtype=bool)
barrel_confined[events_with_endcap] = False
print(f"  barrel-confined muons (no endcap hit at all): {barrel_confined.sum():,} / {n_events:,} "
      f"({100.0*barrel_confined.sum()/n_events:.1f}%)")

px_evt, py_evt = sig["_primary"]["px"], sig["_primary"]["py"]
pT_gen = np.sqrt(np.asarray(px_evt, dtype=float) ** 2 + np.asarray(py_evt, dtype=float) ** 2)
pT_gen_barrel = pT_gen[barrel_confined]

sigma_t_grid = np.geomspace(0.010, 1.0, 10)  # ns
rows = []
for i, sigma_t in enumerate(sigma_t_grid):
    rng_scan_sig = np.random.default_rng(1000 + i)
    rng_scan_bib = np.random.default_rng(2000 + i)

    sig_t_scan = sig["t_corrected_ns_base"] + rng_scan_sig.normal(0.0, sigma_t, len(sig["x"]))
    sig_t_now = np.where(sig_is_scan, sig_t_scan, sig_t_base_smear)

    bib_t_now = bib["t_corrected_ns_base"] + rng_scan_bib.normal(0.0, sigma_t, len(bib["x"]))

    # time cut (98% containment, N-1 = pT cut + the fixed z0 cut)
    time_limit = {}
    for sysid in SCAN_SYSTEMS:
        sel = (sig["system"] == sysid) & pt_only_mask_sig & sig_z0_pass
        vals = sig_t_now[sel]
        vals = vals[np.isfinite(vals)]
        time_limit[sysid] = float(np.percentile(np.abs(vals), CONTAINMENT * 100))

    # --- track efficiency (pT -> inf), full-detector combined cuts ---
    cuts_now = {name: {"per_system": dict(cuts_today[name]["per_system"]),
                        "enabled": cuts_today[name]["enabled"],
                        "zoom_per_system": cuts_today[name]["zoom_per_system"]}
                for name in bc.CUT_NAMES}
    for sysid in SCAN_SYSTEMS:
        cuts_now["z_axis_intercept_mm"]["per_system"][sysid] = z0_limit[sysid]
        cuts_now["t_corrected_ns"]["per_system"][sysid] = time_limit[sysid]

    # --- track efficiency, scoped to the 6-layer IT+OT barrel tower (REVISION) ---
    sig["t_corrected_ns"] = sig_t_now
    combined_mask, _ = bc.apply_cuts(sig, cuts_now, B_FIELD_T=B_FIELD_T)
    count_mask = combined_mask & np.isin(sig["system"], SCAN_SYSTEMS)
    n_hits_survive = np.bincount(sig["event_id"][count_mask], minlength=n_events)
    found_all = n_hits_survive >= MIN_HITS_FOUND_TOWER  # "5 of 6" barrel layers
    found_barrel = found_all[barrel_confined]
    fit_pt_min = bc.pt_inf_fit_min(cuts_now, counted_systems)
    eff_inf, eff_inf_unc, n_fit = bc.efficiency_at_infinite_pt(pT_gen_barrel, found_barrel, fit_pt_min)

    # --- BIB density per layer, at this sigma_t ---
    bib_time_pass = (np.abs(bib_t_now) <= np.where(
        bib["system"] == IT_BARREL, time_limit[IT_BARREL], time_limit[OT_BARREL]))
    bib_pass = bib_z0_pass & bib_time_pass  # pT cut already applied at load (restrict_systems kept all; apply now)
    cuts_bib_pt_only = cuts_pt_only
    bib_pt_pass, _ = bc.apply_cuts(bib, cuts_bib_pt_only, B_FIELD_T=B_FIELD_T)
    bib_pass &= bib_pt_pass

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

    row = dict(sigma_t_ns=sigma_t,
               z0_limit_IT=z0_limit[IT_BARREL], z0_limit_OT=z0_limit[OT_BARREL],
               time_limit_IT=time_limit[IT_BARREL], time_limit_OT=time_limit[OT_BARREL],
               prod_n=prod_n, eff_inf=eff_inf, eff_inf_unc=eff_inf_unc, n_fit=n_fit)
    for name, p in P_TRUE.items():
        row[f"efakes_{name}"] = p * prod_n
    rows.append(row)
    print(f"sigma_t={sigma_t*1000:7.1f} ps  time_cut IT/OT={time_limit[IT_BARREL]:.4f}/{time_limit[OT_BARREL]:.4f} ns  "
          f"prod_n={prod_n:.3e}  eff_inf={eff_inf*100:.2f}+-{eff_inf_unc*100:.2f}%  "
          f"E[fakes] exact={row['efakes_exact_helix']:.3e}   ({time.time()-t0:.0f}s elapsed)")

import csv
out_csv = f"{ANALYSIS}/claude_scratch/time_resolution_scan_results.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("Wrote", out_csv)
