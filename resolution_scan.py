"""
resolution_scan.py - scan of the timing, pointing-angle or hit-position
resolution, with the z0 and time cuts re-derived at each point from signal
containment; and the "headroom" set of three scans (./run_scan --headroom).
Normally run by ./run_scan (see run_scan.sh); the method is described in
templates/__scan_config.txt, and in Sections 15-16 of the paper (this
script replaces resolution_scans/time_resolution_scan.py and
angle_resolution_scan.py, with every setting read from the config files).

Usage:
    python3 resolution_scan.py --check <scan config>
    python3 resolution_scan.py <scan config> <cuts config> <smearing config>
                               <signal ROOT file> <BIB ROOT file> <output folder>
  headroom (three scans, each a multiple of today's resolution, then one PDF):
    python3 resolution_scan.py --check-headroom <scan config>
    python3 resolution_scan.py --headroom-one <time|angle|position> <scan config>
            <cuts config> <smearing config> <signal> <BIB> <run folder>
    python3 resolution_scan.py --headroom-pdf <scan config> <cuts config>
            <smearing config> <run folder>

Writes into <output folder>: scan_results.csv, efakes_vs_<p>.png,
efficiency_vs_<p>.png, cuts_vs_<p>.png and scan_<p>.pdf (p = time / angle /
position); --headroom-one writes them into <run folder>/scan_<p>/, and
--headroom-pdf writes <run folder>/headroom.pdf (+ its plots).

Position resolution enters twice: through the cuts (z0, pT and time are
computed from the smeared hit positions) and through the fake probability
of the fake-rate framework, which scales as sigma_u^4 sigma_v^4 for the
helix (fake_rate.py).
"""
import configparser
import csv
import sys
import time
from pathlib import Path

import numpy as np

import bib_common as bc

B_FIELD_T = 5.0
SYSTEMS = tuple(sorted(bc.SYSTEM_NAMES))          # 1..6
IT_BARREL, IT_ENDCAP, OT_BARREL, OT_ENDCAP = 3, 4, 5, 6
TOWER = (IT_BARREL, OT_BARREL)                    # the 6-layer fake-rate tower
SHORT = {1: "VXD_barrel", 2: "VXD_endcap", 3: "IT_barrel",
         4: "IT_endcap", 5: "OT_barrel", 6: "OT_endcap"}

# fake-rate framework constants and formula: fake_rate.py
from fake_rate import AREA_MM2, P_TRUE, mu_one_fake
import fake_rate

PARAM = {
    "time":  dict(unit="ns", label=r"timing resolution $\sigma_t$", col="sigma_t_ns",
                  show_scale=1000.0, show_unit="ps"),
    "angle": dict(unit="deg", label=r"pointing-angle resolution $\sigma_\theta$",
                  col="sigma_angle_deg", show_scale=1.0, show_unit="deg"),
    "position": dict(unit="mm", label=r"hit position resolution $\sigma_u$",
                     col="sigma_pos_mm", show_scale=1000.0, show_unit="µm"),
}
# default grids of the headroom scans (multiples of today's value), used for
# any resolution missing from the [headroom] section of __scan_config.txt
HEADROOM_DEFAULT = {"time": "factor 1 10 10", "angle": "factor 1 10 10",
                    "position": "factor 1 10 10"}
HEADROOM_CONTAINMENT = 98.0      # used when there is no __scan_config.txt
HEADROOM_ORDER = ("position", "time", "angle")


# ------------------------------------------------------------------ config
def parse_grid(words):
    """words of a grid setting -> (grid array or None, kind, [errors]).
    Kinds: log/linear <first> <last> <n>, list <values...>, and
    factor <first> <last> <n> (linear, in multiples of today's value)."""
    errors, grid, kind = [], None, None
    if not words:
        return None, None, ["grid is missing"]
    kind = words[0].lower()
    try:
        nums = [float(w) for w in words[1:]]
    except ValueError:
        return None, kind, [f"grid = '{' '.join(words)}': the values must be numbers"]
    if kind in ("log", "linear", "factor"):
        if len(nums) != 3 or nums[2] != int(nums[2]) or nums[2] < 2:
            errors.append(f"grid = '{' '.join(words)}': expected "
                          f"'{kind} <first> <last> <number of points (>=2)>'")
        elif kind == "log" and min(nums[0], nums[1]) <= 0:
            errors.append("grid: a log grid needs values > 0 (use linear to include 0)")
        else:
            f = np.geomspace if kind == "log" else np.linspace
            grid = f(nums[0], nums[1], int(nums[2]))
    elif kind == "list":
        if not nums:
            errors.append("grid = list: give at least one value")
        else:
            grid = np.array(nums)
    else:
        errors.append(f"grid = '{' '.join(words)}': must start with log, linear, list or factor")
    if grid is not None and np.any(grid < 0):
        errors.append("grid: resolutions can't be negative")
    return grid, kind, errors


def load_headroom_config(path):
    """The [headroom] section of __scan_config.txt (containment from [scan]).
    Returns {parameter: scan dict} in HEADROOM_ORDER; raises ValueError."""
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if path in (None, "", "-"):
        cp.read_dict({"scan": {"containment": str(HEADROOM_CONTAINMENT)}})
        path = "(defaults)"
    elif not cp.read(path):
        raise ValueError(f"file not found: {path}")
    errors = []
    try:
        containment = float(cp.get("scan", "containment", fallback=""))
        if not 0 < containment < 100:
            errors.append(f"containment = {containment}: must be between 0 and 100 (a %)")
    except ValueError:
        containment = None
        errors.append("[scan] containment: not a number")
    sec = cp["headroom"] if cp.has_section("headroom") else {}
    for k in sorted(set(sec.keys()) - set(PARAM)):
        errors.append(f"[headroom]: unknown setting '{k}' (expected time, angle, position)")
    out = {}
    for p in HEADROOM_ORDER:
        text = sec.get(p, HEADROOM_DEFAULT[p]) if sec else HEADROOM_DEFAULT[p]
        grid, kind, errs = parse_grid(text.split())
        errors += [f"[headroom] {p}: {e}" for e in errs]
        out[p] = dict(parameter=p, containment=containment, grid=grid, grid_kind=kind,
                      grid_text=" ".join(text.split()))
    if errors:
        raise ValueError("\n".join(f"{path}: {e}" for e in errors))
    return out


def load_scan_config(path):
    """Read and check __scan_config.txt. Returns dict(parameter, containment,
    grid (np array), grid_kind, grid_text); raises ValueError listing every
    problem found."""
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    errors = []
    try:
        if not cp.read(path):
            raise ValueError(f"file not found: {path}")
    except configparser.Error as e:
        raise ValueError(f"cannot read {path}: {e}")
    if not cp.has_section("scan"):
        raise ValueError(f"{path}: no [scan] section")
    sec = cp["scan"]
    unknown = set(sec.keys()) - {"parameter", "containment", "grid"}
    for k in sorted(unknown):
        errors.append(f"unknown setting '{k}' (expected parameter, containment, grid)")

    parameter = sec.get("parameter", "").strip().lower()
    if parameter not in PARAM:
        errors.append(f"parameter = '{parameter}': must be time, angle or position")

    containment = None
    try:
        containment = float(sec.get("containment", ""))
        if not 0 < containment < 100:
            errors.append(f"containment = {containment}: must be between 0 and 100 (a %)")
    except ValueError:
        errors.append(f"containment = '{sec.get('containment', '')}': not a number")

    words = sec.get("grid", "").split()
    grid, kind, errs = parse_grid(words)
    errors += errs
    if errors:
        raise ValueError("\n".join(f"{path}: {e}" for e in errors))
    return dict(parameter=parameter, containment=containment, grid=grid,
                grid_kind=kind, grid_text=" ".join(words))


# ------------------------------------------------------------------ helpers
def limit_per_system(values, system, select, containment):
    """{system: containment-% quantile of |values|} over the selected hits."""
    out = {}
    for s in SYSTEMS:
        v = values[select & (system == s)]
        v = v[np.isfinite(v)]
        out[s] = float(np.percentile(np.abs(v), containment)) if len(v) else float("nan")
    return out


def n1_cuts(sig, sig_pt, containment, tol=1e-3, max_iter=20):
    """z0 and time cuts of every subsystem such that each keeps
    `containment` % of the signal hits in its own N-1 distribution:
    z0 on the hits passing pT and time, time on those passing pT and z0
    (as in the N-1 plots of ./run_all). The two depend on each other, so
    they are found by iterating from z0 on the hits passing pT alone,
    until no cut changes by more than `tol` (relative).
    Returns (z0 limits, time limits, number of iterations)."""
    z, t, sy = sig["z_axis_intercept_mm"], sig["t_corrected_ns"], sig["system"]
    z0 = limit_per_system(z, sy, sig_pt, containment)
    tl = limit_per_system(t, sy, sig_pt & pass_mask(z, sy, z0), containment)
    for it in range(1, max_iter + 1):
        z0_new = limit_per_system(z, sy, sig_pt & pass_mask(t, sy, tl), containment)
        tl_new = limit_per_system(t, sy, sig_pt & pass_mask(z, sy, z0_new), containment)
        change = max(max(abs(z0_new[s] / z0[s] - 1), abs(tl_new[s] / tl[s] - 1)) for s in SYSTEMS)
        z0, tl = z0_new, tl_new
        if change <= tol:
            return z0, tl, it
    print(f"  WARNING: z0/time cuts not converged after {max_iter} iterations "
          f"(last change {change:.1e})")
    return z0, tl, max_iter


def pass_mask(values, system, limits):
    lim = np.full(32, np.inf)
    for s, x in limits.items():
        lim[s] = x
    with np.errstate(invalid="ignore"):
        return np.abs(values) <= lim[system]


def whole_tracker_efficiency(hits, passed, pt_gen, n_events, min_hits, excl_vxd, cuts):
    """Track-finding efficiency for pT -> inf over the whole tracker, as
    track_efficiency.py (./run_all) defines it: all muons; found = at least
    min_hits of its hits pass the cuts (VXD hits not counted if excl_vxd).
    Returns (eff, unc, n_fit)."""
    count = passed
    if excl_vxd:
        count = count & ~np.isin(hits["system"], list(bc.VERTEX_SYSTEM_IDS))
    found = np.bincount(hits["event_id"][count], minlength=n_events) >= min_hits
    counted = [s for s in bc.SYSTEM_NAMES if not (excl_vxd and s in bc.VERTEX_SYSTEM_IDS)]
    return bc.efficiency_at_infinite_pt(pt_gen, found, bc.pt_inf_fit_min(cuts, counted))


def no_smear(cfg, *names):
    """Copy of a smearing config with the named quantities switched off."""
    c = {k: (dict(v) if isinstance(v, dict) else v) for k, v in cfg.items()}
    if len(names) == 1 and isinstance(names[0], tuple):
        names = names[0]
    for n in names:
        c[n] = {"sigma": 0.0, "enabled": False}
    return c


def base_values(p, smear):
    """Today's value of the scanned resolution, (u, v) - (sigma, sigma) for
    time - from the smearing config; 0 where that smearing is off."""
    if p == "time":
        t = smear["time"]
        v = t["sigma"] if t["enabled"] else 0.0
        return v, v
    if p == "angle":
        u, v = smear["angle_u"], smear["angle_v"]
        return (u["sigma"] if u["enabled"] else 0.0), (v["sigma"] if v["enabled"] else 0.0)
    pos = smear["position"]
    return ((pos["sigma_u"], pos["sigma_v"]) if pos["enabled"] else (0.0, 0.0))


def fixed_smearing(p, smear):
    """The smearing config without the scanned quantity."""
    return no_smear(smear, {"time": ("time",), "angle": ("angle_u", "angle_v"),
                            "position": ("position",)}[p])


POS_KEYS = ("x", "y", "z", "u", "v")


# ------------------------------------------------------------------ the scan
def run(scan_path, cuts_path, smear_path, signal_file, bib_file, outdir):
    """scan_path: a scan config file, or an already loaded scan dict."""
    t0 = time.time()
    scan = scan_path if isinstance(scan_path, dict) else load_scan_config(scan_path)
    p, containment, grid = scan["parameter"], scan["containment"], scan["grid"]
    cuts_cfg = bc.load_cuts(cuts_path)
    track_params = bc.load_track_params(cuts_path)
    min_hits = track_params["min_hits_found"]
    excl_vxd = track_params["exclude_vertex_hits"]
    smear = bc.load_smearing_config(smear_path)
    outdir = Path(outdir)
    base_u, base_v = base_values(p, smear)
    if scan["grid_kind"] == "factor":
        if not base_u > 0 or (p != "time" and not base_v > 0):
            raise ValueError(f"grid = factor needs today's {p} resolution > 0 in "
                             f"__smearing_config.txt (it is off or zero)")
        factors = np.asarray(scan["grid"], float)
        grid = factors * base_u
        ratio_v = base_v / base_u        # u and v each scaled from their own value
    else:
        factors = None
        ratio_v = 1.0                    # an explicit grid: the same value on u and v
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Scan of the {p} resolution: {scan['grid_text']}  ({len(grid)} points)")
    print(f"Cuts: z0 and time re-derived in every subsystem, each keeping {containment:g}% of the "
          f"signal hits in its N-1 distribution (passing the other two cuts); pT cut as in the "
          f"cuts config")
    print(f"Track found: barrel tower: >= {min_hits} of its IT/OT barrel hits pass the cuts "
          f"(muons with no IT/OT endcap hit)")
    print(f"             whole tracker: >= {min_hits} of its hits pass the cuts"
          f"{', VXD hits not counted' if excl_vxd else ''} (all muons, as in ./run_all)")
    print(f"Fixed smearing: {bc.describe_smearing(fixed_smearing(p, smear))}")

    rng_geom = bc.smearing_rng(smear)

    def prepare(root_file, muon_hits_only, systems):
        hits = bc.load_hits(root_file, muon_hits_only=muon_hits_only)
        if systems is not None:      # restrict early: saves memory
            hits = bc.mask_hits(hits, np.isin(hits["system"], systems))
        if p == "time":
            # position + angle smearing fixed; time smearing added per point
            fixed = no_smear(smear, "time")
            bc.apply_position_time_smearing(hits, fixed, rng_geom)
            bc.add_incidence_angles(hits)
            bc.apply_angle_smearing(hits, fixed, rng_geom)
            bc.add_time_of_flight(hits, B_FIELD_T=B_FIELD_T)
            hits["t_corrected_ns_base"] = hits["t_corrected_ns"].copy()
        elif p == "angle":
            # position + time smearing fixed; angle smearing applied per point
            bc.apply_position_time_smearing(hits, no_smear(smear, "angle_u", "angle_v"), rng_geom)
        else:
            # time smearing fixed; position (then angles, time of flight) per point
            bc.apply_position_time_smearing(hits, no_smear(smear, "position"), rng_geom)
            for k in POS_KEYS:
                hits[k + "_base"] = hits[k].copy()
        return hits

    print("Loading the signal (all subsystems)...")
    sig = prepare(signal_file, True, None)
    n_events = sig["_n_events"]
    print(f"  {len(sig['system']):,} hits, {n_events:,} events  ({time.time()-t0:.0f}s)")
    print("Loading the BIB (IT + OT barrel only)...")
    bib = prepare(bib_file, False, TOWER)
    print(f"  {len(bib['system']):,} hits  ({time.time()-t0:.0f}s)")

    # efficiency population: muons with no hit in the IT/OT endcaps
    ev_endcap = np.unique(sig["event_id"][np.isin(sig["system"], (IT_ENDCAP, OT_ENDCAP))])
    barrel_only = np.ones(n_events, dtype=bool)
    barrel_only[ev_endcap] = False
    pt_gen = np.hypot(np.asarray(sig["_primary"]["px"], float), np.asarray(sig["_primary"]["py"], float))
    print(f"  barrel-confined muons: {barrel_only.sum():,} of {n_events:,}")

    pt_only = {n: dict(cuts_cfg[n]) for n in bc.CUT_NAMES}
    for n in ("t_corrected_ns", "z_axis_intercept_mm"):
        pt_only[n]["enabled"] = False

    rows = []
    for i, value in enumerate(grid):
        if p == "time":
            sig["t_corrected_ns"] = sig["t_corrected_ns_base"] + \
                np.random.default_rng(1000 + i).normal(0.0, value, len(sig["system"]))
            bib["t_corrected_ns"] = bib["t_corrected_ns_base"] + \
                np.random.default_rng(2000 + i).normal(0.0, value, len(bib["system"]))
        elif p == "angle":
            a = {"angle_u": {"sigma": float(value), "enabled": True},
                 "angle_v": {"sigma": float(value) * ratio_v, "enabled": True}}
            for hits, seed in ((sig, 3000 + i), (bib, 4000 + i)):
                bc.add_incidence_angles(hits)
                bc.apply_angle_smearing(hits, a, np.random.default_rng(seed))
                bc.add_time_of_flight(hits, B_FIELD_T=B_FIELD_T)
        else:
            pos = {"position": {"sigma_u": float(value), "sigma_v": float(value) * ratio_v,
                                "enabled": True},
                   "time": {"sigma": 0.0, "enabled": False}}
            # the same angle smearing at every point (fixed seeds): only the position changes
            for hits, seed, aseed in ((sig, 5000 + i, 7000), (bib, 6000 + i, 7001)):
                for k in POS_KEYS:
                    np.copyto(hits[k], hits[k + "_base"])
                bc.apply_position_time_smearing(hits, pos, np.random.default_rng(seed))
                bc.add_incidence_angles(hits)
                bc.apply_angle_smearing(hits, smear, np.random.default_rng(aseed))
                bc.add_time_of_flight(hits, B_FIELD_T=B_FIELD_T)

        sig_pt, _ = bc.apply_cuts(sig, pt_only, B_FIELD_T=B_FIELD_T)
        bib_pt, _ = bc.apply_cuts(bib, pt_only, B_FIELD_T=B_FIELD_T)

        z0, tl, n_iter = n1_cuts(sig, sig_pt, containment)

        cuts_now = {n: {"per_system": dict(cuts_cfg[n]["per_system"]), "enabled": True,
                        "zoom_per_system": cuts_cfg[n]["zoom_per_system"]} for n in bc.CUT_NAMES}
        cuts_now["momentum_gev"]["enabled"] = cuts_cfg["momentum_gev"]["enabled"]
        for s in SYSTEMS:
            cuts_now["z_axis_intercept_mm"]["per_system"][s] = z0[s]
            cuts_now["t_corrected_ns"]["per_system"][s] = tl[s]

        # track-finding efficiency, pT -> infinity, IT/OT barrel tower
        passed, _ = bc.apply_cuts(sig, cuts_now, B_FIELD_T=B_FIELD_T)
        counted = passed & np.isin(sig["system"], TOWER)
        n_found = np.bincount(sig["event_id"][counted], minlength=n_events)
        found = (n_found >= min_hits)[barrel_only]
        eff, eff_unc, n_fit = bc.efficiency_at_infinite_pt(
            pt_gen[barrel_only], found, bc.pt_inf_fit_min(cuts_now, list(TOWER)))
        # whole tracker, defined exactly as track_efficiency.py (./run_all)
        eff_t, eff_t_unc, n_fit_t = whole_tracker_efficiency(
            sig, passed, pt_gen, n_events, min_hits, excl_vxd, cuts_now)

        # BIB density per barrel layer after the cuts -> E[#fakes]
        bib_pass = bib_pt & pass_mask(bib["z_axis_intercept_mm"], bib["system"], z0) \
                          & pass_mask(bib["t_corrected_ns"], bib["system"], tl)
        n_layer = {(s, layer): int(np.sum(bib_pass & (bib["system"] == s) & (bib["layer"] == layer)))
                   for (s, layer) in AREA_MM2}
        if p == "position":
            sig_uv = (float(value), float(value) * ratio_v)
        else:
            sig_uv = fake_rate.position_sigmas(smear)
        prod_n, ef = fake_rate.efakes(n_layer, AREA_MM2, *sig_uv)

        row = {PARAM[p]["col"]: float(value)}
        if factors is not None:
            row["factor"] = float(factors[i])
        if p == "position":
            row["sigma_u_mm"], row["sigma_v_mm"] = sig_uv
        for s in SYSTEMS:
            row[f"z0_cut_mm_{SHORT[s]}"] = z0[s]
        for s in SYSTEMS:
            row[f"time_cut_ns_{SHORT[s]}"] = tl[s]
        row["n1_iterations"] = n_iter
        row.update(eff_inf=eff, eff_inf_unc=eff_unc, n_fit=n_fit,
                   eff_inf_tracker=eff_t, eff_inf_tracker_unc=eff_t_unc, n_fit_tracker=n_fit_t,
                   prod_n=prod_n)
        for m in P_TRUE:
            row[f"efakes_{m}"] = ef[m]
        for (s, layer), n in n_layer.items():
            row[f"bib_hits_after_{SHORT[s]}_L{layer}"] = n
        rows.append(row)
        sh = PARAM[p]
        print(f"{p} res = {value*sh['show_scale']:8.3g} {sh['show_unit']}:  "
              f"z0 IT/OT = {z0[IT_BARREL]:.1f}/{z0[OT_BARREL]:.1f} mm  "
              f"time IT/OT = {tl[IT_BARREL]:.4f}/{tl[OT_BARREL]:.4f} ns  "
              f"eff barrel/tracker = {eff*100:.2f} +- {eff_unc*100:.2f} / "
              f"{eff_t*100:.2f} +- {eff_t_unc*100:.2f} %  "
              f"E[fakes] (exact helix) = {row['efakes_exact_helix']:.3e}  ({time.time()-t0:.0f}s)")

    with open(outdir / "scan_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("Wrote scan_results.csv")
    make_plots_and_pdf(rows, scan, smear, cuts_cfg, min_hits, excl_vxd, outdir)
    print(f"Done ({time.time()-t0:.0f}s)")
    return rows


# ------------------------------------------------------------------ output
def make_plots_and_pdf(rows, scan, smear, cuts_cfg, min_hits, excl_vxd, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    p = scan["parameter"]
    sh = PARAM[p]
    x = np.array([r[sh["col"]] for r in rows]) * sh["show_scale"]
    logx = scan["grid_kind"] == "log" or (np.all(x > 0) and x.max() / x.min() > 50)
    xlabel = f"{sh['label']} ({sh['show_unit']}), all subsystems"
    bu, bv = base_values(p, smear)
    today = bu if (bu > 0 and (scan["grid_kind"] == "factor" or p == "time" or bu == bv)) else None
    sub = f"(z0 and time cuts re-derived at {scan['containment']:g}% signal containment)"

    def mark_today(ax):
        if today is not None and (today > 0 or not logx):
            ax.axvline(today * sh["show_scale"], color="black", lw=1.0, alpha=0.6)
            ax.text(0.02, 0.97, f"vertical line: value in __smearing_config.txt "
                    f"({today*sh['show_scale']:g} {sh['show_unit']})",
                    transform=ax.transAxes, fontsize=8, va="top")

    figs = []
    # 1. E[#fakes]
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.plot(x, [r["efakes_exact_helix"] for r in rows], color="#1b9e77", lw=2, marker="o", ms=4,
            label="exact-helix track model")
    ax.axhline(1.0, color="gray", lw=0.8)
    ax.set_yscale("log")
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$E[\#\mathrm{fakes}]$ (6-layer IT+OT barrel tower)")
    ax.set_title(f"Expected fake tracks vs {p} resolution\n{sub}", fontsize=11)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, which="both", alpha=0.25)
    mark_today(ax)
    fig.tight_layout()
    fig.savefig(outdir / f"efakes_vs_{p}.png", dpi=150)
    figs.append(fig)

    # 2. efficiency
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    vxd_txt = ", VXD hits not counted" if excl_vxd else ""
    for key, color, label in (
            ("eff_inf", "#1b9e77",
             f"IT+OT barrel tower: muons with no IT/OT endcap hit,\n"
             f">= {min_hits} of their 6 IT/OT barrel hits pass"),
            ("eff_inf_tracker", "#d95f02",
             f"whole tracker (as in ./run_all): all muons,\n>= {min_hits} of their hits pass{vxd_txt}")):
        ax.errorbar(x, np.array([r[key] for r in rows]) * 100,
                    yerr=np.array([r[key + "_unc"] for r in rows]) * 100,
                    color=color, lw=2, marker="o", ms=5, capsize=3, label=label)
    ax.legend(fontsize=8, loc="lower left")
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"track-finding efficiency, $p_T\to\infty$ (%)")
    ax.set_title(f"Track-finding efficiency vs {p} resolution\n{sub}", fontsize=11)
    ax.grid(True, which="both", alpha=0.25)
    mark_today(ax)
    fig.tight_layout()
    fig.savefig(outdir / f"efficiency_vs_{p}.png", dpi=150)
    figs.append(fig)

    # 3. the derived cuts
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.8))
    colors = plt.get_cmap("tab10")
    for k, s in enumerate(SYSTEMS):
        a1.plot(x, [r[f"z0_cut_mm_{SHORT[s]}"] for r in rows], marker="o", ms=3, lw=1.5,
                color=colors(k), label=bc.SYSTEM_NAMES[s])
        a2.plot(x, [r[f"time_cut_ns_{SHORT[s]}"] for r in rows], marker="o", ms=3, lw=1.5,
                color=colors(k), label=bc.SYSTEM_NAMES[s])
    for a, yl in ((a1, "z0 cut (mm)"), (a2, "time cut (ns)")):
        if logx:
            a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel(xlabel, fontsize=9)
        a.set_ylabel(yl)
        a.grid(True, which="both", alpha=0.25)
    a1.legend(fontsize=8)
    fig.suptitle(f"Cuts keeping {scan['containment']:g}% of the signal hits, per subsystem", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / f"cuts_vs_{p}.png", dpi=150)
    figs.append(fig)

    # title page + table, then the three plots
    run_name = outdir.name
    pt = cuts_cfg["momentum_gev"]["per_system"]
    pt_txt = ", ".join(f"{bc.SYSTEM_NAMES[s]} {pt[s]:g}" for s in SYSTEMS if pt.get(s)) or "off"
    fixed = fixed_smearing(p, smear)
    has_f = "factor" in rows[0]
    grid_txt = scan["grid_text"]
    if scan["grid_kind"] == "factor":
        grid_txt += f"  (x today's value: {bu * sh['show_scale']:g}"
        if p != "time" and bv != bu:
            grid_txt += f" on u, {bv * sh['show_scale']:g} on v, each scaled"
        grid_txt += f" {sh['show_unit']})"
    lines = [
        f"Resolution scan: {p}      run {run_name}",
        "",
        f"Scanned:        {sh['label'].replace('$', '').replace(chr(92), '')} = {grid_txt} "
        f"({len(rows)} points), applied to all subsystems",
        f"Fixed smearing: {bc.describe_smearing(fixed)}  (from __smearing_config.txt)",
        f"pT cut (GeV/c): {pt_txt}  (from __cuts_config.txt, fixed)",
        f"z0, time cuts:  re-derived at each point, every subsystem: each keeps {scan['containment']:g}% of the "
        f"signal hits in its N-1 distribution",
        f"Efficiency:     pT -> inf;  barrel: muons with no IT/OT endcap hit, found = >= {min_hits} of "
        f"their IT/OT barrel hits pass",
        f"                tracker (as ./run_all): all muons, found = >= {min_hits} of their hits pass"
        f"{', VXD not counted' if excl_vxd else ''}",
        f"E[#fakes]:      6-layer IT+OT barrel tower, BIB density after cuts x P_true (exact-helix track model,\n"
        f"                calibrated at 0.1 mm and scaled to the position resolution: sigma_u^4 sigma_v^4);\n"
        f"                1 fake at: common factor on all six densities that gives one fake",
        "",
        f"{'factor':>7} " * has_f + f"{'res (' + sh['show_unit'] + ')':>10} {'z0 IT':>7} {'z0 OT':>7} {'t IT':>8} {'t OT':>8}"
        f" {'eff barrel (%)':>15} {'eff tracker (%)':>16} {'E[fakes]':>12} {'1 fake at':>10}",
        f"{'':>8}" * has_f + f"{'':>10} {'(mm)':>7} {'(mm)':>7} {'(ns)':>8} {'(ns)':>8}",
    ]
    for r, xv in zip(rows, x):
        lines.append(
            (f"{r['factor']:6.3g}x " if has_f else "") + f"{xv:10.4g} {r['z0_cut_mm_IT_barrel']:7.1f} {r['z0_cut_mm_OT_barrel']:7.1f} "
            f"{r['time_cut_ns_IT_barrel']:8.4f} {r['time_cut_ns_OT_barrel']:8.4f} "
            f"{r['eff_inf']*100:8.2f} +- {r['eff_inf_unc']*100:4.2f} "
            f"{r['eff_inf_tracker']*100:9.2f} +- {r['eff_inf_tracker_unc']*100:4.2f} "
            f"{r['efakes_exact_helix']:12.3e} {mu_one_fake(r['efakes_exact_helix']):9.0f}x")
    lines += ["", "Cuts for every subsystem: cuts_vs_" + p + ".png and scan_results.csv"]
    title = plt.figure(figsize=(11, 8.5))
    title.text(0.05, 0.95, "\n".join(lines), family="monospace", fontsize=8.5, va="top")
    with PdfPages(outdir / f"scan_{p}.pdf") as pdf:
        pdf.savefig(title)
        for f in figs:
            pdf.savefig(f)
    print(f"Wrote efakes_vs_{p}.png, efficiency_vs_{p}.png, cuts_vs_{p}.png, scan_{p}.pdf")


# ------------------------------------------------------------------ headroom
def crossing_factor(factors, efakes, level=1.0):
    """First factor where E[#fakes] reaches `level`, interpolated linearly in
    log E vs log factor; None if it stays below over the scanned range."""
    f = np.asarray(factors, float)
    e = np.log10(np.maximum(np.asarray(efakes, float), 1e-300))
    lv = np.log10(level)
    if e[0] >= lv:
        return float(f[0])
    for k in range(1, len(f)):
        if e[k] >= lv:
            w = (lv - e[k - 1]) / (e[k] - e[k - 1])
            return float(10 ** (np.log10(f[k - 1]) + w * (np.log10(f[k]) - np.log10(f[k - 1]))))
    return None


def read_rows(path):
    with open(path) as fh:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]


def make_headroom_pdf(scan_cfg_path, cuts_path, smear_path, run_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    run_dir = Path(run_dir).resolve()
    scans = load_headroom_config(scan_cfg_path)
    smear = bc.load_smearing_config(smear_path)
    cuts_cfg = bc.load_cuts(cuts_path)
    tp = bc.load_track_params(cuts_path)
    data = {}
    for p in HEADROOM_ORDER:
        f = run_dir / f"scan_{p}" / "scan_results.csv"
        if not f.exists():
            raise FileNotFoundError(f"missing {f} - did the {p} scan fail?")
        data[p] = read_rows(f)
    color = {"position": "#7570b3", "time": "#1b9e77", "angle": "#d95f02"}
    name = {"position": "position", "time": "time", "angle": "angle"}

    def today_txt(p):
        u, v = base_values(p, smear)
        sh = PARAM[p]
        if p == "time" or u == v:
            return f"{u * sh['show_scale']:g} {sh['show_unit']}"
        return f"{u * sh['show_scale']:g}/{v * sh['show_scale']:g} {sh['show_unit']} (u/v)"

    figs = []
    summary = {}
    for p in HEADROOM_ORDER:
        rows = data[p]
        fac = [r["factor"] for r in rows]
        ef = [r["efakes_exact_helix"] for r in rows]
        c1 = crossing_factor(fac, ef)
        u0, _ = base_values(p, smear)
        summary[p] = dict(c1=c1, v1=(c1 * u0 if c1 else None), e0=ef[0], fmax=fac[-1],
                          vmax=fac[-1] * u0, emax=ef[-1],
                          eff0=rows[0]["eff_inf_tracker"], effmax=rows[-1]["eff_inf_tracker"],
                          effb0=rows[0]["eff_inf"], effbmax=rows[-1]["eff_inf"])

    def three_panes(ykey, ylabel, title, out_png, log=True, eff=False):
        fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
        for ax, p in zip(axes, HEADROOM_ORDER):
            rows, sh, sm = data[p], PARAM[p], summary[p]
            x = np.array([r[sh["col"]] for r in rows]) * sh["show_scale"]
            u0 = base_values(p, smear)[0] * sh["show_scale"]
            if eff:
                ax.errorbar(x, [r["eff_inf_tracker"] * 100 for r in rows],
                            yerr=[r["eff_inf_tracker_unc"] * 100 for r in rows], color=color[p],
                            lw=2, marker="o", ms=4, capsize=2, label="whole tracker")
                ax.plot(x, [r["eff_inf"] * 100 for r in rows], color=color[p], lw=1.2, ls="--",
                        marker=".", label="IT+OT barrel tower")
                ax.legend(fontsize=8, loc="lower left")
                ax.set_title(f"{name[p]}:  {sm['eff0']*100:.1f}% -> {sm['effmax']*100:.1f}% "
                             f"at {sm['fmax']:g}x (whole tracker)", fontsize=10)
            else:
                ax.plot(x, [r[ykey] for r in rows], color=color[p], lw=2, marker="o", ms=5)
                ax.axhline(1.0, color="gray", lw=1)
                if sm["c1"]:
                    ax.plot([sm["v1"] * sh["show_scale"]], [1.0], "o", color=color[p], ms=10,
                            mfc="white", mew=2, zorder=5)
                    t = (f"{name[p]}: 1 fake at {sm['c1']:.1f}x "
                         f"({sm['v1'] * sh['show_scale']:.3g} {sh['show_unit']})")
                else:
                    t = f"{name[p]}: {sm['emax']:.1e} fakes at {sm['fmax']:g}x"
                ax.set_title(t, fontsize=10.5)
            ax.axvline(u0, color="black", lw=1, alpha=0.6)
            ax.text(u0, 0.97, "  today", transform=ax.get_xaxis_transform(), fontsize=8,
                    va="top")
            if not eff:
                ax.text(0.98, 1.0, "1 fake ", transform=ax.get_yaxis_transform(), fontsize=8,
                        color="gray", ha="right", va="bottom")
            ax.set_xlabel(f"{sh['label']} ({sh['show_unit']})")
            sec = ax.secondary_xaxis("top", functions=(lambda v, u=u0: v / u,
                                                       lambda f, u=u0: f * u))
            sec.set_xlabel("x today's value", fontsize=8)
            sec.tick_params(labelsize=8)
            if log:
                ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.25)
        axes[0].set_ylabel(ylabel)
        if log:
            lo = min(min(r[ykey] for r in data[p]) for p in HEADROOM_ORDER)
            hi = max(max(r[ykey] for r in data[p]) for p in HEADROOM_ORDER)
            axes[0].set_ylim(lo / 5, max(hi * 5, 10))
        fig.suptitle(title, fontsize=11.5)
        fig.tight_layout()
        fig.savefig(run_dir / out_png, dpi=150)
        figs.append(fig)

    contain = scans["time"]["containment"]
    three_panes("efakes_exact_helix", r"$E[\#\mathrm{fakes}]$ (6-layer IT+OT barrel tower)",
                "Fake-track headroom in the detector resolutions: one resolution scanned at a time, "
                f"the other two at today's values\n(z0 and time cuts re-derived at {contain:g}% signal "
                "containment at each point; exact-helix model)",
                "headroom_efakes_vs_resolution.png")
    three_panes("eff_inf_tracker", r"track-finding efficiency, $p_T\to\infty$ (%)",
                "Track-finding efficiency along the same scans", "headroom_efficiency_vs_resolution.png",
                log=False, eff=True)

    with open(run_dir / "headroom_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["parameter", "unit", "today", "max_factor", "max_value", "efakes_today",
                    "efakes_at_max", "factor_one_fake", "value_one_fake", "eff_tracker_today",
                    "eff_tracker_at_max", "containment"])
        for p in HEADROOM_ORDER:
            sm, sh = summary[p], PARAM[p]
            w.writerow([p, sh["show_unit"], base_values(p, smear)[0] * sh["show_scale"], sm["fmax"],
                        sm["vmax"] * sh["show_scale"], sm["e0"], sm["emax"], sm["c1"] or "",
                        (sm["v1"] * sh["show_scale"]) if sm["v1"] else "", sm["eff0"],
                        sm["effmax"], contain])

    # title page
    pt = cuts_cfg["momentum_gev"]["per_system"]
    pt_txt = ", ".join(f"{bc.SYSTEM_NAMES[s]} {pt[s]:g}" for s in SYSTEMS if pt.get(s)) or "off"
    run_name = run_dir.parent.name if run_dir.name.startswith("step5") else run_dir.name
    lines = [f"Resolution headroom      run {run_name}", "",
             f"Today's resolutions (__smearing_config.txt): {bc.describe_smearing(smear)}",
             "Each scan multiplies ONE resolution by the factors below, for all subsystems; the other",
             "two stay at today's values. At each point the z0 and time cuts of every subsystem are",
             f"re-derived to keep {scans['time']['containment']:g}% of the signal hits in their N-1 distributions; "
             f"pT cut fixed ({pt_txt} GeV/c).",
             "E[#fakes]: 6-layer IT+OT barrel tower, exact-helix model. Position resolution enters both",
             "through the cuts and through the fake probability (~ sigma_u^4 sigma_v^4); time and angle",
             "only through the cuts.", "",
             f"{'resolution':<11}{'today':>22}{'factors':>18}{'E[fakes] today':>16}"
             f"{'1 fake at':>12}{'E at max':>11}{'eff. tracker today -> max':>28}"]
    for p in HEADROOM_ORDER:
        sm = summary[p]
        c1 = f"{sm['c1']:.1f}x" if sm["c1"] else f"> {sm['fmax']:g}x"
        lines.append(f"{name[p]:<11}{today_txt(p):>22}{scans[p]['grid_text'].replace('factor ', ''):>18}"
                     f"{sm['e0']:>16.2e}{c1:>12}{sm['emax']:>11.1e}"
                     f"{sm['eff0']*100:>18.2f}% -> {sm['effmax']*100:.2f}%")
    lines += ["", "1 fake at: factor on that resolution alone at which E[#fakes] reaches 1 "
              "(interpolated in log-log).",
              "Pages: E[#fakes] and efficiency vs each resolution, then for each resolution its own",
              "scan (E[#fakes], efficiencies and cuts vs the resolution itself); every number is in",
              "scan_<resolution>/scan_results.csv."]
    title = plt.figure(figsize=(11, 8.5))
    title.text(0.04, 0.95, "\n".join(lines), family="monospace", fontsize=8.3, va="top")

    with PdfPages(run_dir / "headroom.pdf") as pdf:
        pdf.savefig(title)
        for f in figs:
            pdf.savefig(f)
        for p in HEADROOM_ORDER:
            for png in (f"efakes_vs_{p}.png", f"efficiency_vs_{p}.png", f"cuts_vs_{p}.png"):
                img = plt.imread(run_dir / f"scan_{p}" / png)
                h, w = img.shape[:2]
                fg = plt.figure(figsize=(11, 8.5))
                ax = fg.add_axes([0.03, 0.03, 0.94, 0.94])
                ax.imshow(img)
                ax.axis("off")
                pdf.savefig(fg)
                plt.close(fg)
    plt.close("all")
    print("Headroom in the detector resolutions (E[#fakes], exact helix; one resolution at a time):")
    for p in HEADROOM_ORDER:
        sm = summary[p]
        print(f"  {name[p]:<9} today {today_txt(p):<18} E = {sm['e0']:.2e};  "
              + (f"1 fake at {sm['c1']:.1f}x" if sm["c1"] else
                 f"below 1 up to {sm['fmax']:g}x (E = {sm['emax']:.1e})"))
    print("Wrote headroom.pdf, headroom_efakes_vs_resolution.png, "
          "headroom_efficiency_vs_resolution.png, headroom_summary.csv")


def main(argv):
    if len(argv) == 2 and argv[0] == "--check-headroom":
        try:
            load_headroom_config(argv[1])
        except ValueError as e:
            print(e, file=sys.stderr)
            return 1
        print("headroom")
        return 0
    if len(argv) == 8 and argv[0] == "--headroom-one":
        p = argv[1]
        scan = load_headroom_config(argv[2])[p]
        run(scan, argv[3], argv[4], argv[5], argv[6], Path(argv[7]) / f"scan_{p}")
        return 0
    if len(argv) == 5 and argv[0] == "--headroom-pdf":
        make_headroom_pdf(*argv[1:])
        return 0
    if len(argv) == 2 and argv[0] == "--check":
        try:
            s = load_scan_config(argv[1])
        except ValueError as e:
            print(e, file=sys.stderr)
            return 1
        print(s["parameter"])
        return 0
    if len(argv) != 6:
        print(__doc__, file=sys.stderr)
        return 2
    run(*argv)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
