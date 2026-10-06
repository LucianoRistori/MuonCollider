"""
resolution_scan.py - scan of the timing or pointing-angle resolution, with
the z0 and time cuts re-derived at each point from signal containment.
Normally run by ./run_scan (see run_scan.sh); the method is described in
templates/__scan_config.txt, and in Sections 15-16 of the paper (this
script replaces resolution_scans/time_resolution_scan.py and
angle_resolution_scan.py, with every setting read from the config files).

Usage:
    python3 resolution_scan.py --check <scan config>
    python3 resolution_scan.py <scan config> <cuts config> <smearing config>
                               <signal ROOT file> <BIB ROOT file> <output folder>

Writes into <output folder>: scan_results.csv, efakes_vs_<p>.png,
efficiency_vs_<p>.png, cuts_vs_<p>.png and scan_<p>.pdf (p = time / angle).
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

# ---- fake-rate framework constants (Sections 14-16 of the paper) --------
# Real per-layer sensitive areas (mm^2) of the IT/OT barrel layers
AREA_MM2 = {
    (IT_BARREL, 0): 1043723.52, (IT_BARREL, 1): 2203416.32, (IT_BARREL, 2): 5167881.04,
    (OT_BARREL, 0): 14003290.56, (OT_BARREL, 1): 19482839.04, (OT_BARREL, 2): 28615419.84,
}
# Plane areas W_i^2 (mm^2) of the projective tower (1 m^2 outermost plane)
Y_TRUE = {(IT_BARREL, 0): 164.0, (IT_BARREL, 1): 354.0, (IT_BARREL, 2): 554.0,
          (OT_BARREL, 0): 819.0, (OT_BARREL, 1): 1153.0, (OT_BARREL, 2): 1486.0}
W2_SYNTH = {k: (1000.0 * y / 1486.0) ** 2 for k, y in Y_TRUE.items()}
# Calibrated fake probability P_true(cut) for this tower
# (fake_rate_framework/real_tower_calibration.py)
P_TRUE = {"exact_helix": 7.3054e-26, "conservative_quad": 6.2186e-26, "line_B0": 3.6895e-28}
MODEL_LABEL = {"exact_helix": "exact helix", "conservative_quad": "conservative quad",
               "line_B0": "line (B=0)"}

PARAM = {
    "time":  dict(unit="ns", label=r"timing resolution $\sigma_t$", col="sigma_t_ns",
                  show_scale=1000.0, show_unit="ps"),
    "angle": dict(unit="deg", label=r"pointing-angle resolution $\sigma_\theta$",
                  col="sigma_angle_deg", show_scale=1.0, show_unit="deg"),
}


# ------------------------------------------------------------------ config
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
        errors.append(f"parameter = '{parameter}': must be time or angle")

    containment = None
    try:
        containment = float(sec.get("containment", ""))
        if not 0 < containment < 100:
            errors.append(f"containment = {containment}: must be between 0 and 100 (a %)")
    except ValueError:
        errors.append(f"containment = '{sec.get('containment', '')}': not a number")

    grid, kind = None, None
    words = sec.get("grid", "").split()
    if not words:
        errors.append("grid is missing")
    else:
        kind = words[0].lower()
        try:
            nums = [float(w) for w in words[1:]]
        except ValueError:
            nums = None
            errors.append(f"grid = '{' '.join(words)}': the values must be numbers")
        if nums is not None:
            if kind in ("log", "linear"):
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
                errors.append(f"grid = '{' '.join(words)}': must start with log, linear or list")
    if grid is not None and np.any(grid < 0):
        errors.append("grid: resolutions can't be negative")
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


def pass_mask(values, system, limits):
    lim = np.full(32, np.inf)
    for s, x in limits.items():
        lim[s] = x
    with np.errstate(invalid="ignore"):
        return np.abs(values) <= lim[system]


def no_smear(cfg, *names):
    """Copy of a smearing config with the named quantities switched off."""
    c = {k: (dict(v) if isinstance(v, dict) else v) for k, v in cfg.items()}
    for n in names:
        c[n] = {"sigma": 0.0, "enabled": False}
    return c


# ------------------------------------------------------------------ the scan
def run(scan_path, cuts_path, smear_path, signal_file, bib_file, outdir):
    t0 = time.time()
    scan = load_scan_config(scan_path)
    p, containment, grid = scan["parameter"], scan["containment"], scan["grid"]
    cuts_cfg = bc.load_cuts(cuts_path)
    min_hits = bc.load_track_params(cuts_path)["min_hits_found"]
    smear = bc.load_smearing_config(smear_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Scan of the {p} resolution: {scan['grid_text']}  ({len(grid)} points)")
    print(f"Cuts: z0 and time re-derived at {containment:g}% signal containment in every "
          f"subsystem; pT cut as in the cuts config")
    print(f"Track found: >= {min_hits} of its IT/OT barrel hits pass the cuts")
    print(f"Fixed smearing: {bc.describe_smearing(no_smear(smear, 'time') if p == 'time' else no_smear(smear, 'angle_u', 'angle_v'))}")

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
        else:
            # position + time smearing fixed; angle smearing applied per point
            bc.apply_position_time_smearing(hits, no_smear(smear, "angle_u", "angle_v"), rng_geom)
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

    z0_fixed = None
    rows = []
    for i, value in enumerate(grid):
        if p == "time":
            sig["t_corrected_ns"] = sig["t_corrected_ns_base"] + \
                np.random.default_rng(1000 + i).normal(0.0, value, len(sig["system"]))
            bib["t_corrected_ns"] = bib["t_corrected_ns_base"] + \
                np.random.default_rng(2000 + i).normal(0.0, value, len(bib["system"]))
        else:
            a = {"angle_u": {"sigma": float(value), "enabled": True},
                 "angle_v": {"sigma": float(value), "enabled": True}}
            for hits, seed in ((sig, 3000 + i), (bib, 4000 + i)):
                bc.add_incidence_angles(hits)
                bc.apply_angle_smearing(hits, a, np.random.default_rng(seed))
                bc.add_time_of_flight(hits, B_FIELD_T=B_FIELD_T)

        sig_pt, _ = bc.apply_cuts(sig, pt_only, B_FIELD_T=B_FIELD_T)
        bib_pt, _ = bc.apply_cuts(bib, pt_only, B_FIELD_T=B_FIELD_T)

        if p == "angle" or z0_fixed is None:
            z0_fixed = limit_per_system(sig["z_axis_intercept_mm"], sig["system"], sig_pt, containment)
        z0 = z0_fixed
        sig_z0 = pass_mask(sig["z_axis_intercept_mm"], sig["system"], z0)
        tl = limit_per_system(sig["t_corrected_ns"], sig["system"], sig_pt & sig_z0, containment)

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

        # BIB density per barrel layer after the cuts -> E[#fakes]
        bib_pass = bib_pt & pass_mask(bib["z_axis_intercept_mm"], bib["system"], z0) \
                          & pass_mask(bib["t_corrected_ns"], bib["system"], tl)
        prod_n = 1.0
        n_layer = {}
        for (s, layer), area in AREA_MM2.items():
            n_after = int(np.sum(bib_pass & (bib["system"] == s) & (bib["layer"] == layer)))
            n_layer[(s, layer)] = n_after
            prod_n *= n_after / area * W2_SYNTH[(s, layer)]

        row = {PARAM[p]["col"]: float(value)}
        for s in SYSTEMS:
            row[f"z0_cut_mm_{SHORT[s]}"] = z0[s]
        for s in SYSTEMS:
            row[f"time_cut_ns_{SHORT[s]}"] = tl[s]
        row.update(eff_inf=eff, eff_inf_unc=eff_unc, n_fit=n_fit, prod_n=prod_n)
        for m, pt in P_TRUE.items():
            row[f"efakes_{m}"] = pt * prod_n
        for (s, layer), n in n_layer.items():
            row[f"bib_hits_after_{SHORT[s]}_L{layer}"] = n
        rows.append(row)
        sh = PARAM[p]
        print(f"{p} res = {value*sh['show_scale']:8.3g} {sh['show_unit']}:  "
              f"z0 IT/OT = {z0[IT_BARREL]:.1f}/{z0[OT_BARREL]:.1f} mm  "
              f"time IT/OT = {tl[IT_BARREL]:.4f}/{tl[OT_BARREL]:.4f} ns  "
              f"eff = {eff*100:.2f} +- {eff_unc*100:.2f} %  "
              f"E[fakes] (exact helix) = {row['efakes_exact_helix']:.3e}  ({time.time()-t0:.0f}s)")

    with open(outdir / "scan_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("Wrote scan_results.csv")
    make_plots_and_pdf(rows, scan, smear, cuts_cfg, min_hits, outdir)
    print(f"Done ({time.time()-t0:.0f}s)")


# ------------------------------------------------------------------ output
def make_plots_and_pdf(rows, scan, smear, cuts_cfg, min_hits, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    p = scan["parameter"]
    sh = PARAM[p]
    x = np.array([r[sh["col"]] for r in rows]) * sh["show_scale"]
    logx = scan["grid_kind"] == "log" or (np.all(x > 0) and x.max() / x.min() > 50)
    xlabel = f"{sh['label']} ({sh['show_unit']}), all subsystems"
    if p == "time":
        today = smear["time"]["sigma"] if smear["time"]["enabled"] else None
    else:
        u, v = smear["angle_u"], smear["angle_v"]
        today = u["sigma"] if (u["enabled"] and v["enabled"] and u["sigma"] == v["sigma"]) else None
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
    styles = {"exact_helix": ("#1b9e77", "-"), "conservative_quad": ("#d95f02", "--"),
              "line_B0": ("#7570b3", ":")}
    for m, (c, ls) in styles.items():
        ax.plot(x, [r[f"efakes_{m}"] for r in rows], color=c, ls=ls, lw=2, marker="o", ms=4,
                label=MODEL_LABEL[m])
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
    eff = np.array([r["eff_inf"] for r in rows]) * 100
    unc = np.array([r["eff_inf_unc"] for r in rows]) * 100
    ax.errorbar(x, eff, yerr=unc, color="#1b9e77", lw=2, marker="o", ms=5, capsize=3)
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"track-finding efficiency, $p_T\to\infty$ (%)")
    ax.set_title(f"Track-finding efficiency vs {p} resolution\n"
                 f"(barrel-confined muons, >= {min_hits} of their IT/OT barrel hits pass)\n{sub}",
                 fontsize=10)
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
    fixed = no_smear(smear, "time") if p == "time" else no_smear(smear, "angle_u", "angle_v")
    lines = [
        f"Resolution scan: {p}      run {run_name}",
        "",
        f"Scanned:        {sh['label'].replace('$', '').replace(chr(92), '')} = {scan['grid_text']} "
        f"({len(rows)} points), applied to all subsystems",
        f"Fixed smearing: {bc.describe_smearing(fixed)}  (from __smearing_config.txt)",
        f"pT cut (GeV/c): {pt_txt}  (from __cuts_config.txt, fixed)",
        f"z0, time cuts:  re-derived at each point, {scan['containment']:g}% signal containment, every subsystem",
        f"Efficiency:     pT -> inf, barrel-confined muons, found = >= {min_hits} IT/OT barrel hits pass",
        f"E[#fakes]:      6-layer IT+OT barrel tower, BIB density after cuts x calibrated P_true",
        "",
        f"{'res (' + sh['show_unit'] + ')':>10} {'z0 IT':>7} {'z0 OT':>7} {'t IT':>8} {'t OT':>8}"
        f" {'eff (%)':>14} {'E[fakes] helix':>15} {'quad':>10} {'line':>10}",
        f"{'':>10} {'(mm)':>7} {'(mm)':>7} {'(ns)':>8} {'(ns)':>8}",
    ]
    for r, xv in zip(rows, x):
        lines.append(
            f"{xv:10.4g} {r['z0_cut_mm_IT_barrel']:7.1f} {r['z0_cut_mm_OT_barrel']:7.1f} "
            f"{r['time_cut_ns_IT_barrel']:8.4f} {r['time_cut_ns_OT_barrel']:8.4f} "
            f"{r['eff_inf']*100:7.2f} +- {r['eff_inf_unc']*100:4.2f} {r['efakes_exact_helix']:15.3e} "
            f"{r['efakes_conservative_quad']:10.2e} {r['efakes_line_B0']:10.2e}")
    lines += ["", "Cuts for every subsystem: cuts_vs_" + p + ".png and scan_results.csv"]
    title = plt.figure(figsize=(11, 8.5))
    title.text(0.05, 0.95, "\n".join(lines), family="monospace", fontsize=8.5, va="top")
    with PdfPages(outdir / f"scan_{p}.pdf") as pdf:
        pdf.savefig(title)
        for f in figs:
            pdf.savefig(f)
    print(f"Wrote efakes_vs_{p}.png, efficiency_vs_{p}.png, cuts_vs_{p}.png, scan_{p}.pdf")


def main(argv):
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
