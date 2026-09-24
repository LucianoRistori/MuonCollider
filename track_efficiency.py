"""
Step 4 (part 2): tracking efficiency vs. generated transverse momentum.

For each simulated muon-gun track (one event = one generated muon), this
script applies the selection cuts (cuts_config.txt, via
bib_common.load_cuts/apply_cuts) and then counts how many of that
track's own hits survive the cuts. A track is counted as "found" if that
count is >= min_hits_found (cuts_config.txt, [track] section, default 5).

The available muon-gun sample is not a single fixed pT - it's generated
FLAT IN 1/pT, spanning pT = 1.5 GeV/c up to a long high-pT tail
(confirmed on data: bins evenly spaced in 1/pT come out with equal
track counts per bin). So tracking efficiency vs. pT is built directly
from this one file: each track's own generated pT (from part_px/part_py
at the primary vertex) is computed, tracks are grouped into bins evenly
spaced in 1/pT (matching the generation - this keeps roughly equal
statistics per bin all the way out into the high-pT tail, unlike
equal-width bins in pT itself, which would leave the tail nearly
empty), and efficiency = (# found)/(# examined) is computed per bin.
Multiple input files are still supported (their tracks are simply
pooled before binning), in case additional muon-gun samples covering a
different pT range are added later.

The plot's x-axis is linear in 1/pT (not pT, and not logarithmic) -
the same convention used for the curvature/pT axis in
incidence_angle_plots.py: tick *positions* stay evenly spaced in 1/pT,
and only the printed tick *labels* are the nonlinear reciprocal pT
values (in GeV/c). This keeps the axis on the same natural coordinate
as the binning itself (evenly spaced in 1/pT), while still reading off
directly in physical pT.

Usage:
    python3 track_efficiency.py <signal1.root> [signal2.root ...] [--cuts cuts_config.txt] [--out output_dir] [--bins N]
"""
import csv
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bib_common import (
    load_hits, add_incidence_angles, add_time_of_flight,
    load_cuts, apply_cuts, load_track_params, VERTEX_SYSTEM_IDS,
    prepare_output_dir, _display_path,
    load_smearing_config, smearing_rng, apply_position_time_smearing,
    apply_angle_smearing,
)

N_BINS_DEFAULT = 150  # evenly spaced in 1/pT; x10 finer than the initial 15


def parse_args(argv):
    signal_files = []
    cuts_config = str(Path(__file__).resolve().parent / "cuts_config.txt")
    outdir = "../output_track_efficiency"
    n_bins = N_BINS_DEFAULT
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--cuts":
            cuts_config = argv[i + 1]
            i += 2
        elif a == "--out":
            outdir = argv[i + 1]
            i += 2
        elif a == "--bins":
            n_bins = int(argv[i + 1])
            i += 2
        else:
            signal_files.append(a)
            i += 1
    if not signal_files:
        print("Usage: python3 track_efficiency.py <signal1.root> [signal2.root ...] "
              "[--cuts cuts_config.txt] [--out output_dir] [--bins N]")
        sys.exit(1)
    return signal_files, cuts_config, Path(outdir), n_bins


def _thin_round_pt_candidates(candidates, pt_lo, pt_hi, x_range, min_gap_frac):
    """Shared thinning logic: from a set of candidate "round" pT values
    (ascending), keep a candidate only once its 1/pT position is at
    least `min_gap_frac` of the visible axis range away from the last
    kept one. Walking from the low-pT (large 1/pT) end outward means
    the high-pT candidates that would otherwise pile up on top of each
    other near 1/pT=0 get thinned out, while the low-pT/turn-on region
    keeps dense, evenly-spaced ticks. Returns a list of (x, pt_value)
    pairs, ascending in pT."""
    candidates = sorted(c for c in set(candidates) if pt_lo <= c <= pt_hi)
    min_gap = min_gap_frac * x_range
    kept = []
    last_x = None
    for c in candidates:
        x = 1.0 / c
        if last_x is None or (last_x - x) >= min_gap:
            kept.append((x, c))
            last_x = x
    return kept


def _nice_minor_step(gap):
    """Round a raw (hi - lo) / n_subdivisions step to the nearest of
    1/2/5 x its order of magnitude, so minor ticks land on round
    numbers too (e.g. a raw step of 1 -> 1, of 2 -> 2, of 6 -> 5)."""
    if gap <= 0:
        return gap
    magnitude = 10 ** np.floor(np.log10(gap))
    return min((1 * magnitude, 2 * magnitude, 5 * magnitude), key=lambda s: abs(s - gap))


def pt_ticks_for_inv_pt_axis(inv_pt_min, inv_pt_max,
                              major_min_gap_frac=0.07, n_minor_per_gap=5):
    """Major and minor tick positions (in 1/pT, 1/(GeV/c)) chosen so
    the printed pT labels land on round-ish numbers, rather than being
    evenly spaced in 1/pT (which prints arbitrary values like 6.75,
    3.37, ...). The axis/binning itself is still linear in 1/pT (see
    module docstring); only where the ticks land is chosen this way.

    Major ticks are labeled 1/2/5 x a power of ten (5, 10, 20, 50,
    100, ...), thinned outward from the low-pT/turn-on end so they
    don't pile up on top of each other near 1/pT=0 (the high-pT end,
    see `_thin_round_pt_candidates`).

    Minor ticks fill in each *gap between consecutive kept major
    ticks* with unlabeled marks at a round step (1/2/5 x an order of
    magnitude, picked to divide that gap into roughly `n_minor_per_gap`
    parts - e.g. 5-10 GeV/c gets minor ticks at 6,7,8,9; 10-20 GeV/c
    at 12,14,16,18; 50-100 GeV/c at 60,70,80,90). Building them gap by
    gap (rather than from a single global 1-9 x 10^k candidate set)
    is what makes them land inside every major interval, including
    ones like 10-20 that a plain "tens" sequence would skip, and keeps
    them from spilling out past the last kept major tick."""
    pt_lo = 1.0 / inv_pt_max if inv_pt_max > 0 else 0.0
    pt_hi = 1.0 / inv_pt_min if inv_pt_min > 0 else float("inf")
    x_range = inv_pt_max - inv_pt_min

    major_candidates = {base * 10 ** decade for decade in range(-1, 7) for base in (1, 2, 5)}
    major = _thin_round_pt_candidates(major_candidates, pt_lo, pt_hi, x_range, major_min_gap_frac)
    major_ticks = [x for x, c in major]
    major_labels = [f"{c:g}" for x, c in major]
    major_pts = [c for x, c in major]  # ascending

    minor_ticks = []
    for lo, hi in zip(major_pts[:-1], major_pts[1:]):
        step = _nice_minor_step((hi - lo) / n_minor_per_gap)
        if step <= 0:
            continue
        v = lo + step
        while v < hi - 1e-9 * hi:
            minor_ticks.append(1.0 / v)
            v += step

    return major_ticks, major_labels, minor_ticks


def load_one_file(signal_file, cuts, min_hits_found, exclude_vertex_hits=False,
                   smear_cfg=None, smear_rng=None):
    """Returns (pT_gen [n_events], found [n_events] bool) for one file's tracks.

    If exclude_vertex_hits is True, hits in the vertex detector (VXD
    barrel/endcap, system in VERTEX_SYSTEM_IDS) don't count toward a
    track's surviving-hit total even when they pass the selection cuts -
    i.e. "found" then requires >=min_hits_found surviving hits *outside*
    the vertex detector."""
    print(f"\nLoading signal file: {_display_path(signal_file)}")
    hits = load_hits(signal_file)
    if smear_cfg is not None:
        apply_position_time_smearing(hits, smear_cfg, smear_rng)
    n_events = hits["_n_events"]
    n_hit = len(hits["x"])
    print(f"  {hits['_tree_name']}: {n_events:,} event(s)/track(s), n_hit={n_hit:,}")

    add_incidence_angles(hits)
    if smear_cfg is not None:
        apply_angle_smearing(hits, smear_cfg, smear_rng)
    add_time_of_flight(hits)
    combined_mask, _ = apply_cuts(hits, cuts)

    px, py = hits["_primary"]["px"], hits["_primary"]["py"]
    pT_gen = np.sqrt(px ** 2 + py ** 2)

    count_mask = combined_mask
    if exclude_vertex_hits:
        system = hits["system"]
        is_vertex = np.isin(system, list(VERTEX_SYSTEM_IDS))
        count_mask = count_mask & ~is_vertex

    event_id = hits["event_id"]
    n_hits_survive = np.bincount(event_id[count_mask], minlength=n_events)
    found = n_hits_survive >= min_hits_found
    return pT_gen, found


def main():
    signal_files, cuts_config, outdir, n_bins = parse_args(sys.argv[1:])
    outdir = prepare_output_dir(outdir)

    smearing_config = os.environ.get("SMEARING_CONFIG", "").strip()
    smear_cfg = None
    smear_rng = None
    if smearing_config:
        smear_cfg = load_smearing_config(smearing_config)
        smear_rng = smearing_rng(smear_cfg)
        applied = [f"{name} sigma={smear_cfg[name]['sigma']}"
                   for name in ("position", "time", "angle_long", "angle_trans")
                   if smear_cfg[name]["enabled"] and smear_cfg[name]["sigma"] > 0]
        joined = ", ".join(applied) if applied else "all disabled/zero"
        print(f"Smearing config: {_display_path(smearing_config)} ({joined})")

    print(f"Loading cuts: {_display_path(cuts_config)}")
    cuts = load_cuts(cuts_config)
    for name, spec in cuts.items():
        state = "enabled" if spec["enabled"] else "DISABLED"
        print(f"  {name:20s}  center={spec['center']:+.4g}  halfwidth={spec['halfwidth']:.4g}  [{state}]")
    track_params = load_track_params(cuts_config)
    min_hits_found = track_params["min_hits_found"]
    exclude_vertex_hits = track_params["exclude_vertex_hits"]
    print(f"  min_hits_found       {min_hits_found}")
    print(f"  exclude_vertex_hits  {exclude_vertex_hits}")

    pT_all, found_all = [], []
    for f in signal_files:
        pT_gen, found = load_one_file(f, cuts, min_hits_found, exclude_vertex_hits,
                                       smear_cfg=smear_cfg, smear_rng=smear_rng)
        pT_all.append(pT_gen)
        found_all.append(found)
    pT_gen = np.concatenate(pT_all)
    found = np.concatenate(found_all)
    n_tracks = len(pT_gen)
    print(f"\nTotal tracks pooled from {len(signal_files)} file(s): {n_tracks:,}, "
          f"pT range [{pT_gen.min():.4g}, {pT_gen.max():.4g}] GeV/c")

    # Bin evenly spaced in 1/pT (matches the flat-in-1/pT generation - see
    # module docstring), so bin statistics stay roughly equal all the way
    # into the high-pT tail instead of thinning out there.
    inv_pT = 1.0 / pT_gen
    edges_inv = np.linspace(inv_pT.min(), inv_pT.max(), n_bins + 1)
    bin_idx = np.clip(np.digitize(inv_pT, edges_inv[1:-1]), 0, n_bins - 1)
    bin_width_inv = edges_inv[1] - edges_inv[0]

    rows = []
    for b in range(n_bins):
        m = bin_idx == b
        n_examined = int(m.sum())
        n_found = int(found[m].sum())
        eff = n_found / n_examined if n_examined else float("nan")
        unc = float(np.sqrt(eff * (1.0 - eff) / n_examined)) if n_examined else float("nan")
        inv_pT_lo, inv_pT_hi = edges_inv[b], edges_inv[b + 1]
        inv_pT_center = 0.5 * (inv_pT_lo + inv_pT_hi)
        # 1/pT bin edges are ascending -> pT edges are descending
        pT_lo, pT_hi = 1.0 / inv_pT_hi, 1.0 / inv_pT_lo
        pT_mean = float(pT_gen[m].mean()) if n_examined else float("nan")
        rows.append({
            "bin": b,
            "inv_pT_lo_per_gev": inv_pT_lo, "inv_pT_hi_per_gev": inv_pT_hi,
            "inv_pT_center_per_gev": inv_pT_center,
            "pT_lo_gev": pT_lo, "pT_hi_gev": pT_hi, "pT_mean_gev": pT_mean,
            "n_tracks_examined": n_examined, "n_tracks_found": n_found,
            "efficiency": eff, "efficiency_unc": unc,
        })

    csv_path = outdir / "track_efficiency_vs_pt.csv"
    fields = ["bin", "inv_pT_lo_per_gev", "inv_pT_hi_per_gev", "inv_pT_center_per_gev",
              "pT_lo_gev", "pT_hi_gev", "pT_mean_gev",
              "n_tracks_examined", "n_tracks_found", "efficiency", "efficiency_unc"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    vertex_note = ", excl. vertex hits" if exclude_vertex_hits else ""
    print(f"\n-- Track-finding efficiency vs. generated pT "
          f"({min_hits_found}+ surviving hits{vertex_note}) --")
    for r in rows:
        print(f"  pT=[{r['pT_lo_gev']:9.3g}, {r['pT_hi_gev']:9.3g}] GeV/c  "
              f"n={r['n_tracks_examined']:>6,}  found={r['n_tracks_found']:>6,}  "
              f"eff={r['efficiency']*100:6.2f}% +/- {r['efficiency_unc']*100:4.2f}%")

    # x-axis: linear in 1/pT (matches the binning), with tick labels
    # relabeled as the reciprocal pT in GeV/c - same convention as the
    # curvature/pT axis in incidence_angle_plots.py.
    inv_pT_center = np.array([r["inv_pT_center_per_gev"] for r in rows])
    eff_arr = np.array([r["efficiency"] * 100 for r in rows])
    unc = [r["efficiency_unc"] * 100 for r in rows]
    xerr = 0.5 * bin_width_inv

    # Crop the x-range to where it's informative: increasing 1/pT (i.e.
    # decreasing pT) runs from the high-efficiency plateau down through
    # the turn-on to the long near-zero tail below threshold. That tail
    # carries no information (efficiency is consistent with 0 all the
    # way down to pT=1.5 GeV/c), so the range is cut just past the last
    # bin (in increasing 1/pT / decreasing pT) with efficiency > 1%,
    # with one bin's margin so that point and its error bar stay clear
    # of the axis edge.
    above_1pct = np.where(eff_arr > 1.0)[0]
    if len(above_1pct) > 0:
        x_hi = inv_pT_center[above_1pct[-1]] + 1.5 * bin_width_inv
    else:
        x_hi = edges_inv[-1]
    x_lo = edges_inv[0]

    tick_pos, tick_labels, minor_tick_pos = pt_ticks_for_inv_pt_axis(x_lo, x_hi)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.errorbar(inv_pT_center, eff_arr, yerr=unc, xerr=xerr,
                fmt="o", markersize=4, capsize=2, color="#2ca858", ecolor="#8fcaa4")
    ax.set_xlim(x_lo, x_hi)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels)
    ax.set_xticks(minor_tick_pos, minor=True)  # unlabeled intermediate ticks
    ax.tick_params(axis="x", which="major", length=7)
    ax.tick_params(axis="x", which="minor", length=3.5)
    ax.set_xlabel("generated transverse momentum pT (GeV/c)  "
                  "[axis linear in 1/pT; cropped where efficiency drops below 1%]")
    ax.set_ylabel("tracking efficiency (%)")
    ax.set_ylim(0, 108)
    ax.set_title(f"Track-finding efficiency vs. pT ({n_bins} bins evenly spaced in 1/pT)\n"
                 f"(found = ≥{min_hits_found} surviving hits after cuts, per track{vertex_note})")
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.12)
    plt.tight_layout()
    plt.savefig(outdir / "track_efficiency_vs_pt.png", dpi=140)
    plt.close(fig)

    print(f"\nWrote track_efficiency_vs_pt.csv and track_efficiency_vs_pt.png "
          f"to {_display_path(outdir.resolve())}")


if __name__ == "__main__":
    main()
