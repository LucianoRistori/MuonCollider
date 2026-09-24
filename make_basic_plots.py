"""
Step 1: read the BIB n-tuple and produce a few simple diagnostic plots +
a region summary table. All spatial plots show hit DENSITY (hits/mm^2),
not raw hit counts - see bib_common.region_table() for how the sensitive
area of each (system, side, layer) region is estimated from the hits
themselves.

Usage:
    python3 make_basic_plots.py <input.root> [output_dir]
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bib_common import (
    load_hits, region_table, add_peak_density, subsystem_density_table,
    write_logfile, prepare_output_dir, _display_path, SYSTEM_NAMES,
)
import geometry as geom_mod

PEAK_PERCENTILE = 99.0
PEAK_KEY = f"peak_density_p{PEAK_PERCENTILE:.0f}_hits_per_mm2"
# Peak-density bin size is chosen per region (not fixed globally): mean
# hit density spans 5+ orders of magnitude across regions, so one fixed
# bin size is either too fine for sparse regions (peak collapses to
# small-integer, noise-dominated bin counts) or too coarse for dense
# ones. See bib_common.add_peak_density() for the adaptive-sizing rule.
PEAK_TARGET_HITS_PER_BIN = 20.0

GEOMETRY_FILES = (
    "MuSIC_v2.xml", "Vertex_o2_v06_01.xml",
    "InnerTracker_o2_v07_01.xml", "OuterTracker_o2_v07_01.xml",
)


def region_label(row):
    side_tag = "" if row["side"] == 0 else (" +z" if row["side"] == 1 else " -z")
    return f"L{row['layer']}{side_tag}"


def main():
    root_file = sys.argv[1] if len(sys.argv) > 1 else \
        "/mnt/user-data/uploads/ntu_bib_plus_1evt.root"
    outdir = Path(sys.argv[2] if len(sys.argv) > 2 else "../output")
    # geometry dir: 3rd CLI arg, else same directory as the input file
    # (where MuSIC_v2.xml and its includes are expected to live)
    geom_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(root_file).resolve().parent
    outdir = prepare_output_dir(outdir)

    hits = load_hits(root_file)
    p = hits["_primary"]
    n_hit = len(hits["x"])
    pdg_unique = np.unique(p["pdg"])
    pdg_str = str(pdg_unique[0]) if pdg_unique.size == 1 else f"{list(pdg_unique)} (mixed)"
    print(f"Loaded {hits['_tree_name']}: {hits['_n_events']} event(s), "
          f"primary pdg={pdg_str}, n_hit={n_hit:,}")

    rows = region_table(hits)

    geometry_used = False
    if all((geom_dir / f).exists() for f in GEOMETRY_FILES):
        try:
            area_lookup = geom_mod.build_area_lookup(geom_dir)
            geom_mod.annotate_rows_with_geometry_area(rows, area_lookup)
            geometry_used = True
            print(f"Using true geometry-based area from {_display_path(geom_dir)} "
                  f"({len(area_lookup)}/{len(rows)} regions matched).")
        except Exception as e:
            print(f"WARNING: geometry parsing failed ({e}); "
                  f"falling back to hit-inferred area for all regions.")
    else:
        missing = [f for f in GEOMETRY_FILES if not (geom_dir / f).exists()]
        print(f"NOTE: geometry file(s) not found in {_display_path(geom_dir)} "
              f"({', '.join(missing)}); using hit-inferred area (less accurate, "
              f"see bib_common.region_table docstring).")

    add_peak_density(hits, rows, bin_size_mm=None, percentile=PEAK_PERCENTILE,
                      target_hits_per_bin=PEAK_TARGET_HITS_PER_BIN)
    sys_rows = subsystem_density_table(rows)
    sys_ids = sorted(SYSTEM_NAMES.keys())

    # ---- 1. hit density per subsystem: mean vs peak (grouped bar chart) --
    sys_by_id = {r["system"]: r for r in sys_rows}
    mean_d = [sys_by_id[s]["density_hits_per_mm2"] for s in sys_ids]
    peak_d = [sys_by_id[s][PEAK_KEY] for s in sys_ids]
    labels = [SYSTEM_NAMES[s] for s in sys_ids]
    xpos = np.arange(len(sys_ids))
    w = 0.38

    fig, ax = plt.subplots(figsize=(9, 5))
    b1 = ax.bar(xpos - w / 2, mean_d, w, label="mean (total hits / total area)",
               color="#3b7dd8")
    b2 = ax.bar(xpos + w / 2, peak_d, w,
               label=f"peak (p{PEAK_PERCENTILE:.0f}, adaptive bin size, "
                     f"~{PEAK_TARGET_HITS_PER_BIN:.0f} hits/bin target)",
               color="#e05252")
    ax.set_ylabel(r"hit density (hits / mm$^2$)")
    ax.set_yscale("log")
    ax.set_title(f"Hit density per subsystem: mean vs. peak ({n_hit:,} hits total)")
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation=30)
    ax.legend()
    for bars in (b1, b2):
        for b in bars:
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                     f"{b.get_height():.3g}", ha="center", va="bottom", fontsize=7)
    plt.tight_layout()
    plt.savefig(outdir / "density_per_subsystem.png", dpi=130)
    plt.close(fig)

    # ---- 2. hit density per region, barrels: mean vs peak (grouped bars) -
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    for ax, s in zip(axes, [1, 3, 5]):  # VXD/IT/OT barrel
        rs = sorted([r for r in rows if r["system"] == s], key=lambda r: r["mean_r"])
        x = [region_label(r) for r in rs]
        xp = np.arange(len(x))
        ax.bar(xp - w / 2, [r["density_hits_per_mm2"] for r in rs], w,
              label="mean", color="#3b7dd8")
        ax.bar(xp + w / 2, [r[PEAK_KEY] for r in rs], w,
              label="peak", color="#e05252")
        ax.set_xticks(xp)
        ax.set_xticklabels(x)
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_ylabel(r"hits / mm$^2$")
        ax.set_yscale("log")
        ax.set_xlabel("layer (ordered by r)")
        ax.legend(fontsize=8)
    plt.suptitle("Hit density per layer, mean vs. peak - barrel subsystems")
    plt.tight_layout()
    plt.savefig(outdir / "density_per_layer_barrels.png", dpi=130)
    plt.close(fig)

    # ---- 3. hit density per region, endcaps: mean vs peak (grouped bars) -
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    for ax, s in zip(axes, [2, 4, 6]):  # VXD/IT/OT endcap
        rs = sorted([r for r in rows if r["system"] == s], key=lambda r: r["mean_z"])
        x = [region_label(r) for r in rs]
        xp = np.arange(len(x))
        ax.bar(xp - w / 2, [r["density_hits_per_mm2"] for r in rs], w,
              label="mean", color="#3b7dd8")
        ax.bar(xp + w / 2, [r[PEAK_KEY] for r in rs], w,
              label="peak", color="#e05252")
        ax.set_xticks(xp)
        ax.set_xticklabels(x, rotation=45)
        ax.set_title(SYSTEM_NAMES[s])
        ax.set_ylabel(r"hits / mm$^2$")
        ax.set_yscale("log")
        ax.set_xlabel("disk (ordered by z)")
        ax.legend(fontsize=8)
    plt.suptitle("Hit density per disk, mean vs. peak - endcap subsystems")
    plt.tight_layout()
    plt.savefig(outdir / "density_per_disk_endcaps.png", dpi=130)
    plt.close(fig)

    # ---- 4. hit time distribution (raw counts - not an area density) -----
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(hits["t"], bins=200, color="#3b7dd8")
    ax.set_xlabel("hit time t (ns)")
    ax.set_ylabel("hits")
    ax.set_yscale("log")
    ax.set_title("Hit time distribution (all subsystems, raw counts)")
    plt.tight_layout()
    plt.savefig(outdir / "hit_time_distribution.png", dpi=130)
    plt.close(fig)

    # ---- 5. r-z map colored by subsystem (spatial reference, not density)
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = {1: "#e6194B", 2: "#f58231", 3: "#3cb44b",
              4: "#4363d8", 5: "#911eb4", 6: "#42d4f4"}
    rng = np.random.default_rng(0)
    system = hits["system"]
    for s in sys_ids:
        mask = system == s
        idx = np.where(mask)[0]
        if len(idx) > 40000:
            idx = rng.choice(idx, 40000, replace=False)
        ax.scatter(hits["z"][idx], hits["r"][idx], s=0.5,
                   color=colors[s], label=SYSTEM_NAMES[s])
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("r (mm)")
    ax.set_title("Hit map (r vs z), colored by subsystem (subsampled, spatial reference only)")
    ax.legend(markerscale=15, loc="upper right")
    plt.tight_layout()
    plt.savefig(outdir / "rz_map_by_subsystem.png", dpi=140)
    plt.close(fig)

    # ---- region summary table (csv), now with area + density -------------
    # union of keys across rows (in case geometry area is missing for some
    # regions, keeping the csv well-formed either way)
    csv_fields = list(dict.fromkeys(k for r in rows for k in r.keys()))
    csv_path = outdir / "region_summary.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=csv_fields, restval="")
        w.writeheader()
        w.writerows(rows)

    sys_csv_path = outdir / "subsystem_summary.csv"
    with open(sys_csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sys_rows[0].keys()))
        w.writeheader()
        w.writerows(sys_rows)

    # ---- run log ------------------------------------------------------
    log_path = outdir / "run_log.txt"
    write_logfile(log_path, root_file, hits, sys_rows,
                  script_name=Path(__file__).name, outdir=str(outdir.resolve()),
                  area_source="geometry" if geometry_used else "hit-inferred")

    print(f"\nWrote plots, CSVs and {log_path.name} to {_display_path(outdir.resolve())}")
    print(f"{len(rows)} (system, side, layer) regions found.")
    print(f"\nSubsystem densities (hits/mm^2) - mean vs peak (p{PEAK_PERCENTILE:.0f}, "
          f"adaptive bin size targeting ~{PEAK_TARGET_HITS_PER_BIN:.0f} hits/bin):")
    for r in sys_rows:
        print(f"  {r['system_name']:12s} mean={r['density_hits_per_mm2']:.4g}  "
              f"peak={r[PEAK_KEY]:.4g}  (bin={r.get('peak_bin_size_mm', 0):.3g}mm)")


if __name__ == "__main__":
    main()
