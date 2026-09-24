"""
Merge the "plus" and "minus" BIB ROOT files (each a single-entry HTAtree,
one simulated mega-event of beam-induced background hits) into one
combined ROOT file with 2 entries - one per original file - so that
downstream analysis can be pointed at a single input file instead of
always having to run twice and add results together by hand.

All branches (per-event scalars and per-hit jagged vectors alike) are
carried over unchanged; nothing is recomputed or reshaped. bib_common's
load_hits() reads and flattens hits across ALL entries of a tree by
default, so existing scripts automatically see the combined plus+minus
hit sample when pointed at the merged file, with no other code changes
needed - and the same load_hits() still works unmodified on the original
single-entry plus/minus files too.

Usage:
    python3 merge_bib_files.py <plus.root> <minus.root> <output.root>
"""
import sys

import awkward as ak
import uproot


def main():
    if len(sys.argv) != 4:
        print(f"Usage: python3 {sys.argv[0]} <plus.root> <minus.root> <output.root>")
        sys.exit(1)
    plus_path, minus_path, out_path = sys.argv[1:4]

    def read_all_branches(path):
        f = uproot.open(path)
        candidates = sorted(
            [k for k in f.keys() if k.split(";")[0] == "HTAtree"],
            key=lambda k: int(k.split(";")[1]),
        )
        tree_name = candidates[-1]
        tree = f[tree_name]
        print(f"  {path}: {tree_name}, {tree.num_entries} entr{'y' if tree.num_entries == 1 else 'ies'}, "
              f"{sum(int(n) for n in tree['n_hit'].array(library='np')):,} hits")
        return tree.arrays(library="ak")

    print("Reading input files...")
    arrs_plus = read_all_branches(plus_path)
    arrs_minus = read_all_branches(minus_path)

    print("Concatenating...")
    combined = ak.concatenate([arrs_plus, arrs_minus], axis=0)
    print(f"  combined: {len(combined)} entries")

    print(f"Writing {out_path} ...")
    with uproot.recreate(out_path) as fout:
        fout["HTAtree"] = combined

    # Sanity check: re-open and verify entry count + total hits match input
    f_check = uproot.open(out_path)
    t_check = f_check["HTAtree"]
    n_hit_check = t_check["n_hit"].array(library="np")
    print(f"Wrote {out_path}: {t_check.num_entries} entries, "
          f"{int(n_hit_check.sum()):,} total hits "
          f"(per-entry: {list(n_hit_check)})")


if __name__ == "__main__":
    main()
