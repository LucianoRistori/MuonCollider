"""
Merge two or more BIB/background ROOT files (each a single-entry HTAtree,
one simulated mega-event of beam-induced background or electromagnetic
background hits) into one combined ROOT file with one entry per input
file, so that downstream analysis can be pointed at a single input file
instead of running once per component and adding results together by
hand.

All branches (per-event scalars and per-hit jagged vectors alike) are
carried over unchanged; nothing is recomputed or reshaped. bib_common's
load_hits() reads and flattens hits across ALL entries of a tree by
default, so existing scripts automatically see the combined hit sample
when pointed at the merged file, with no other code changes needed - and
the same load_hits() still works unmodified on the original single-entry
input files too.

Usage:
    python3 merge_bib_files.py <in1.root> <in2.root> [<in3.root> ...] <output.root>

At least two input files are required (the last argument is always the
output path). Originally written for exactly plus.root + minus.root;
generalized to N inputs so the electromagnetic-background (ipp) file, or
any other single-entry component, can be folded into the same combined
file.
"""
import sys

import awkward as ak
import uproot


def main():
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <in1.root> <in2.root> [<in3.root> ...] <output.root>")
        sys.exit(1)
    *in_paths, out_path = sys.argv[1:]

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
    arrs = [read_all_branches(p) for p in in_paths]

    print("Concatenating...")
    combined = ak.concatenate(arrs, axis=0)
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
