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


def _latest_tree(f):
    candidates = sorted(
        [k for k in f.keys() if k.split(";")[0] == "HTAtree"],
        key=lambda k: int(k.split(";")[1]),
    )
    return candidates[-1]


def merge(in_paths, out_path):
    """
    Merge the single-entry input files `in_paths` (in that order) into
    `out_path`, one entry per input. The output is written to a temporary
    file next to it and only renamed into place once complete, so an
    interrupted merge never leaves a half-written file behind.
    """
    import os

    def read_all_branches(path):
        f = uproot.open(path)
        tree_name = _latest_tree(f)
        tree = f[tree_name]
        print(f"  {os.path.basename(str(path))}: {tree.num_entries} "
              f"entr{'y' if tree.num_entries == 1 else 'ies'}, "
              f"{sum(int(n) for n in tree['n_hit'].array(library='np')):,} hits")
        return tree.arrays(library="ak")

    print("Reading input files...")
    arrs = [read_all_branches(p) for p in in_paths]
    combined = ak.concatenate(arrs, axis=0)
    del arrs
    tmp_path = str(out_path) + ".partial"
    print(f"Writing {os.path.basename(str(out_path))} ...")
    with uproot.recreate(tmp_path) as fout:
        fout["HTAtree"] = combined
    del combined

    # Sanity check: re-open and verify the entry count and hits
    t_check = uproot.open(tmp_path)["HTAtree"]
    n_hit_check = t_check["n_hit"].array(library="np")
    os.replace(tmp_path, str(out_path))
    print(f"Wrote {os.path.basename(str(out_path))}: {t_check.num_entries} entries, "
          f"{int(n_hit_check.sum()):,} total hits (per entry: "
          f"{', '.join(f'{int(n):,}' for n in n_hit_check)})")


def main():
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <in1.root> <in2.root> [<in3.root> ...] <output.root>")
        sys.exit(1)
    *in_paths, out_path = sys.argv[1:]
    merge(in_paths, out_path)


if __name__ == "__main__":
    main()
