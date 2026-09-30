"""
Step 1 of rebuilding docs/efficiency_pt_infinity.pdf (see README.txt).

For every signal hit and every signal muon, the generated pT and whether it
passes - a hit: all three cuts of its own subsystem; a muon: "found" - computed
exactly as apply_cuts.py and track_efficiency.py compute them in a run, so the
fits reproduce the run's numbers. Writes fitdata.npz next to this script.

    python3 export_fitdata.py [analysis folder] [--label "run 2026-09-30_145059_pdf"]

The analysis folder (default ~/Dropbox/Documents/MuonColliderSimulation/Analysis)
holds the three settings files; the data are in the folder above it (Data/).
--label names the settings in the note's subtitle (default: the analysis
folder's settings on today's date).
"""
import datetime
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]                      # docs/src/efficiency_pt_infinity -> the repo
sys.path.insert(0, str(REPO))
import bib_common as bc                     # noqa: E402
import input_files                          # noqa: E402


def normals_per_hit(cfg):
    """How many Gaussian numbers the smearing draws per hit (see
    bib_common.apply_position_time_smearing / apply_angle_smearing)."""
    pos, t = cfg["position"], cfg["time"]
    k = 0
    if pos["enabled"] and (pos["sigma_u"] > 0 or pos["sigma_v"] > 0):
        k += (pos["sigma_u"] > 0) + (pos["sigma_v"] > 0)
    if t["enabled"] and t["sigma"] > 0:
        k += 1
    for name in ("angle_u", "angle_v"):
        if cfg[name]["enabled"] and cfg[name]["sigma"] > 0:
            k += 1
    return k


def smeared_signal(signal_file, cfg, rng, muon_hits_only):
    hits = bc.load_hits(str(signal_file), muon_hits_only=muon_hits_only)
    bc.apply_position_time_smearing(hits, cfg, rng)
    bc.add_incidence_angles(hits)
    bc.apply_angle_smearing(hits, cfg, rng)
    bc.add_time_of_flight(hits)
    return hits


def main(argv):
    label = None
    if "--label" in argv:
        i = argv.index("--label")
        label = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    analysis = Path(argv[0] if argv else "~/Dropbox/Documents/MuonColliderSimulation/Analysis").expanduser()
    cuts_path = str(analysis / "__cuts_config.txt")
    cfg = bc.load_smearing_config(str(analysis / "__smearing_config.txt"))
    paths = input_files.input_paths(analysis / "__input_files_config.txt", analysis.parent)
    cuts = bc.load_cuts(cuts_path)
    muon_hits_only = bc.signal_muon_hits_only(cuts_path)
    today = datetime.date.today()
    if label is None:
        label = f"the settings in {analysis.name}/ on {today:%Y-%m-%d}"

    # Hits, as in apply_cuts.py: its random generator first smears the combined
    # BIB file, so skip the numbers that uses before smearing the signal.
    n_bib = sum(n for _, _, n in input_files._entry_signature(paths["bib_combined"]))
    rng = bc.smearing_rng(cfg)
    left = n_bib * normals_per_hit(cfg)
    while left:
        k = min(left, 4_000_000)
        rng.standard_normal(k)
        left -= k
    sig = smeared_signal(paths["signal"], cfg, rng, muon_hits_only)
    passed, _ = bc.apply_cuts(sig, cuts)
    pt_event = np.sqrt(sig["_primary"]["px"] ** 2 + sig["_primary"]["py"] ** 2)
    event = sig["event_id"]

    # Muons, as in track_efficiency.py: its own generator, the signal only.
    trk = smeared_signal(paths["signal"], cfg, bc.smearing_rng(cfg), muon_hits_only)
    trk_passed, _ = bc.apply_cuts(trk, cuts)
    params = bc.load_track_params(cuts_path)
    counted = trk_passed
    if params["exclude_vertex_hits"]:
        counted = counted & ~np.isin(trk["system"], list(bc.VERTEX_SYSTEM_IDS))
    n_surviving = np.bincount(trk["event_id"][counted], minlength=trk["_n_events"])
    found = n_surviving >= params["min_hits_found"]
    trk_pt = np.sqrt(trk["_primary"]["px"] ** 2 + trk["_primary"]["py"] ** 2)

    out = HERE / "fitdata.npz"
    np.savez_compressed(out, hit_pt=pt_event[event], hit_passed=passed,
                        hit_system=sig["system"].astype(np.int8), hit_ev=event.astype(np.int32),
                        trk_pt=trk_pt, trk_found=found, trk_q=trk["_primary"]["q"],
                        label=label, date=f"{today:%Y-%m-%d}")
    print(f"{len(passed):,} signal hits and {len(found):,} muons -> {out}")
    for s in (1, 5):
        sel = sig["system"] == s
        e, u, n = bc.efficiency_at_infinite_pt(pt_event[event][sel], passed[sel], 10.0, groups=event[sel])
        print(f"  check: {bc.SYSTEM_NAMES[s]} hit efficiency for pT -> inf = {e*100:.3f} +- {u*100:.3f} %")
    e, u, n = bc.efficiency_at_infinite_pt(trk_pt, found, 10.0)
    print(f"  check: track-finding efficiency for pT -> inf = {e*100:.3f} +- {u*100:.3f} %")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
