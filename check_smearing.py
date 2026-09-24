"""
Investigate whether hit_u/hit_v (and hit_x/y/z) and hit_t already carry
measurement smearing in the current input files, or are exact/unsmeared.

Checks:
1. hit_du / hit_dv: are they populated, all zero, or a small set of fixed
   per-subsystem "design resolution" values?
2. Self-consistency: within a single sensor module (same system/side/
   layer/module/sensor cellID), for pairs of hits, does
       (global_pos_2 - global_pos_1)
   agree with
       (u2-u1)*u_axis + (v2-v1)*v_axis
   using that module's own local axes? If hit_u/hit_v and hit_x/y/z were
   derived consistently from the same (possibly smeared) position, this
   residual should be ~0 (floating-point only). If one representation
   was smeared independently of the other, the residual should be on
   the order of the smearing sigma.
3. hit_t: look at its value distribution/precision for a coarse sense of
   whether it looks like a raw continuous simulation output or has some
   telltale discretization/rounding.
"""
import sys

import numpy as np
import uproot


def decode_id0(id0):
    system = id0 & 0x1F
    side = (id0 >> 5) & 0x3
    layer = (id0 >> 7) & 0x3F
    module = (id0 >> 13) & 0x7FF
    sensor = (id0 >> 24) & 0xFF
    return system, side, layer, module, sensor


def main():
    path = sys.argv[1]
    f = uproot.open(path)
    candidates = sorted(
        [k for k in f.keys() if k.split(";")[0] == "HTAtree"],
        key=lambda k: int(k.split(";")[1]),
    )
    tree = f[candidates[-1]]
    print(f"Reading {path} ...")
    arrs = tree.arrays(
        ["hit_id0", "hit_x", "hit_y", "hit_z", "hit_u", "hit_v",
         "hit_du", "hit_dv", "hit_ux", "hit_uy", "hit_uz",
         "hit_vx", "hit_vy", "hit_vz", "hit_t"],
        library="np",
    )

    def flat(key):
        a = arrs[key]
        if a.dtype == object:
            return np.concatenate([np.asarray(x) for x in a])
        return a

    id0 = flat("hit_id0").astype(np.int64)
    x, y, z = flat("hit_x"), flat("hit_y"), flat("hit_z")
    u, v = flat("hit_u"), flat("hit_v")
    du, dv = flat("hit_du"), flat("hit_dv")
    ux, uy, uz = flat("hit_ux"), flat("hit_uy"), flat("hit_uz")
    vx, vy, vz = flat("hit_vx"), flat("hit_vy"), flat("hit_vz")
    t = flat("hit_t")

    n = len(id0)
    print(f"Total hits: {n:,}")

    print("\n=== hit_du / hit_dv ===")
    print(f"  hit_du: min={du.min():.6g} max={du.max():.6g} mean={du.mean():.6g} "
          f"n_unique(rounded 1e-6)={len(np.unique(np.round(du, 6))):,} "
          f"n_exact_zero={(du == 0).sum():,}")
    print(f"  hit_dv: min={dv.min():.6g} max={dv.max():.6g} mean={dv.mean():.6g} "
          f"n_unique(rounded 1e-6)={len(np.unique(np.round(dv, 6))):,} "
          f"n_exact_zero={(dv == 0).sum():,}")

    system, side, layer, module, sensor = decode_id0(id0)
    print("\n  hit_du by subsystem (system,side->name approx), first 12 distinct (system,layer) groups:")
    subsys_key = system * 10 + (side > 0).astype(int)
    uniq_keys, first_idx = np.unique(subsys_key, return_index=True)
    for k in uniq_keys[:12]:
        mask = subsys_key == k
        print(f"    system={k//10} side_endcap={k%10}: "
              f"hit_du mean={du[mask].mean():.5g} std={du[mask].std():.5g}  "
              f"hit_dv mean={dv[mask].mean():.5g} std={dv[mask].std():.5g}  n={mask.sum():,}")

    print("\n=== Self-consistency: global vs. local(u,v) within a module ===")
    # group by full cellID (system,side,layer,module,sensor) - hits sharing
    # this id are on the exact same rigid sensor element, so u_axis/v_axis
    # should be constant within the group (check that too).
    order = np.argsort(id0, kind="stable")
    id0_sorted = id0[order]
    boundaries = np.flatnonzero(np.diff(id0_sorted)) + 1
    starts = np.concatenate(([0], boundaries))
    ends = np.concatenate((boundaries, [len(id0_sorted)]))

    rng = np.random.default_rng(0)
    group_sizes = ends - starts
    multi_hit_groups = np.flatnonzero(group_sizes >= 2)
    print(f"  {len(multi_hit_groups):,} sensor elements have >=2 hits (out of {len(starts):,} total elements)")

    sample_groups = rng.choice(multi_hit_groups, size=min(20000, len(multi_hit_groups)), replace=False)
    residuals = []
    axis_consistency = []
    for gi in sample_groups:
        s, e = starts[gi], ends[gi]
        idx = order[s:e]
        if len(idx) < 2:
            continue
        i0 = idx[0]
        for i1 in idx[1:]:
            # check local axes are the same for both hits (rigid module assumption)
            axis_diff = abs(ux[i0] - ux[i1]) + abs(uy[i0] - uy[i1]) + abs(uz[i0] - uz[i1])
            axis_consistency.append(axis_diff)

            dglobal = np.array([x[i1] - x[i0], y[i1] - y[i0], z[i1] - z[i0]])
            u_axis = np.array([ux[i0], uy[i0], uz[i0]])
            v_axis = np.array([vx[i0], vy[i0], vz[i0]])
            duv_projected = (u[i1] - u[i0]) * u_axis + (v[i1] - v[i0]) * v_axis
            resid = dglobal - duv_projected
            residuals.append(np.linalg.norm(resid))

    residuals = np.array(residuals)
    axis_consistency = np.array(axis_consistency)
    print(f"  local-axis consistency within a module: max diff = {axis_consistency.max():.3g} "
          f"(should be ~0 if module axes are rigid per sensor element)")
    print(f"  residual |dglobal - d(u,v)_projected| over {len(residuals):,} hit pairs:")
    print(f"    mean={residuals.mean():.6g} mm  median={np.median(residuals):.6g} mm  "
          f"p95={np.percentile(residuals, 95):.6g} mm  max={residuals.max():.6g} mm")
    print(f"    (for reference, mean hit_du={du.mean():.4g} mm, mean hit_dv={dv.mean():.4g} mm)")

    print("\n=== hit_t ===")
    print(f"  n={len(t):,}  min={t.min():.6g} max={t.max():.6g} mean={t.mean():.6g} std={t.std():.6g}")
    # look at the fractional part distribution as a crude discretization check
    frac = np.abs(t - np.round(t, 6))
    print(f"  fraction of values that round cleanly to 6 decimals: "
          f"{(frac < 1e-9).mean()*100:.3f}%  (near 100% would suggest coarse rounding, not noise)")
    print(f"  sample of 10 raw hit_t values: {t[:10]}")


if __name__ == "__main__":
    main()
