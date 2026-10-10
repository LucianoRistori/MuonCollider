"""Helix validation figure: (a) direct P vs K_geom c'^(d/2); (b) counting identity per plane set."""
import csv, json, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
g = json.load(open("helix_kgeom.json"))
KG = {4: g["N4"]["K"], 6: g["0,1,2,3,4,5"]["K"]}
fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
ax = axes[0]
for N, col, mk in ((4, "#1b9e77", "o"), (6, "#d95f02", "s")):
    pts = []
    for fn in [f"helix_extrapolation_N{N}.csv"] + sorted(f for f in os.listdir(".") if f.startswith(f"helix_extrap_deep_N{N}")):
        if os.path.exists(fn):
            for r in csv.DictReader(open(fn)):
                d = int(r["ndof"]); c = float(r["c_prime"])
                Pg = KG[N] * c ** (d / 2)
                pts.append((c, float(r["P_direct"]) / Pg, float(r["P_direct_err"]) / Pg, float(r["P_direct"])))
    pts = np.array(sorted(pts))
    ax.errorbar(pts[:, 0], pts[:, 1], yerr=pts[:, 2], fmt=mk, color=col, ms=5, capsize=2,
                label=f"helix, {N} planes (ndof {2*N-4}): P from {pts[:,3].min():.0e} to {pts[:,3].max():.0e}")
# straight lines, for comparison (geometric K)
line_K = {4: 0.019038588737145246, 6: 4.624183737779514e-06}
for N, mk in ((4, "o"), (6, "s")):
    fn = f"../extrapolation_N{N}.csv"
    if os.path.exists(fn):
        rs = list(csv.DictReader(open(fn)))
        c = np.array([float(r["c_prime"]) for r in rs]); f = np.array([float(r["K"]) for r in rs]) / line_K[N]
        ax.errorbar(c, np.array([float(r["ratio"]) for r in rs]) * f, yerr=np.array([float(r["ratio_err"]) for r in rs]) * f,
                    fmt=mk, mfc="none", color="0.55", ms=4, capsize=1.5, lw=0.8,
                    label=f"straight line, {N} planes" if N == 4 else "straight line, 6 planes")
ax.axhline(1, color="gray", lw=1)
ax.set_xscale("log"); ax.set_xlim(3e-3, 10); ax.set_ylim(0.3, 1.25)
ax.set_xlabel(r"$c' = c/\lambda^2$  (the cut, mapped onto the reference detector)")
ax.set_ylabel(r"$P_\mathrm{direct}\,/\,K_\mathrm{geom}\,c'^{\,n_\mathrm{dof}/2}$")
ax.set_title("Passing probability: direct Monte Carlo vs. the small-$\\chi^2$ law\nwith K from the track-manifold volume", fontsize=10)
ax.legend(fontsize=7.5, loc="lower left"); ax.grid(True, which="both", alpha=0.25)

ax = axes[1]
import itertools, math
files = sorted(f for f in os.listdir(".") if f.startswith("helix_count_w"))
groups = {}
for fn in files:
    rs = list(csv.DictReader(open(fn)))
    groups.setdefault((float(rs[0]["w_min"]), int(rs[0]["hits"])), []).append(rs)
colors = ["#7570b3", "#e7298a", "#66a61e"]
keys = [tuple(range(6))] + list(itertools.combinations(range(6), 5))
labels = ["all 6"] + [f"w/o {6 - i}" for i in range(6)] + ["sum"]
for gi, ((wm, n), runs) in enumerate(sorted(groups.items())):
    pp = json.load(open(f"helix_pS_w{wm:g}.json"))
    ev = np.array([float(r[0]["events"]) for r in runs])
    pred = [pp[",".join(map(str, k))]["P"] * float(n) ** len(k) for k in keys]
    perr = [pp[",".join(map(str, k))]["P_err"] * float(n) ** len(k) for k in keys]
    pred.append(sum(pred)); perr.append(math.sqrt(sum(e * e for e in perr)))
    q, qe = [], []
    for i in range(8):
        c = np.average([float(r[i]["counted"]) for r in runs], weights=ev)
        ce = 1 / math.sqrt(sum(1 / max(float(r[i]["counted_err"]), 1e-9) ** 2 for r in runs))
        q.append(c / pred[i]); qe.append(math.hypot(ce / pred[i], c * perr[i] / pred[i] ** 2))
    eta = np.average([float(r[-1]["counted"]) for r in runs], weights=ev)
    x = np.arange(8) + 0.2 * (gi - 1)
    ax.errorbar(x, q, yerr=qe, fmt="o", color=colors[gi % 3], ms=5, capsize=2,
                label=f"$w_{{min}}$={wm:g}, {n} hits/plane, {int(ev.sum())} events ($\\eta$={eta:.2f}): sum ratio {q[-1]:.3f} $\\pm$ {qe[-1]:.3f}")
ax.set_xticks(np.arange(8), labels, rotation=30, fontsize=8)
ax.axhline(1, color="gray", lw=1)
ax.set_ylim(0.4, 1.6); ax.set_ylabel("counted / ($P_S \\prod n_i$)")
ax.set_title("Exact helix, noise-only events: every combination fitted\n(6 of 6 and each 5-plane subset, own 3$\\sigma$ cut)", fontsize=10)
ax.legend(fontsize=7, loc="upper left"); ax.grid(True, alpha=0.25)

ax = axes[2]
import math
# helix: combine seeds of the same configuration
cfg = {}
for fn in files:
    rs = {r["set"]: r for r in csv.DictReader(open(fn))}
    key = (float(rs["sum"]["w_min"]), int(rs["sum"]["hits"]))
    ev = int(rs["sum"]["events"])
    d = [k for k in rs if k.startswith("distinct")][0]
    cfg.setdefault(key, []).append((ev, float(rs["sum"]["counted"]), float(rs["sum"]["counted_err"]),
                                    float(rs[d]["counted"]), float(rs[d]["counted_err"]), float(rs["eta"]["counted"])))
for key, L in sorted(cfg.items()):
    w = np.array([l[0] for l in L], float)
    s_ = np.average([l[1] for l in L], weights=w); d_ = np.average([l[3] for l in L], weights=w)
    se = 1 / math.sqrt(sum(1 / l[2] ** 2 for l in L)); de = 1 / math.sqrt(sum(1 / l[4] ** 2 for l in L))
    eta = np.average([l[5] for l in L], weights=w)
    ax.errorbar(eta, d_ / s_, yerr=d_ / s_ * math.hypot(se / s_, de / d_), fmt="o", color="#d95f02", ms=5, capsize=2,
                label="exact helix" if key == sorted(cfg)[0] else None)
    ax.annotate(f"{key[1]}/plane", (eta, d_ / s_ + 0.04), fontsize=7, ha="center", color="0.4")
lr = list(csv.DictReader(open("../kofn_results.csv")))
by = {}
for r in lr:
    by.setdefault(int(r["hits"]), []).append(r)
first = True
for n, rs in sorted(by.items()):
    w = np.array([int(r["events"]) for r in rs], float)
    f = lambda k: np.average([float(r[k]) for r in rs], weights=w)
    ax.plot(f("eta"), f("distinct4") / f("sum_counted"), "s", mfc="none", color="0.5", ms=5,
            label="straight line" if first else None); first = False
x = np.linspace(0, 2.6, 100)
ax.plot(x, np.exp(-x / 2), color="k", lw=1, ls="--", label=r"$e^{-\eta/2}$")
ax.set_xlabel(r"window occupancy $\eta$"); ax.set_ylabel("distinct fakes / passing combinations")
ax.set_title("Merging candidates that share $\\geq$4 hits", fontsize=10)
ax.set_ylim(0, 1.1); ax.set_xlim(0, 2.6); ax.legend(fontsize=8); ax.grid(True, alpha=0.25)
fig.tight_layout(); fig.savefig("validation_helix.png", dpi=150); print("wrote validation_helix.png")
