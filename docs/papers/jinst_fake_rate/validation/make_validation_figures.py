"""Validation figure, straight lines (B = 0): the small-chi2 law with K from the
track-manifold volume (kgeom_line in helix/), and the k-of-N counting test."""
import csv, json, math, itertools, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent / "helix"))

def K_line(Y, W, s, nq=200001):
    N = len(Y); d = 2 * N - 4
    qmax = 1.5 * (W / Y).max(); q = np.linspace(-qmax, qmax, nq)
    hi = np.min(W / 2 - q[:, None] * Y, 1); lo = np.max(-W / 2 - q[:, None] * Y, 1)
    area = np.trapezoid(np.maximum(hi - lo, 0), q)
    detD = math.sqrt(N * (Y ** 2).sum() - Y.sum() ** 2) / s ** 2
    return math.pi ** (d / 2) / math.gamma(d / 2 + 1) * (area * detD) ** 2 / np.prod((W / s) ** 2)

s = 0.1
KG = {N: K_line(np.linspace(100, 100 * N, N), 2 * s * np.linspace(100, 100 * N, N) / 100, s) for N in (4, 6)}
Y6 = np.linspace(100, 600, 6); W6 = 2 * s * Y6 / 100
cal = json.load(open("calib_sets.json"))
r6 = K_line(Y6, W6, s) / cal["0,1,2,3,4,5"]["K"]
num = den = 0.0
for S in itertools.combinations(range(6), 5):
    num += K_line(Y6[list(S)], W6[list(S)], s); den += cal[",".join(map(str, S))]["K"]
r5 = num / den

plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "dejavuserif"})
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axes[0]
for N, col, mk in ((4, "#1b9e77", "o"), (6, "#d95f02", "s")):
    rows = list(csv.DictReader(open(f"extrapolation_N{N}.csv")))
    c = np.array([float(r["c_prime"]) for r in rows])
    f = np.array([float(r["K"]) for r in rows]) / KG[N]
    rt = np.array([float(r["ratio"]) for r in rows]) * f
    er = np.array([float(r["ratio_err"]) for r in rows]) * f
    Pd = [float(r["P_direct"]) for r in rows]
    ax.errorbar(c, rt, yerr=er, fmt=mk, color=col, ms=5, capsize=2,
                label=f"{N} planes ($n_{{dof}}$ = {2*N-4}): P from {min(Pd):.0e} to {max(Pd):.0e}")
ax.axhline(1, color="gray", lw=1)
ax.set_xscale("log"); ax.set_xlim(3e-3, 10); ax.set_ylim(0.4, 1.25)
ax.set_xlabel(r"$c' = c/\lambda^2$ (the cut, mapped onto the reference detector)")
ax.set_ylabel(r"$P_\mathrm{direct}\,/\,K\,c'^{\,n_\mathrm{dof}/2}$")
ax.set_title(r"Passing probability: direct Monte Carlo vs. $K\,c^{n_\mathrm{dof}/2}$, $K$ from the geometry", fontsize=10)
ax.legend(fontsize=8, loc="lower left"); ax.grid(True, which="both", alpha=0.25)

ax = axes[1]
by = {}
for r in csv.DictReader(open("kofn_results.csv")):
    by.setdefault(int(r["hits"]), []).append(r)
eta, q, qe, cm, cme = [], [], [], [], []
for n, rs in sorted(by.items()):
    w = np.array([int(r["events"]) for r in rs], float)
    f = lambda k: np.average([float(r[k]) for r in rs], weights=w)
    sm, d = f("sum_counted"), f("distinct4")
    se = 1 / math.sqrt(sum(1 / float(r["sum_counted_err"]) ** 2 for r in rs))
    de = 1 / math.sqrt(sum(1 / float(r["distinct4_err"]) ** 2 for r in rs))
    pred = f("E6_pred") * r6 + f("E5_pred") * r5
    eta.append(f("eta")); q.append(d / sm); qe.append(d / sm * math.hypot(de / d, se / sm))
    cm.append(sm / pred); cme.append(se / pred)
ax.errorbar(eta, cm, yerr=cme, fmt="s", color="#7570b3", ms=5, capsize=2, label="passing combinations, counted / predicted")
ax.errorbar(eta, q, yerr=qe, fmt="o", color="#1b9e77", ms=5, capsize=2, label="distinct fakes / passing combinations")
x = np.linspace(0, 2.4, 100)
ax.plot(x, np.exp(-x / 2), color="#1b9e77", lw=1, ls="--", label=r"$e^{-\eta/2}$")
ax.axhline(1, color="gray", lw=1)
for e_, n in zip(eta, sorted(by)):
    ax.annotate(f"{n}", (e_, 1.10), fontsize=7, ha="center", color="0.4")
ax.text(0.02, 1.16, "hits per plane:", fontsize=7, color="0.4")
ax.set_xlabel(r"window occupancy $\eta$"); ax.set_ylabel("ratio"); ax.set_ylim(0, 1.22)
ax.set_title("At least 5 of 6 planes: noise-only events, every passing combination", fontsize=10)
ax.legend(fontsize=8, loc="lower left"); ax.grid(True, alpha=0.25)
fig.tight_layout(); fig.savefig("validation_line.png", dpi=150); fig.savefig("validation_line.pdf")
print("K_geom", KG, "set factors", r6, r5, "counted/pred", np.round(cm, 3))
