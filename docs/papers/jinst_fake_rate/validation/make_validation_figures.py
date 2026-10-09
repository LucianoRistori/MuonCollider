"""Draft figures of the end-to-end validation (straight lines, B = 0)."""
import csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import json
corr = json.load(open("ref_corrections.json"))
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axes[0]
for N, col, mk in ((4, "#1b9e77", "o"), (6, "#d95f02", "s")):
    rows = list(csv.DictReader(open(f"extrapolation_N{N}.csv")))
    c = np.array([float(r["c_prime"]) for r in rows])
    rt = np.array([float(r["ratio"]) for r in rows])
    er = np.array([float(r["ratio_err"]) for r in rows])
    cc = np.geomspace(3e-3, 8, 200)
    k = corr[str(N)]
    ax.plot(cc, k["K"] / float(rows[0]["K"]) * (1 + k["beta"] * cc + k["gamma"] * cc ** 2),
            color=col, lw=1, ls="--")
    ax.errorbar(c, rt, yerr=er, fmt=mk, color=col, ms=5, capsize=2,
                label=f"{N} planes (ndof {2*N-4}): P from {min(float(r['P_direct']) for r in rows):.0e} to {max(float(r['P_direct']) for r in rows):.0e}")
ax.axhline(1, color="gray", lw=1)
ax.axvspan(0.3, 4, color="0.92", zorder=0)
ax.text(0.33, 1.19, "calibration Monte Carlo\nhas its statistics here", fontsize=7.5, color="0.35")
ax.text(0.0035, 0.66, "dashed: the power law with the\nfitted corrections $(1+\\beta c'+\\gamma c'^2)$", fontsize=7.5, color="0.35")
ax.set_xscale("log")
ax.set_xlabel(r"$c' = c/\lambda^2$  (the cut, mapped onto the calibration detector)")
ax.set_ylabel(r"$P_\mathrm{direct}\,/\,P_\mathrm{method}$")
ax.set_title("Passing probability: direct Monte Carlo vs. the method", fontsize=10.5)
ax.set_ylim(0.5, 1.25)
ax.set_xlim(3e-3, 8)
ax.legend(fontsize=8, loc="lower left")
ax.grid(True, which="both", alpha=0.25)

ax = axes[1]
rows = sorted(csv.DictReader(open("kofn_results.csv")), key=lambda r: int(r["hits"]))
# combine repeated configurations
by = {}
for r in rows:
    by.setdefault(int(r["hits"]), []).append(r)
eta, q, qe, cm, cme = [], [], [], [], []
for n, rs in sorted(by.items()):
    w = np.array([int(r["events"]) for r in rs], float)
    f = lambda k: np.average([float(r[k]) for r in rs], weights=w)
    s, d = f("sum_counted"), f("distinct4")
    se = 1 / math.sqrt(sum(1 / float(r["sum_counted_err"]) ** 2 for r in rs))
    de = 1 / math.sqrt(sum(1 / float(r["distinct4_err"]) ** 2 for r in rs))
    eta.append(f("eta")); q.append(d / s); qe.append(d / s * math.hypot(de / d, se / s))
    cm.append(s / f("sum_pred")); cme.append(se / f("sum_pred"))
eta = np.array(eta)
ax.errorbar(eta, cm, yerr=cme, fmt="s", color="#7570b3", ms=5, capsize=2,
            label="passing combinations, counted / method")
ax.errorbar(eta, q, yerr=qe, fmt="o", color="#1b9e77", ms=5, capsize=2,
            label="distinct fakes / passing combinations")
x = np.linspace(0, 2.4, 100)
ax.plot(x, np.exp(-x / 2), color="#1b9e77", lw=1, ls="--", label=r"$e^{-\eta/2}$")
ax.axhline(1, color="gray", lw=1)
ax.set_xlabel(r"window occupancy $\eta$ (other hits that could replace one of the candidate's)")
ax.set_ylabel("ratio")
ax.set_title("At least 5 of 6 planes: noise-only events, every passing combination", fontsize=10.5)
ax.set_ylim(0, 1.2)
ax.legend(fontsize=8, loc="lower left")
ax.grid(True, alpha=0.25)
for e_, n in zip(eta, sorted(by)):
    ax.annotate(f"{n}", (e_, 1.08 if n != 400 else 1.13), fontsize=7, ha="center", color="0.4")
ax.text(0.02, 1.15, "hits per plane:", fontsize=7, color="0.4")
fig.tight_layout()
fig.savefig("validation_line.png", dpi=150)
print("wrote validation_line.png")
