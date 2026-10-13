"""
simple_example.py - the simple case of paper 1, section 2: six parallel planes,
straight tracks (no field), uniform noise. Projective tower: planes at
Y = 250, 500, ..., 1500 mm, square, side W = Y/1.5 (1 m at 1.5 m). sigma in
both coordinates. Prints K, P and E, and draws E against the hit density for
several resolutions (simple_example.pdf/.png).
"""
import itertools, json, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
P3 = 1.3499e-3
Y = np.array([250, 500, 750, 1000, 1250, 1500.0]); W = Y / 1.5
SIG0 = 0.1


def K_line(Y, W, s, nq=400001):
    """K = V_d Vol(M)/|B| for straight lines (two views), parallel square planes."""
    N = len(Y); d = 2 * N - 4
    qmax = 1.5 * (W / Y).max(); q = np.linspace(-qmax, qmax, nq)
    hi = np.min(W / 2 - q[:, None] * Y, 1); lo = np.max(-W / 2 - q[:, None] * Y, 1)
    area = np.trapezoid(np.maximum(hi - lo, 0), q)                  # polygon of (intercept, slope)
    vol = (area * math.sqrt(N * (Y ** 2).sum() - Y.sum() ** 2) / s ** 2) ** 2
    return math.pi ** (d / 2) / math.gamma(d / 2 + 1) * vol / np.prod((W / s) ** 2), vol


sets = [tuple(range(6))] + list(itertools.combinations(range(6), 5))
P = {}
for S in sets:
    d = 2 * len(S) - 4
    K, vol = K_line(Y[list(S)], W[list(S)], SIG0)
    P[S] = dict(K=K, P=K * chi2.isf(P3, d) ** (d / 2), ndof=d, vol=vol)


def E(rho, sig=SIG0, k=5):
    n = rho * W ** 2; f = sig / SIG0
    return sum(v["P"] * f ** v["ndof"] * np.prod(n[list(S)]) for S, v in P.items() if len(S) >= k)


out = dict(Y=Y.tolist(), W=W.tolist(), sets={",".join(map(str, S)): v for S, v in P.items()})
out["one_fake_density"] = {f"{s*1000:g}um": brentq(lambda r: E(r, s) - 1, 1e-7, 10) for s in (0.025, 0.05, 0.1, 0.2)}
out["E_at_1e-3"] = {"6of6": E(1e-3, k=6), "ge5": E(1e-3)}
json.dump(out, open(HERE / "simple_example.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "sets"}, indent=1))
print("6 of 6: K %.3e  P %.3e  distinguishable tracks %.2e" % (P[sets[0]]["K"], P[sets[0]]["P"], P[sets[0]]["vol"]))

plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "dejavuserif"})
fig, ax = plt.subplots(figsize=(6.0, 4.4))
rho = np.geomspace(1e-5, 1e-1, 300)
cols = {0.025: "#1b9e77", 0.05: "#7570b3", 0.1: "#d95f02", 0.2: "#e7298a"}
for s, c in cols.items():
    ax.plot(rho, [E(r, s) for r in rho], color=c, lw=2, label=fr"$\sigma$ = {s*1000:g} $\mu$m, at least 5 of 6")
    ax.plot(rho, [E(r, s, k=6) for r in rho], color=c, lw=1, ls="--")
ax.plot([], [], color="0.4", lw=1, ls="--", label="6 of 6")
ax.axhline(1, color="0.5", lw=0.8)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-12, 1e4); ax.set_xlim(1e-5, 1e-1)
ax.set_xlabel(r"noise hit density on every plane (hits/mm$^2$)")
ax.set_ylabel(r"expected fake tracks, $E$")
ax.grid(True, which="major", alpha=0.3); ax.legend(fontsize=8, loc="lower right", frameon=False)
fig.tight_layout(); fig.savefig(HERE / "simple_example.pdf"); fig.savefig(HERE / "simple_example.png", dpi=150)
