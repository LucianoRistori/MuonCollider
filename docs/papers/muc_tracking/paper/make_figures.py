"""make_figures.py - figures of paper 2 from the outputs of run 2026-10-14_162239_pdf
(100 um, 100 ps, 1.5 deg) and of the resolution scans around it.
Inputs (from claude_scratch/op100 of the Analysis folder): hists_op100.npz,
density_per_layer_before_after_cuts.csv, track_efficiency_vs_pt.csv,
track_efficiency_limit.csv, fake_rate.csv, scan_<time|position|angle>_results.csv,
and the whole-barrel results barrel_z0_4.5_s1.json, barrel_z0_1500.0_s1.json.
Usage: python3 make_figures.py <input folder>"""
import csv, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.size": 10,
                     "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True, "legend.frameon": False})
IN = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
BIB, SIG, GREY = "#1f4e8c", "#c0392b", "0.45"
CUTS = {"IT": dict(t=0.233, z0=52, pt=None), "OT": dict(t=0.245, z0=137, pt=2.0),
        "ITE": dict(t=0.242, z0=474, pt=None), "OTE": dict(t=0.241, z0=349, pt=2.0)}
NAMES = {"IT": "IT barrel", "OT": "OT barrel", "ITE": "IT endcap", "OTE": "OT endcap"}


def rows(name):
    return list(csv.DictReader(open(IN / name)))


def save(fig, name):
    fig.savefig(f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{name}.png", dpi=180, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


# ---------------------------------------------------------------- 1. observables (N-1)
H = np.load(IN / "hists_op100.npz")
_rc = plt.rcParams.copy(); plt.rcParams.update({"font.size": 8})
fig, axs = plt.subplots(4, 3, figsize=(7.2, 8.6))
ylabels = {"t": r"$\frac{1}{N}\,\mathrm{d}N/\mathrm{d}t_\mathrm{corr}$ (ns$^{-1}$)",
           "z0": r"$\frac{1}{N}\,\mathrm{d}N/\mathrm{d}z_0$ (mm$^{-1}$)",
           "k": r"$\frac{1}{N}\,\mathrm{d}N/\mathrm{d}(1/R)$ (m)"}
labels = {"t": r"$t_\mathrm{corr}$ (ns)", "z0": r"$z_0$ (mm)", "k": r"curvature $1/R$ (m$^{-1}$)"}
for r, sub in enumerate(("IT", "OT", "ITE", "OTE")):
    for c, var in enumerate(("t", "z0", "k")):
        ax = axs[r, c]; e = H[f"edges_{sub}_{var}"]; x = 0.5 * (e[1:] + e[:-1]); w = np.diff(e)
        for name, col, lab in (("bib", BIB, "BIB"), ("sig", SIG, "signal")):
            h = H[f"{name}_{sub}_{var}_n1"].astype(float)
            ax.step(x, h / h.sum() / w, where="mid", color=col, lw=0.9, label=lab)
        ax.set_yscale("log")
        cut = CUTS[sub]
        if var == "t":
            lim = cut["t"]
        elif var == "z0":
            lim = cut["z0"]
        else:
            lim = 0.3 * 5.0 / cut["pt"] if cut["pt"] else None
        if lim is not None:
            for s in (-lim, lim):
                ax.axvline(s, color="0.2", lw=0.7, ls=(0, (3, 2)))
        ax.set_xlim(e[0], e[-1])
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi * 25)               # headroom for the panel label
        ax.set_xlabel(labels[var], labelpad=1)
        ax.text(0.03, 0.95, NAMES[sub], transform=ax.transAxes, ha="left", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
        if r == 0 and c == 0:
            ax.legend(loc="upper right", fontsize=7)
fig.tight_layout(h_pad=0.5, w_pad=0.6)
save(fig, "fig_observables")
plt.rcParams.update(_rc)

# ---------------------------------------------------------------- 2. densities per layer
dd = rows("density_per_layer_before_after_cuts.csv")
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, systems, xkey, xl in ((axs[0], (("3", "IT"), ("5", "OT")), "mean_r_mm", "barrel layer radius (mm)"),
                              (axs[1], (("4", "IT"), ("6", "OT")), "mean_abs_z_mm", r"endcap disk $|z|$ (mm)")):
    for sysid, nm in systems:
        d = [r for r in dd if r["system"] == sysid]
        X = [float(r[xkey]) for r in d]
        mk = "o" if nm == "IT" else "s"
        ax.plot(X, [float(r["bib_mean_density_before"]) for r in d], mk + "-", color=GREY, lw=1.0, ms=4,
                label=f"{nm}, before the cuts")
        ax.plot(X, [float(r["bib_mean_density_after"]) for r in d], mk + "-", color=BIB, lw=1.0, ms=4,
                label=f"{nm}, after the cuts")
    ax.set_yscale("log"); ax.set_ylim(5e-5, 5); ax.set_xlabel(xl)
    ax.legend(loc="upper right", fontsize=7.5)
axs[0].set_ylabel(r"BIB hit density (mm$^{-2}$)")
fig.tight_layout(); save(fig, "fig_density")

# ---------------------------------------------------------------- 3. efficiency vs pT
e = rows("track_efficiency_vs_pt.csv")
lim = rows("track_efficiency_limit.csv")[0]
x = np.array([float(r["inv_pT_center_per_gev"]) for r in e])
y = np.array([float(r["efficiency"]) for r in e]); ey = np.array([float(r["efficiency_unc"]) for r in e])
fig, ax = plt.subplots(figsize=(4.6, 3.1))
ax.errorbar(x, 100 * y, yerr=100 * ey, fmt="o", ms=2.5, color=BIB, lw=0.7, capsize=0)
einf, c2 = float(lim["efficiency_pt_inf"]), float(lim["fit_c_gev2"])
xx = np.linspace(0, 0.1, 50)
ax.plot(xx, 100 * (einf + c2 * xx ** 2), color=SIG, lw=1.2,
        label=rf"$\varepsilon_\infty + c/p_T^2$, $\varepsilon_\infty = {100*einf:.1f}\%$")
ax.set_xlabel(r"$1/p_T$ (GeV$^{-1}c$)"); ax.set_ylabel("track-finding efficiency (%)")
ax.set_xlim(0, x.max() * 1.02)
sec = ax.secondary_xaxis("top", functions=(lambda v: v, lambda v: v))
ticks = [2, 3, 5, 10, 20]
sec.set_xticks([1 / t for t in ticks]); sec.set_xticklabels([str(t) for t in ticks])
sec.set_xlabel(r"$p_T$ (GeV/$c$)")
ax.legend(loc="lower left", fontsize=8.5)
fig.tight_layout(); save(fig, "fig_efficiency")

# ---------------------------------------------------------------- 4. E vs density factor
fr = rows("fake_rate.csv")[0]
E5, E6 = float(fr["efakes_5_layers"]), float(fr["efakes_6_layers"])
B = json.load(open(IN / "barrel_z0_4.5_s1.json"))
B5, B6 = B["by_m"]["5"], B["by_m"]["6"]
RG = json.load(open(IN / "regions_op100_z4.5_s1.json"))
T5, T6 = RG["by_n_layers"]["5"], RG["by_n_layers"]["6"]
lam = np.geomspace(0.1, 10, 200)
fig, ax = plt.subplots(figsize=(4.6, 3.2))
ax.plot(lam, T5 * lam ** 5 + T6 * lam ** 6, color="k", lw=1.4, label=r"IT + OT, $|z_0|<4.5$ mm")
ax.plot(lam, B5 * lam ** 5 + B6 * lam ** 6, color=SIG, lw=1.2, ls=(0, (5, 2)), label=r"barrel only")
ax.plot(lam, E5 * lam ** 5 + E6 * lam ** 6, color=BIB, lw=1.4, label="one projective tower")
ax.plot(lam, E6 * lam ** 6, color=BIB, lw=0.9, ls=(0, (3, 2)), label="tower, 6 of 6 only")
ax.axhline(1, color="0.3", lw=0.7); ax.axvline(1, color="0.6", lw=0.7, ls=":")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-8, 1e3)
ax.set_xlabel("BIB density / simulated density"); ax.set_ylabel("expected fakes per bunch crossing")
ax.legend(loc="upper left", fontsize=8.5)
fig.tight_layout(); save(fig, "fig_density_factor")

# ---------------------------------------------------------------- 5. fakes vs polar angle, three regions
fig, ax = plt.subplots(figsize=(4.6, 3.1))
tb = np.array(RG["t_bins"]); ang = np.degrees(np.arctan2(1.0, tb))       # cot(theta) -> theta (deg)
for reg, col, ls in (("barrel", SIG, "-"), ("intermediate", BIB, "-"), ("endcap", "0.35", "-")):
    h = np.array(RG[reg]["t_hist"])
    dens = h / np.abs(np.diff(ang))                                       # per degree, both ends together
    ax.step(np.r_[ang[:-1], ang[-1]], np.r_[dens, dens[-1]], where="post", color=col, lw=1.2, ls=ls,
            label=f"{reg} ({RG[reg]['E']:.1e})")
ax.set_yscale("log"); ax.set_ylim(1e-10, 1e-2); ax.set_xlim(5, 90)
ax.set_xlabel(r"polar angle $\vartheta$ of the fake track (deg, both ends folded)")
ax.set_ylabel(r"$\mathrm{d}E/\mathrm{d}\vartheta$ per BX (deg$^{-1}$)")
ax.legend(loc="upper left", fontsize=8)
fig.tight_layout(); save(fig, "fig_polar")

# ---------------------------------------------------------------- 6. resolution scans
OP = {"position": 0.1, "time": 0.1, "angle": 1.5}
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6), sharey="row")
for c, (p, col, unit, scale, xl) in enumerate((("position", "sigma_u_mm", "µm", 1000, r"position resolution ($\mu$m)"),
                                                ("time", "sigma_t_ns", "ps", 1000, "time resolution (ps)"),
                                                ("angle", "sigma_angle_deg", "deg", 1, r"angle resolution ($^\circ$)"))):
    rr = rows(f"scan_{p}_results.csv")
    x = np.array([float(r[col]) for r in rr]) * scale
    E = np.array([float(r["efakes_exact_helix"]) for r in rr])
    eff = np.array([float(r["eff_inf"]) for r in rr]); de = np.array([float(r["eff_inf_unc"]) for r in rr])
    ax = axs[0, c]
    ax.plot(x, E, "o-", color=BIB, ms=3, lw=1.0, label="one tower")
    ax.plot(x, E * B["E"] / E5 if False else E * (B["E"] / (E5 + E6)), "s-", color=SIG, ms=2.5, lw=1.0,
            label="whole barrel (scaled)")
    ax.axhline(1, color="0.3", lw=0.7); ax.axvline(OP[p] * scale, color="0.6", lw=0.7, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-9, 1e4)
    if c == 0:
        ax.set_ylabel("expected fakes per BX"); ax.legend(loc="upper left", fontsize=7.5)
    ax = axs[1, c]
    ax.errorbar(x, 100 * eff, yerr=100 * de, fmt="o-", ms=3, lw=0.9, color=BIB, capsize=0)
    ax.axvline(OP[p] * scale, color="0.6", lw=0.7, ls=":")
    ax.set_xscale("log"); ax.set_ylim(95, 100); ax.set_xlabel(xl)
    if c == 0:
        ax.set_ylabel(r"efficiency, $p_T\to\infty$ (%)")
fig.tight_layout(h_pad=0.6, w_pad=0.6)
save(fig, "fig_scans")
print("done")
