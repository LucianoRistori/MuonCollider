import json, glob, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "dejavuserif"})
rows = {}
for f in glob.glob("v_*.log"):
    for l in open(f):
        if l.startswith("{"):
            r = json.loads(l); rows.setdefault(r["geom"], {})[r["L"]] = r
lab = {"cylinders": "4 coaxial cylinders (barrel-like)", "tilted45z": "4 planes tilted 45° (transverse)",
       "tilted60_endcaplike": "4 planes tilted 60° towards the beam", "mixed": "2 cylinders + 2 endcap disks"}
fig, ax = plt.subplots(figsize=(7, 4.6))
for (g, mk, col) in (("cylinders", "o", "#1b9e77"), ("tilted45z", "s", "#d95f02"), ("tilted60_endcaplike", "^", "#7570b3"), ("mixed", "D", "#e7298a")):
    rs = [rows[g][L] for L in sorted(rows[g])]
    ax.errorbar([r["P_tube"] for r in rs], [r["ratio"] for r in rs], yerr=[r["ratio_err"] for r in rs],
                fmt=mk + "-", color=col, ms=5, capsize=2, lw=0.8, label=lab[g])
ax.axhline(1, color="gray", lw=1); ax.set_xscale("log"); ax.invert_xaxis()
ax.set_xlabel("passing probability from the formula (detectors grow to the right)")
ax.set_ylabel(r"$P_\mathrm{direct}\,/\,P_\mathrm{formula}$"); ax.set_ylim(0.2, 1.35)
ax.set_title("General surfaces: direct Monte Carlo vs. the track-manifold formula", fontsize=10)
ax.legend(fontsize=8, loc="lower right"); ax.grid(True, which="both", alpha=0.25)
fig.tight_layout(); fig.savefig("validation_general_surfaces.png", dpi=150); fig.savefig("validation_general_surfaces.pdf")
