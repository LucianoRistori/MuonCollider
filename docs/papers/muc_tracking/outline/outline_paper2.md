---
title: "Paper 2 (Muon Collider application) — proposed outline"
subtitle: "For discussion before writing. Draft 1, 12 October 2026"
---

# Working title

*Background rejection with hit timing and direction in a muon-collider tracker, and the resolutions it requires*

Alternatives:

- *How good must a muon-collider tracker be? Hit-level background rejection and the lower limit on fake tracks*
- *Per-hit timing and direction cuts against the beam-induced background at a 10 TeV muon collider*

**Target:** to be decided (JINST, or Eur. Phys. J. C, where the recent muon-collider detector papers are). About 20–25 pages, 10–12 figures.

**Readers:** the muon-collider detector and physics community, and tracker designers. Paper 1 is cited for the method: no derivations repeated, one paragraph of summary.

# The message in three sentences

1. With a time and a direction measured on every hit, three cuts that test whether a hit is consistent with a particle coming from the interaction point remove the beam-induced background (BIB) by a factor of about 130 in the inner and outer trackers. The cuts are on the time-of-flight-corrected time, on the z-axis intercept and on a curvature estimate. They keep about 97% of the signal hits and a track-finding efficiency of 98%.
2. After these cuts, the algorithm-independent lower limit on fake tracks (paper 1) is far below one per bunch crossing at present-day resolutions. This shows that the BIB does not make track reconstruction impossible in principle.
3. Scans of the position, time and angle resolutions show how much margin there is and which resolution matters most. They give targets for detector R&D, including what better future detectors would buy.

# Sections

## 1. Introduction (2 pp)

- The muon collider, the BIB, and its hit densities in the tracker; the tracking challenge (citations: IMCC reports, Bartosik 2020, Collamati 2021, Casarsa–Lucchesi–Sestini review, the MUSIC and MAIA concepts).
- Existing BIB mitigation: timing windows, double-layer direction matching (Di Benedetto 2018), cluster shape, and their use in full reconstruction (ACTS-based).
- What is new here:
    - Hit-level cuts using the full hit information (time and 3-D direction), with the cuts tied to signal containment rather than to one reconstruction algorithm.
    - The resulting background turned into a lower limit on fakes valid for any algorithm (paper 1).
    - The resolution requirements that follow.
- Plan of the paper.

## 2. Simulated samples and detector (2.5 pp)

- **Detector geometry:** the xml files describing the detector; the subsystems used (VXD, IT and OT, barrel and endcaps), with layers, radii and areas in a table. Which concept (MUSIC-like? to be confirmed and cited).
- **BIB sample:**
    - one bunch crossing: μ⁺ and μ⁻ beam samples plus incoherent pair production, 17.3 M hits;
    - provenance: FLUKA, then detector simulation; who produced it, and the version, are needed;
    - hit content: position, local coordinates, time, momentum (hence direction);
    - the October 2026 fix of the hit momenta.
- **Signal sample:** single-muon gun, flat in 1/p_T above 1.5 GeV/c, θ from 10° to 170°, vertex spread ±1.5 mm, 100 k events; the muon's own hits only (secondaries removed).
- **Normalization:** one BIB bunch crossing plus one signal muon per "collision".
- **Detector resolutions** as a smearing model: position (u, v), time, and hit direction (two angles), applied to both samples. The assumed values and what they correspond to physically (e.g. the direction from a doublet layer or a cluster shape) need a paragraph; see the questions at the end.

## 3. Per-hit observables (3 pp)

- **z-axis intercept z₀:** extrapolation of the hit direction to the beam axis in the meridian plane.
- **Curvature estimate:** the circle through the origin tangent to the hit's transverse direction, giving 1/R and p_T = 0.3 B R; the cut is on |1/R| (continuous at infinite p_T).
- **Time-of-flight-corrected time:** the time minus the flight time of a muon from the IP, along the helix through the hit, with the speed from the curvature estimate.
- **Distributions for BIB and signal**, per subsystem, before cuts. Figures: z₀, p_T and corrected time, full range and zoomed, with signal overlaid (per collision, signal scaled by 10ⁿ).
- **Features to describe:**
    - signal peaks sharply at the IP values;
    - the BIB is broad;
    - the BIB corrected time peaks at about +0.5 ns in the IT and OT barrels (not yet understood; Monte Carlo truth from a colleague may explain it);
    - late second crossings of low-p_T signal muons.

## 4. Selection (3 pp)

- **Cut construction:**
    - per subsystem, a symmetric window on time and z₀, plus a minimum p_T (OT only, 2 GeV/c);
    - time and z₀ are set to keep 98% of the signal hits in each cut's N−1 distribution, found together by iteration;
    - why the order matters (the late-hit population).
- **The cuts at the operating point** (table): per subsystem.
- **Results:**
    - BIB rejection factor and signal hit efficiency for p_T → ∞, per subsystem and in total: VXD 31, IT+OT 129, all 57; hit efficiency 96–98%;
    - each cut alone and on top of the other two;
    - N−1 figures.
- **Track-finding efficiency:** at least 5 surviving hits, VXD not counted. Efficiency against p_T and the extrapolation to p_T → ∞ (method in an appendix): 98.0 ± 0.2%.
- **BIB density per layer before and after the cuts** (table and figure).

## 5. Fake tracks after the cuts (2.5 pp)

- **Paper 1 in one paragraph:**
    - E = Σ_S P_S ∏ n_i, with fakes counted on at least 5 of 6 layers;
    - P_S from the volume of the accepted track manifold;
    - the lower-limit statement.
- **The IT+OT barrel tower** (six layers, projective, 1 m² outermost plane), with the densities after the cuts:
    - E per tower and bunch crossing (3.2×10⁻⁵ at 100 µm, 50 ps, 1.5°);
    - the 6-of-6 and 5-of-6 terms;
    - the window occupancy.
- **From one tower to the whole barrel:** the number of towers that tile the barrel, and the variation of the density with z and φ. *New work.*
- **Headroom against the BIB level:** one fake at 7.9× today's densities at the current operating point. Which uncertainties in the BIB simulation this covers.

## 6. Resolution requirements (4 pp)

- **Scans:**
    - one resolution at a time (position, time, angle);
    - at each point the cuts are re-derived at 98% containment and the efficiency and fake limit recomputed;
    - each scan extended both ways around the operating point: degraded up to 10×, and improved to future-detector values, e.g. 10–25 µm, 10–20 ps, 0.1–0.5° (*new runs*).
- **Results (current, for the degraded side):** one fake at 5.6× (position, 560 µm), 7.9× (time, 390 ps) and 4.6× (angle, 6.9°). Efficiency is flat by construction in the barrel tower.
- **Why the curves differ:**
    - position enters through the cuts and through P (σ⁶ for the 5-layer terms);
    - time and angle enter only through the cuts;
    - the shape of the BIB time distribution;
    - the step in the angle scan, where late hits enter the time-cut distribution;
    - the effect of the p_T cut.
- **Two-dimensional scans** (time × angle, or time × position) to show trade-offs. *New work; optional.*
- **R&D targets:** the resolutions at which the fake limit stays below, say, 10⁻² per bunch crossing for the whole barrel, with a factor of a few in margin for the BIB level.

## 7. Discussion and limitations (1.5 pp)

- **The BIB sample:** a single bunch crossing; the plus/minus consistency as a check; the statistical uncertainty on the densities.
- **What is not modelled:**
    - endcap towers;
    - clusters and merged hits;
    - inefficiency and noise of the sensors;
    - partly-fake tracks (which only add);
    - correlated BIB (curling low-momentum particles).
- **Open items:**
    - the whole-tracker efficiency drop with angle resolution;
    - the +0.5 ns BIB feature;
    - the non-monotonic gain of the p_T cut against angle resolution (minimum at 2–2.5°).
- **Cut optimization against the BIB** rather than by signal containment: open, and how much it could gain.

## 8. Conclusions (0.5 pp)

## Appendices

- A. The efficiency for p_T → ∞: the fit ε∞ + c/p_T², its justification and checks (from docs/efficiency_pt_infinity.pdf).
- B. The selection settings and every number in the tables: the configuration files of the runs used.
- C. Code and data availability: the analysis package (public), the run archive, the input samples (on request, or a public location).

# Figures (proposed)

1. r–z map of BIB hits with the signal overlaid.
2. BIB hit density per layer and disk, before and after the cuts.
3. z₀, p_T and corrected time for BIB and signal, per subsystem (one figure each, zoomed).
4. N−1 distributions with the cuts (three figures, or one combined).
5. Track-finding efficiency against p_T, with the p_T → ∞ extrapolation.
6. Expected fakes against a common density factor (barrel tower, 5 of 6).
7. Resolution scans: expected fakes against position, time and angle resolution (three panes, extended range).
8. Efficiency along the same scans.
9. Optional: a 2-D resolution map.

# New work needed

1. Decide the operating point and the reference values for the scans; rerun the analysis there (minutes).
2. Extend the three scans to improved resolutions; add finer points near "one fake".
3. Tower to whole barrel: tile the barrel, use the local (z, φ) densities after the cuts, and sum.
4. Statistical uncertainties on densities, rejection factors and E (bootstrap over BIB primaries, which the per-primary entries allow).
5. Optional:
    - the 2-D scan;
    - the endcap equivalent;
    - understand the +0.5 ns BIB peak, with the colleague's Monte Carlo truth.
6. Figures redrawn for print from the run outputs, by scripts in docs/papers/.

# Questions for you

1. **Journal:** JINST again, or EPJC?
2. **Authors:** the collaboration and authorship rules for this paper, and who produced the BIB and signal samples, to cite and acknowledge.
3. **Detector concept:** which concept the geometry files correspond to (MUSIC?), and how to cite it.
4. **Hit direction:** how should the per-hit direction measurement (1.5°) be justified physically — doublet layers, cluster shape, or a generic future capability?
5. **Operating point:** keep today's settings (100 µm, 50 ps, 1.5°), or use paper 1's conventional 100 ps and treat 50 ps as an improvement in the scans?
6. **Scope:**
    - the whole barrel (tiling) in this paper, or the tower only;
    - the endcaps left out?
