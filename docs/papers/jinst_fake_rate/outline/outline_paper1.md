---
title: "Paper 1 (JINST) — proposed outline"
subtitle: "For approval before writing. Draft 2, 9 October 2026"
---

# Title

*An algorithm-independent lower limit on combinatorial fake tracks in high-occupancy tracking detectors*

**Target:** JINST regular article, about 20 pages, 7–8 figures. **Readers:** physicists who know track reconstruction at large collider experiments; no tutorial material on fits or Kalman filters.

# The message in three sentences

1. When noise hits become dense enough compared with the detector resolution, some random combinations of hits fit a track as well as real tracks do; no algorithm that keeps real tracks at the same efficiency can reject them, because the information is not in the data.
2. Their expected number is exactly E = P × ∏ n~i~ (P: probability that one random combination passes the χ² cut; n~i~: noise hits per plane), and P, which can be 10^-26^, is computed with an exact rescaling of the detector plus ordinary Monte Carlo.
3. E is therefore a lower limit on the fake rate of any algorithm; we show how it scales with hit density and resolution and apply it to a six-layer barrel tower with Muon Collider–like backgrounds.

# Sections

## 1. Introduction (1.5 pp)

- High-occupancy tracking: HL-LHC, heavy ions, the Muon Collider beam-induced background. Fake rates are normally obtained by full simulation of one specific algorithm.
- The question asked here: at what density does tracking become impossible *for any algorithm*?
- Precursor: Casarsa, Jindariani, Ristori, JINST 20 (2025) P04030, Sec. 8.5 (inflated-resolution Monte Carlo scaled back with a product-of-resolutions law). This paper makes that argument exact and general.
- Related ideas outside tracking, cited and contrasted: random alignments of points (Kendall & Kendall; Broadbent), the "number of false alarms" of a contrario detection (Desolneux, Moisan, Morel), rare-event methods that inflate a scale and extrapolate (Bucher; Naess et al.). Accidental coincidences in counter systems and fake roads in associative-memory triggers as the familiar product-of-occupancies heuristics.
- Plan of the paper.

## 2. Definitions and the lower-limit statement (2 pp)

- Detector model: N planes, noise hits uniformly distributed with density ρ~i~ on plane i, measurement resolution σ (per coordinate, per plane).
- Track hypothesis and acceptance: least-squares fit, χ² < cut, with the cut fixed by the efficiency for genuine tracks (one-sided 3σ: P(χ² > cut) = 1.35×10^-3^).
- The counting identity E = P ∏ n~i~: exact for the mean by linearity of expectation, whatever the correlations between overlapping combinations.
- **The lower-limit statement** (wording as agreed): any algorithm that uses only the hit positions and keeps genuine tracks passing this χ² cut must also accept, on average, at least E pure-noise combinations. State the class of algorithms it covers; what can lower it (an algorithm that drops candidates through hit sharing or ambiguity resolution also loses real tracks) and what it does not count (partly-fake tracks, which only add to it).
- k of N planes: E~k-of-N~ = Σ over subsets of k planes; bounds (largest subset ≤ distinct fakes ≤ sum).

## 3. The probability that one random combination passes (4 pp)

Developed throughout on the B = 0 straight line, the simplest case, which shows the whole method; Section 5 adds the field.

- **Whitened coordinates from the start:** every measured coordinate divided by its resolution. The fit becomes an unweighted least-squares fit, χ² is the squared length of the residual vector, and χ² < c means that this vector lies inside a **sphere** of radius √c in the ndof-dimensional residual space. Planes of different size and different resolutions need no separate treatment: in these coordinates each plane is simply a box of side W~i~/σ~i~.
- Small-χ² power law P(χ² < c) ≈ K c^ndof/2^: the density of the residual vector near zero is constant, so P is that density times the volume of the sphere. The straight-line case reproduces the (w/L)^k−2^ alignment law.
- Exact rescaling: in whitened coordinates shrinking the planes by λ (or, equivalently, enlarging every σ by λ) multiplies χ² by λ^2^ exactly, since the residuals are linear in the data; so P~true~(c) = P~shrunk~(c/λ^2^), and the rare event becomes common in the shrunken detector.
- Calibration: ordinary Monte Carlo on the shrunken detector; fit K on the plateau of P/c^ndof/2^ with the exponent fixed; how to check the plateau; choice of λ (smallest plane 5–20 resolutions wide).
- Consequence: P ∝ σ^ndof^, the resolution scaling, is the same symmetry read the other way.
- Relation to rare-event methods and to extreme-value theory: here the scaling is an exact symmetry and the exponent is known, so only K is fitted.

## 4. Validation (2.5 pp)

- The plateau of P/c^ndof/2^ over several decades of c, for the straight line and the helix (figure).
- The resolution law checked numerically (the 10^ndof/N^ ratio of crossing densities).
- **New:** an end-to-end test in small, tractable detectors where fakes are frequent enough to count directly: generate noise-only events, fit *all* combinations, count those passing the cut, and compare with P ∏ n~i~ — for the straight line and for the helix, for 6 of 6 and k of N (figure).

## 5. Tracks in a solenoidal field (3.5 pp)

- Tracks from the beam line: the fit splits into a bending view (circle through the beam line, parameters b, κ) and a depth view (a line), ndof = 2N − 4.
- The exact helix fit (Levenberg–Marquardt, several starting points; one sentence on why a single start would underestimate fakes).
- Scale equivariance of the circle model: the rescaling stays exact provided the curvature limit scales with λ (κ → λκ).
- Curvature acceptance (the p~T~ threshold): the accepted fraction grows linearly with κ~max~ and does not depend on σ.
- Comparison with the B = 0 line of Section 3: what the field costs and gains (one more fitted parameter, but a curvature acceptance).

## 6. Worked example: a six-layer barrel tower (3 pp)

- Geometry and densities *inspired by* a Muon Collider tracker (citations): six layers at radii ~160–1500 mm, projective tower with 1 m^2^ outermost plane, densities after background-suppression cuts of 10^-4^–10^-2^ hits/mm^2^, all rounded to two significant figures; σ = 50 µm; B = 5 T; p~T~ > 10 GeV/c.
- Results: P for the exact helix; E for 6 of 6 and for 5 of 6 (six subsets, compared).
- Headroom: E vs a common density factor (∝ μ^6^) and vs resolution (∝ σ^8^), the factor at which one fake is reached.

## 7. Uncertainties and limitations (1.5 pp)

- Statistical uncertainty of K; systematic from the plateau range; fit convergence.
- Assumptions: Gaussian resolution without tails, multiple scattering neglected (high-p~T~ tracks), projective idealization, noise uniform within a plane (local density is what enters).
- What would raise the true fake rate above the limit: partly-fake tracks, correlated noise (clusters, loopers).

## 8. Conclusions (0.5 pp)

## Appendices

- A. Notation.
- B. Calibration constants (K, ndof, λ, statistics) for every configuration used.
- C. Code and data availability (git repository; archived version).
- Acknowledgements, including the AI disclosure statement required by IOP.

# What is left out of the compendium

The history of the work and the "pitfalls" narratives (one sentence each where useful); the discussion on choosing λ (one rule of thumb); the conservative quadratic model; the angular-acceptance and impact-parameter cuts; the first-order validity study (one remark); everything that refers to simulation data, which goes to paper 2.

# New work needed before writing results

1. Recalibrate on the rounded geometry: exact helix for 6 of 6 and for each of the six 5-plane subsets (the slow runs; time to be estimated first).
2. Uncertainties: statistical (bootstrap) and plateau-range systematic for every K.
3. The end-to-end validation of Section 4 (new Monte Carlo, small detectors; straight line and helix; 6 of 6 and k of N).
4. All figures regenerated by scripts in docs/papers/jinst_fake_rate/.

# Decided (9 October 2026)

1. Title: the one above.
2. Whitened coordinates from the start; spheres, not ellipsoids.
3. Only the exact helix fit; the B = 0 line as the pedagogical case. No conservative quadratic model.
4. No angular-acceptance or impact-parameter cuts.
5. The end-to-end validation in small, tractable detectors is included.
