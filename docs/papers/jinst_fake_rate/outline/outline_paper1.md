---
title: "Paper 1 (JINST) — proposed outline"
subtitle: "For approval before writing. Draft 4, 10 October 2026"
---

# What changed since draft 3

- **K is now computed from the geometry**, not fitted to a Monte Carlo (Section 3). K is the volume of the set of accepted tracks in whitened hit space, times the volume of the unit sphere, divided by the volume of the hit space. This is exact in the small-χ² limit. The Monte Carlo is kept, but only to measure the small finite-cut corrections and as a check. Reason: the end-to-end validation showed that, for the helix, a fit of K to calibration Monte Carlo is reliable only to about ±10%, while the geometric value is exact to 10^-4^ and agrees with the direct Monte Carlo to 2–3% down to P = 3.5×10^-7^.
- **Subsection 2.x:** the bound "largest term ≤ distinct fakes ≤ sum" was wrong (a single subset can also contain overlapping candidates). It is replaced by the window-occupancy relation, distinct ≈ sum × e^−η/2^, measured in the validation.
- **Section 4** now describes the validation as done, with its results.
- **Literature:** the 1994 ATLAS notes are dropped (internal, not publicly accessible, not citable). The volume-of-tubes literature (Hotelling 1939; Weyl 1939; Johansen & Johnstone 1990) is added as the mathematical root of the geometric K.

# Title

*An algorithm-independent lower limit on combinatorial fake tracks in high-occupancy tracking detectors*

**Target:** JINST regular article, about 20 pages, 7–8 figures. **Readers:** physicists who know track reconstruction at large collider experiments; no tutorial material on fits or Kalman filters.

# The message in three sentences

1. When noise hits become dense enough compared with the detector resolution, some random combinations of hits fit a track as well as real tracks do; no algorithm that keeps real tracks at the same efficiency can reject them, because the information is not in the data.
2. Their expected number is exactly E = P × ∏ n~i~ (P: probability that one random combination passes the χ² cut; n~i~: noise hits per plane), and P, which can be 10^-26^, follows from the geometry: it is the number of distinguishable tracks in the acceptance times the volume of a χ² sphere, relative to the volume of the hit space.
3. E is therefore a lower limit on the fake rate of any algorithm; we show how it scales with hit density and resolution and apply it to a six-layer barrel tower with Muon Collider–like backgrounds.

# Sections

## 1. Introduction (1.5 pp)

- High-occupancy tracking: HL-LHC, heavy ions, the Muon Collider beam-induced background. Fake rates are normally obtained by full simulation of one specific algorithm.
- The question asked here: at what density does tracking become impossible *for any algorithm*?
- Precursor: Casarsa, Jindariani, Ristori, JINST 20 (2025) P04030, Sec. 8.5 (inflated-resolution Monte Carlo scaled back with a product-of-resolutions law). This paper makes that argument exact and general.
- Related ideas outside tracking, cited and contrasted: random alignments of points (Kendall & Kendall; Broadbent); the "number of false alarms" of a contrario detection (Desolneux, Moisan, Morel); the volume of tubes (Hotelling; Weyl), used in statistics for significance tests on nonlinear models; rare-event methods that inflate a scale and extrapolate (Bucher; Naess et al.). Accidental coincidences in counter systems and fake roads in associative-memory triggers as the familiar product-of-occupancies heuristics.
- Plan of the paper.

## 2. Definitions and the lower-limit statement (2 pp)

- Detector model: N planes, noise hits uniformly distributed with density ρ~i~ on plane i, measurement resolution σ (per coordinate, per plane).
- Track hypothesis and acceptance: least-squares fit, χ² < cut, with the cut fixed by the efficiency for genuine tracks (one-sided 3σ: P(χ² > cut) = 1.35×10^-3^).
- The counting identity E = P ∏ n~i~: exact for the mean by linearity of expectation, whatever the correlations between overlapping combinations.
- **The lower-limit statement** (wording as agreed): any algorithm that uses only the hit positions and keeps genuine tracks passing this χ² cut must also accept, on average, at least E pure-noise combinations. State the class of algorithms it covers; what can lower it (an algorithm that drops candidates through hit sharing or ambiguity resolution also loses real tracks) and what it does not count (partly-fake tracks, which only add to it).

### 2.x Requiring at least k of N planes

- **Acceptance rule:** a candidate is accepted if its hits on some set S of at least k planes pass the χ² cut for that set's ndof; each ndof has its own cut, at the same efficiency for genuine tracks (one-sided 3σ). For 5 of 6: ndof = 8 for the full set, 6 for each 5-plane subset.
- **Expected number of accepted combinations**, by the same linearity argument: E(≥ k of N) = Σ P~S~ ∏~i∈S~ n~i~, the sum running over the sets S of at least k planes — for 5 of 6, E~6~ plus six terms E~5~(S). Each P~S~ is computed for its own geometry, since removing a plane changes the lever arm and the acceptance.
- **From combinations to distinct fake tracks:** the sum counts some fakes more than once — a 6-hit fake usually contains 5-hit subsets that also pass, and two 5-hit fakes can share 4 hits; any algorithm would merge them. What controls the double counting is the **window occupancy** η: for an accepted candidate, the number of other hits that could replace one of its hits and still pass, summed over its planes. For a plane i it is η~i~ ≈ (n~i~ − 1) π (c/2) (σ~i~^2^ + V~i~)/A~i~, with V~i~ the variance of the fit's prediction on that plane and A~i~ the plane area, so η is known from the geometry.
- **The relation:** distinct fakes ≈ sum × e^−η/2^. It is exact as η → 0, where the sum itself is the answer.
- **Measured, not just argued (Section 4):** e^−η/2^ holds within 2–5% for straight lines over 0.17 < η < 2.2. For the helix, the measured ratio falls below it by a fraction of about 0.1 η: at η = 2.4 it is 0.23 against 0.30. In the regime where E matters for the lower limit, η much smaller than 1, the correction is a few per cent at most, and the statement quotes the sum together with the occupancy correction.

## 3. The probability that one random combination passes (4 pp)

Developed throughout on the B = 0 straight line, the simplest case, which shows the whole method; Section 5 adds the field.

- **Whitened coordinates from the start:** every measured coordinate is divided by its resolution. The fit becomes an unweighted least-squares fit, χ² is the squared distance from the measurement point to the set of points a track can produce, and χ² < c means that the point lies within distance √c of that set. Planes of different size and different resolutions need no separate treatment: in these coordinates each plane is simply a box of side W~i~/σ~i~.
- **The small-χ² power law, derived for readers who have not met it** (about one page):
    - Each track with parameters θ (p of them) produces an ideal point x(θ) in the m-dimensional whitened hit space. As θ runs over the accepted tracks, these points form a p-dimensional surface, the **track manifold** M. χ² of a combination is its squared distance from M; ndof = m − p.
    - The combinations with χ² < c fill a thin "tube" of radius √c around M. Locally, in the ndof directions perpendicular to M, the tube is a **sphere** of radius √c. Its volume is therefore the area of M times the sphere volume V~ndof~ c^ndof/2^, with V~d~ = π^d/2^/Γ(d/2 + 1). This is the volume-of-tubes result of Hotelling and Weyl (1939), to leading order in c.
    - Random noise hits are uniform in the box, whose whitened volume is ∏ (W~i~/σ~i~)^2^. Hence **P ≈ K c^ndof/2^, with K = V~ndof~ × Vol(M) / ∏ (W~i~/σ~i~)^2^.** The exponent comes from geometry alone, and K is computed as well, not fitted.
    - Vol(M) = ∫ √det(J^T^J) dθ over the accepted track parameters (J: derivatives of the whitened ideal hits with respect to θ). For straight lines it reduces to the area of a polygon in (intercept, slope) times a constant, one per view. Physical reading: **Vol(M) is the number of distinguishable tracks in the acceptance**, counted in units of the resolution. E is then (number of distinguishable tracks) × (probability that random hits fall within the χ² sphere of one of them).
    - The leading corrections are of relative order c, from the edges of the acceptance and the curvature of M: P = K c^ndof/2^ (1 + β c + …). They matter only when the sphere is not small compared with the planes.
    - The same power appears in the familiar small-argument form of the χ² distribution, P ≈ (c/2)^ndof/2^/Γ(ndof/2 + 1) for Gaussian residuals (incomplete-gamma series; DLMF §8.7), with a different constant because here the hits are uniform noise. For straight lines it reproduces the (w/L)^k−2^ law for chance alignments of k random points (Kendall & Kendall 1980; Broadbent 1980).
- **Exact rescaling:** shrinking the planes by λ (equivalently, enlarging every σ by λ) multiplies χ² by λ^2^ exactly, so P~true~(c) = P~shrunk~(c/λ^2^). The shrunken detector makes the event common enough to measure P directly with ordinary Monte Carlo. Here it serves to measure the corrections β and to check the geometric K. In any realistic configuration the cut maps to c′ = c/λ^2^ much smaller than 1, where the corrections are negligible.
- **Consequence:** P ∝ σ^ndof^, the resolution scaling, is the same symmetry read the other way, and also follows directly from the formula for K.
- **Relation to rare-event methods:** the scale is an exact symmetry and both the exponent and the constant are known, so nothing is extrapolated by a fit.

## 4. Validation (2.5 pp)

Small, tractable towers (4 and 6 planes, σ = 0.1 mm), for the straight line and for the exact helix (two figures, already produced by scripts in the repository).

- **The small-χ² law with the geometric K:** direct Monte Carlo P divided by K c′^ndof/2^ tends to 1 as c′ → 0.
    - Straight line: within ±5% down to P ~ 10^-7^.
    - Helix, 4 planes: 0.97, 0.99, 1.02, 1.01 (±0.02–0.03) at c′ = 0.07–0.008, P down to 3.5×10^-7^.
    - Helix, 6 planes: P down to 4×10^-8^, consistent with the same 1 + β c′ trend.
    - The figure also shows the finite-c′ corrections.
- **Counting identity, by brute force:** noise-only events, every combination fitted (6 of 6 and the six 5-plane subsets, each with its own cut), counted against P~S~ ∏ n~i~.
    - Straight line: 0.990 ± 0.008 with the geometric K (1.006 ± 0.008 with the fitted one) over six configurations, deep in the small-χ² regime (P~S~ down to 10^-15^, using exact pruning of the combinatorics).
    - Helix: 0.97 ± 0.03 over three configurations.
- **Distinct fakes against window occupancy** (figure): e^−η/2^, as discussed in 2.x.
- **The geometric K against the Monte Carlo calibration:** straight lines agree within 1–4%. For the helix the calibration fit depends on its range and order by ±10%, which is why the paper uses the geometric value.

## 5. Tracks in a solenoidal field (3.5 pp)

- Tracks from the beam line: the fit splits into a bending view (circle through the beam line, parameters b, κ) and a depth view (a line in the path length), ndof = 2N − 4.
- The exact helix fit (Levenberg–Marquardt, several starting points; one sentence on why a single start would underestimate fakes).
- **K for the helix:** Vol(M) is a two-dimensional integral over (b, κ) of the bending-view Jacobian, times the depth-view acceptance area (a polygon), computed on a grid to 10^-4^ in about a minute.
    - The depth fit uses the measured, not the fitted, bending coordinate. This only shears the space and does not change the volume, so the factorization still holds; one paragraph.
- **Scale equivariance of the circle model:** the rescaling stays exact provided the curvature limit scales with λ (κ → λκ).
- **Curvature acceptance (the p~T~ threshold):** it enters K through the κ range of the integral. The accepted fraction grows linearly with κ~max~ and does not depend on σ.
- **Comparison with the B = 0 line of Section 3:** what the field costs and gains (one more fitted parameter, but a curvature acceptance).

## 6. Worked example: a six-layer barrel tower (3 pp)

- **Geometry and densities *inspired by* a Muon Collider tracker** (citations):
    - six layers at radii ~160–1500 mm, projective tower with a 1 m^2^ outermost plane;
    - densities after background-suppression cuts of 10^-4^–10^-2^ hits/mm^2^;
    - all rounded to two significant figures;
    - σ = 100 µm (both coordinates) and 100 ps time resolution: conventional, present-day assumptions, suited to a paper about a general method. The time resolution enters only through the hit densities left after the timing cut. B = 5 T; p~T~ > 10 GeV/c.
- **Results:** K from the geometry and P for the exact helix; E for 6 of 6 and for 5 of 6 (six subsets, compared); η, and the distinct-fake correction. The 5-of-6 terms dominate by about four orders of magnitude, so the k-of-N treatment is essential, not a refinement.
- **Headroom:** E vs a common density factor (∝ μ^6^ for 6 of 6, dominated by μ^5^ terms for 5 of 6) and vs resolution (∝ σ^ndof^), and the factor at which one fake is reached.
- Futuristic detector assumptions (better resolutions for a detector 30–50 years from now) and the resolution scans that could set goals for detector R&D are left to paper 2.

## 7. Uncertainties and limitations (1.5 pp)

- **Uncertainty of K:** numerical integration (10^-4^); finite-c corrections at the example's c′ (bounded with the measured β); fit convergence of the multi-start helix fit.
- **Assumptions:**
    - Gaussian resolution without tails;
    - multiple scattering neglected (high-p~T~ tracks);
    - projective idealization;
    - noise uniform within a plane: the local density is what enters, and a non-uniform density enters K as a weight in the same integral.
- **What would raise the true fake rate above the limit:** partly-fake tracks; correlated noise (clusters, loopers).

## 8. Conclusions (0.5 pp)

## Appendices

- A. Notation.
- B. The track-manifold volume: derivation of the tube formula to leading order; the straight-line and helix integrals; the shear argument for the depth view.
- C. Constants (K, ndof, the correction β, η) for every configuration used.
- D. Code and data availability (git repository; archived version).
- Acknowledgements, including the AI disclosure statement required by IOP.

# What is left out of the compendium

The history of the work and the "pitfalls" narratives (one sentence each where useful); the conservative quadratic model; the angular-acceptance and impact-parameter cuts; the calibration-fit procedure, which survives only as a cross-check in Section 4; everything that refers to simulation data, which goes to paper 2.

# New work needed before writing results

Done on 10 October 2026 (example/worked_example.py, 13 s):

- Rounded tower: radii 160, 350, 550, 820, 1200, 1500 mm; densities after cuts at 100 µm and 100 ps, rounded: 0.039, 0.011, 0.0036, 0.0013, 0.00083, 0.00053 hits/mm^2^ (444–599 hits per tower plane).
- Expected fakes per tower: 6 of 6: 5.2×10^-7^; the six 5-of-6 terms: 7.8×10^-5^ to 1.4×10^-3^; **at least 5 of 6: 3.8×10^-3^**. One fake is reached at 3.0 times these densities.
- Window occupancy η = 0.01–0.03, so distinct fakes ≈ 0.986 × the sum. At the example's cut, c′ ≈ 10^-4^ on the reference detector, so the finite-c corrections are negligible (about 10^-5^) and need no separate Monte Carlo.
- In the process, an error in the earlier calibration (depth-view residuals in r instead of z) was found and corrected in the analysis package; the earlier draft carries a correction note (v1.8).

Still to do:

1. All figures regenerated by scripts in docs/papers/jinst_fake_rate/ (the two validation figures exist; the example needs E vs density factor and vs position resolution).
2. Before quoting them, read Hotelling 1939, Weyl 1939 and Johansen & Johnstone 1990, to cite the tube formula in its standard form.

# Decided

1. Title: the one above.
2. Whitened coordinates from the start; spheres, not ellipsoids.
3. Only the exact helix fit; the B = 0 line as the pedagogical case. No conservative quadratic model.
4. No angular-acceptance or impact-parameter cuts.
5. The end-to-end validation in small, tractable detectors is included.
6. k of N gets its own subsection (2.x); distinct fakes from the window occupancy, checked in the validation.
7. The small-χ² power law is derived in the paper, with references, not assumed known.
8. K is computed from the track-manifold volume; Monte Carlo is used for the corrections and as a check (10 October 2026).
9. The internal 1994 ATLAS notes are not cited.
10. Example operating point: 100 µm position and 100 ps time resolution; futuristic resolutions and R&D-goal scans go to paper 2 (10 October 2026).
