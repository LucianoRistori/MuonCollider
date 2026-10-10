# Formalize your own precursor, then cite three neighbours

The method's individual ingredients are almost all known somewhere. What appears to be new is putting them together. Nothing found in the HEP literature computes expected combinatorial fakes as E = P_true × ∏ n_i with a χ²-based acceptance, a small-cut law whose exponent is fixed by geometry, an exact shrink-and-rescale Monte Carlo calibration, and an algorithm-independent lower-limit interpretation. Every HEP fake-rate study checked in full text (ATLAS FTK, CDF SVT, CMS L1 tracklet, the Ryd–Skinnari review) treats fakes by simulation or qualitatively. The most important prior art is your own paper: Casarsa, Jindariani, Ristori, JINST 20 (2025) P04030, Sec. 8.5. It already extrapolates an inflated-resolution Monte Carlo back to nominal using a "product of resolutions" scaling (7.1^10 for ten degrees of freedom). The new paper should present itself as the rigorous formalization and generalization of that section. Outside HEP, a referee could raise three neighbours, and the paper should cite and contrast them itself:

- **Random-alignment geometry** (Kendall & Kendall 1980; Broadbent 1980; Edmunds & George 1981) already contains the (w/L)^(k−2) law, which is the straight-line case of σ^ndof.
- **The a contrario "number of false alarms"** (Desolneux, Moisan, Morel 2000) is conceptually the same as "expected fakes = number of trials × per-trial probability".
- **Asymptotic Sampling** (Bucher 2009) and **Enhanced Monte Carlo** (Naess et al. 2009) both make a rare event common with an artificial scale parameter and then extrapolate back. Your version differs because the scaling is an exact identity and the exponent is known rather than fitted.

The lower-limit interpretation was the main open risk. The one possible precursor, a 1994 ATLAS internal note by S. Haywood, was read by the authors (October 2026): it studies only ghost tracks from projection ambiguities, not random combinations of noise hits, and it is not publicly accessible, so it is neither a precursor nor citable. On AI assistance, IOP Publishing, whose policy JINST follows, requires a separate AI disclosure statement in the Acknowledgements. It must name the tool and version, say how the tool was used, and say how the output was verified. AI cannot be an author.

## Novelty per element: new as a whole, with three elements at medium risk

The table below gives a verdict for each element of the method. "Risk" means the risk that a referee judges the element not new, or finds that it has been mis-credited.

| Method element | Closest prior art found | Verdict | Risk |
|---|---|---|---|
| E = P_true × ∏ n_i as an exact expectation over all N-plane combinations | Back-of-envelope folklore in trigger design; the accidental-coincidence formula for counters (Leo, ch. 15); the a contrario NFA (Desolneux et al. 2000); the count "product of hits per layer" in Casarsa et al. 2025 Sec. 8.5 | Not new as an idea. Your contribution is to state it exactly, and to note that it needs no independence between overlapping combinations, only linearity of expectation (the same argument Desolneux et al. use) | **Medium**: a referee may call it "well known". Pre-empt this by crediting it explicitly |
| k-of-N extension | Casarsa et al. 2025 drops layers one at a time; radar M-of-N track initiation (Hu et al. 1997, IEEE TSP 45, 445; full text not read) | New in HEP form. Acknowledge the radar analogue | Low |
| Small-cut law P ≈ K·cut^(ndof/2) | Series expansion of the χ² CDF, (x/2)^(ν/2)/Γ(ν/2+1); the generic small-ball argument (density at the origin × ellipsoid volume); the (w/L)^(k−2) alignment law for straight lines | The mathematics is standard, and the exponent for straight lines is known. New: holding the exponent fixed at ndof/2 for any linear model and both views, fitting only K on a plateau, and extending to the helix | **Medium** for the exponent itself; low for its use as an extrapolation anchor |
| Exact shrink-and-rescale calibration, P_true(c) = P_shrunk(c/λ²) | Asymptotic Sampling (Bucher 2009) inflates σ and extrapolates with an asymptotic form; Enhanced MC (Naess et al. 2009) scales the failure domain and fits a four-parameter tail; Casarsa et al. 2025 ran MC at inflated resolution and scaled back | New in its exactness. For linear least squares χ² is homogeneous of degree 2, so the scaling step itself makes no extrapolation error. The only approximation is the small-cut law, which the plateau tests | **Medium**: cite AS and EMC and contrast with them |
| Resolution scaling P ∝ σ^ndof | Casarsa et al. 2025 Sec. 8.5 (heuristic "product of resolutions"); the alignment law (w/L)^(k−2) | Your formalization of your own heuristic. It also settles the ambiguity left in Sec. 8.5 between 15 resolutions and 10 degrees of freedom | Medium-high, but the precursor is your own paper, so the fix is to cite it |
| Helix fit in a solenoid | Circle-fit literature (Karimäki; Chernov–Ososkov; Crawford). No reference found for small-ball behaviour of nonlinear least squares | Apparently new. Rests on the scale-equivariance of the circle model and on the fit being locally regular near the zero-residual manifold. This should be presented as a physics argument checked numerically | Low on novelty; **moderate on rigour** |
| Interpretation as an algorithm-independent lower limit | None found in reviews (Mankel; Strandlie & Frühwirth; Frühwirth & Strandlie). The 1994 ATLAS internal notes turned out to be unrelated (ghosts from projection ambiguities) and are not citable | Appears new | Low to medium; rests on stating precisely which algorithms are bounded |
| Muon Collider–inspired 6-layer barrel example | Muon Collider tracking papers quantify BIB fakes only with specific algorithms and full simulation | Fills a real gap: no algorithm-independent fake floor has been published for the Muon Collider | Low |

Several sources back the HEP side of this table. The FTK associative-memory work describes fake roads only through simulation, for example "the number of fake matched roads increases nearly linearly with the bank size" ([Annovi et al., PoS RD13 014](https://pos.sissa.it/189/014/pdf)). The CDF SVT paper says only that "tracks passing programmable goodness-of-fit cuts propagate downstream" ([arXiv:physics/0306169](https://arxiv.org/pdf/physics/0306169)). The CMS tracklet papers give only a pair-count estimate of about 3600 candidate tracklets per seeding combination, and otherwise measure fakes with an emulator ([arXiv:1706.09225](https://arxiv.org/pdf/1706.09225); [arXiv:1910.09970](https://arxiv.org/pdf/1910.09970)). The HL-LHC track-trigger review treats combinatorial stubs only qualitatively ([Ryd & Skinnari, arXiv:2010.13557](https://arxiv.org/pdf/2010.13557)). INSPIRE full-text searches for "fake track rate" + "occupancy" + "analytic", and for "random combinations of hits" + "chi2", turned up only TDRs, theses and ATLAS performance documents, with no dedicated analytic paper ([INSPIRE query 1](https://inspirehep.net/api/literature?q=fulltext:%22fake%20track%20rate%22%20and%20fulltext:%22occupancy%22%20and%20fulltext:%22analytic%22&size=25&fields=titles,arxiv_eprints,publication_info); [INSPIRE query 2](https://inspirehep.net/api/literature?q=fulltext:%22random%20combinations%20of%20hits%22%20and%20fulltext:%22chi2%22&size=25&fields=titles,arxiv_eprints,publication_info)). This is good evidence but not exhaustive. The FTK TDR, the CDF XFT papers, the CMS Phase-2 Tracker TDR, the Amstutz Hough-transform papers, and the Mu3e, heavy-ion and Mu2e TDRs were not read.

Two technical points came out of this review. They are my own inferences, not findings from the literature, so check them against your derivation. First, in a fixed solenoid field, shrinking the geometry by λ also shrinks the radius of curvature by λ. Any p_T threshold or parameter-range constraint in the fit therefore has to be scaled consistently in the shrunk Monte Carlo, or the identity is no longer exact. Second, the lower-limit claim holds only for a stated class of algorithms. An algorithm that keeps every combination passing the same χ² cut sees all the fakes counted in E. An algorithm that discards some of them, through ambiguity resolution, hit sharing or extra cuts, can produce fewer. A referee will ask for the precise wording. One defensible form is "a lower limit for any algorithm that is fully efficient for genuine tracks passing this χ² cut and applies no information beyond the hits used in the fit". Another is to call E "the irreducible rate at a given χ² acceptance".

## Your 2025 JINST Sec. 8.5 is the direct precursor; the paper should say what it adds

Casarsa, Jindariani and Ristori's Sec. 8.5 ran 10,000 BIB-only events. Zero tracks passed the χ² cut at nominal resolution, against 576 with resolutions inflated to between 100 µm and a few mm. The section took the probability of a "seemingly good track" to be "proportional to the product of the fifteen resolutions of the fifteen coordinates involved". It then extrapolated with an average degradation factor of 7.1 as 7.1^10 = 3.3 × 10^8 "for a fit with five hits and ten degrees of freedom". That gives about 1.7 × 10⁻⁶ fakes per 10⁴ events. The section also stated that the number of combinations "is equal to the product of the number of hits in each layer", with a minimum of 5 layers ([arXiv:2412.14136](https://arxiv.org/pdf/2412.14136); [INSPIRE](https://inspirehep.net/api/literature?q=arxiv:2412.14136&fields=titles,authors.full_name,publication_info,dois)). That is σ^ndof in embryonic form, together with an inflate-and-scale-back calibration. The section has no explicit E = P·∏n_i formula, no small-cut law, and no quoted χ² cut value.

Suggested positioning sentences:

> "In ref. [Casarsa 2025, Sec. 8.5] the fake rate of a specific algorithm was extrapolated from a Monte Carlo sample with inflated resolutions by assuming that the probability of a random combination passing the fit scales as the product of the resolutions. Here we show that this scaling follows exactly from the homogeneity of the least-squares χ², fix its exponent to the number of degrees of freedom, replace the extrapolation in resolution by an exact geometric rescaling, and turn the result into an algorithm-independent lower limit."

> "The power law P ∝ cut^(ndof/2) generalizes the classic result that chance alignments of k random points within a strip of width w scale as (w/L)^(k−2) [Kendall & Kendall 1980; Broadbent 1980; Edmunds & George 1981] to arbitrary linear and helical track models with a χ² acceptance."

> "In the language of a contrario detection [Desolneux et al. 2000], E is the number of false alarms of the track selection; as there, linearity of expectation makes it exact even though overlapping combinations are not independent."

> "Our calibration is related to rare-event techniques that make the event frequent through an auxiliary scale parameter and then extrapolate, such as Asymptotic Sampling [Bucher 2009] and Enhanced Monte Carlo [Naess et al. 2009]. Here, however, the scaling is an exact symmetry of the least-squares χ² and the asymptotic exponent is known in closed form, so only the prefactor K is fitted."

> "Rates of accidental N-fold coincidences in counter systems [Leo] and fake roads in associative-memory triggers [Annovi et al.] follow the same product-of-occupancies logic. To our knowledge no algorithm-independent treatment with a χ²-defined acceptance has been published."

The contrasts with the reliability-engineering methods are concrete. Asymptotic Sampling "progressively inflate[s]" input standard deviations, regresses the scaled reliability index on support points, and extrapolates. It is known to be sensitive to where the support points are placed ([Bayrak & Acar 2021](https://link.springer.com/article/10.1007/s00158-021-03057-0); [Sichani, Nielsen & Bucher 2011](https://vbn.aau.dk/da/publications/applications-of-asymptotic-sampling-on-high-dimensional-structura)). Enhanced Monte Carlo fits p_f(λ) ≈ q·exp{−a(λ−b)^c} with all four parameters free ([secondary source](https://www.academia.edu/57816630/Efficient_System_Reliability_Analysis_by_Finite_Element_Structural_Models)). Your method has no scaling-model error, fixes the exponent, and assumes uniform rather than Gaussian inputs. You can also add a sentence on extreme-value theory: the cut^(ndof/2) law is a regularly varying lower tail. Fitting the plateau of P/cut^(ndof/2) is therefore a peaks-over-threshold estimate with the shape parameter fixed by geometry, which has lower variance than a free-ξ fit (Coles 2001).

One optional analytic cross-check: K = f_r(0) × Vol(unit ellipsoid). Here f_r(0) is the density of the least-squares residual vector at zero. For uniform hits it is set by the volume of the box's intersection with the column space of the design matrix, that is, the fraction of hit space occupied by exactly fitting tracks. This is my own derivation, and no reference was found for it.

## Muon Collider citations and the geometry discrepancy

The coordinator confirmed that the example radii (164, 354, 554, 819, 1153, 1486 mm) come from InnerTracker_o2_v07_01.xml and OuterTracker_o2_v07_01.xml in your simulation, which sit next to MuSIC_v2.xml. The outer three radii match the documented MuColl_v1/MAIA Outer Tracker exactly. The documented Inner Tracker radii are 127, 340 and 554 mm ([MC detector wiki](https://mcd-wiki.web.cern.ch/detector/tracking/)), and MAIA quotes IT 12.7–55.4 cm and OT 81.9–148.6 cm ([MAIA, arXiv:2502.00181](https://arxiv.org/html/2502.00181)). The geometry is therefore most likely that of the newer MUSIC concept ([Andreetto et al., EPJC 86 (2026) 554](https://arxiv.org/abs/2511.23273)). Confirm this with the collaboration. Then cite MUSIC as the source of the example geometry, and cite MAIA and the IMCC overviews for context. Rounding to two significant figures and calling the geometry "inspired by" is appropriate.

For the densities: MAIA quotes about 30,000 BIB hits/cm² per crossing in the innermost pixel layer before timing cuts. No published per-layer IT/OT density after cuts was found, so attribute the example densities to your own analysis of the public BIB samples.

## Must-check before submission

| # | Item | Why | Who / how |
|---|---|---|---|
| 1–2 | (Removed 10 October 2026: the 1994 ATLAS internal notes are not publicly accessible and therefore not citable; the one read does not apply) | | |
| 3 | ATLAS FTK TDR (CERN-LHCC-2013-007), ANIMMA 2011 and NSS/MIC 2012 AM papers | Possible explicit ∏(occupancy × road width) formula | Full text not reached here |
| 4 | CDF XFT papers; CMS Phase-2 Tracker TDR (CERN-LHCC-2017-009); Mu3e, Mu2e, heavy-ion TDRs | Short analytic estimates may be hidden in TDRs | Not read |
| 5 | Kendall & Kendall full text; Edmunds & George | The exact form of the k-point alignment law; the (w/L)^(k−2) form here is reconstructed from a secondary source | Read before quoting |
| 6 | Bucher 2009 regression form (β(f) = A f + B/f) and σ/f vs f·σ notation; Naess margin definition | Needed only if the contrast is stated in formula form | Read the PDFs |
| 7 | Detector concept behind InnerTracker_o2_v07_01 (MUSIC vs MAIA) | The radii disagree with the MuColl_v1 Inner Tracker | Ask the collaboration |
| 8 | Lower-limit wording and the shrink rules for the helix (p_T or curvature range scales with λ) | Referee robustness (own inference) | Authors |
| 9 | Exact wording of the IOP language-editing exemption; whether a 2026 PDG edition exists; whether arXiv:2504.21417 and MAIA are now published | Fetched text was truncated or the status may have changed | Check pages at submission |

## Reference list

Status key: **V** means verified against INSPIRE, Crossref, a publisher page, or an independent reference list during this research. **P** means partly verified, with the unverified fields named. **U** means unverified and must be checked before citing.

### HEP fake rates, track triggers and coincidences

| Citation | Status | Cite for |
|---|---|---|
| M. Casarsa, S. Jindariani, L. Ristori, "A space-time tracking algorithm for high occupancy events at future colliders", JINST 20 (2025) P04030, doi:10.1088/1748-0221/20/04/P04030, arXiv:2412.14136 | V | The direct precursor (Sec. 8.5); the example environment |
| W. Ashmanskas et al., "The CDF Silicon Vertex Trigger", Nucl. Instrum. Meth. A 518 (2004) 532–536, doi:10.1016/j.nima.2003.11.078, arXiv:physics/0306169 | V | χ²-cut track triggering at hadron colliders |
| A. Annovi, S. Amerio, M. Beretta et al., "A new 'Variable Resolution Associative Memory' for High Energy Physics", ANIMMA 2011, doi:10.1109/ANIMMA.2011.6172856 | V (metadata only) | Fake roads vs road width |
| A. Annovi, A. Castegnaro, P. Giannetti et al., "Variable resolution Associative Memory optimization and simulation for the ATLAS FastTracker project", PoS(RD13)014 (2013), doi:10.22323/1.189.0014 | V | Fake roads are treated by simulation only |
| T. Iizawa, "Fast tracker performance using the new 'variable resolution associative memory' for ATLAS", IEEE NSS/MIC 2012, pp. 1392–1395, doi:10.1109/NSSMIC.2012.6551339 | V (metadata only) | Optional AM reference |
| ATLAS Collaboration, "Fast TracKer (FTK) Technical Design Report", CERN-LHCC-2013-007, ATLAS-TDR-021 | P (report number only) | Optional |
| A. Ryd, L. Skinnari, "Tracking Triggers for the HL-LHC", Annu. Rev. Nucl. Part. Sci. 70 (2020) 1–26, arXiv:2010.13557 | V | Combinatorial stubs at HL-LHC, treated qualitatively |
| E. Bartz et al., "FPGA-based tracking for the CMS Level-1 trigger using the tracklet algorithm", arXiv:1910.09970 | P (journal not identified) | Fakes measured by emulation |
| W. R. Leo, *Techniques for Nuclear and Particle Physics Experiments*, 2nd ed., Springer, 1994, doi:10.1007/978-3-642-57920-2, ch. 15 (pp. 303–316) | P (exact page of the accidental-coincidence formula not checked) | The accidental-coincidence analogue |
| Z. Hu et al., "Statistical performance analysis of track initiation techniques", IEEE Trans. Signal Process. 45(2) (1997) 445–456 | P (full author list not seen) | Radar M-of-N analogue of k-of-N |

### Statistics, rare events, a contrario and random alignments

| Citation | Status | Cite for |
|---|---|---|
| D. G. Kendall, W. S. Kendall, "Alignments in two-dimensional random sets of points", Adv. Appl. Prob. 12(2) (1980) 380–424, doi:10.2307/1426603 | V | The chance-alignment ancestor of σ^ndof |
| S. Broadbent, "Simulating the ley hunter", J. R. Stat. Soc. A 143(2) (1980) 109–140, doi:10.2307/2981985 | P (end page 140 not seen) | Same |
| M. G. Edmunds, G. H. George, "Random alignment of quasars", Nature 290 (1981) 481–483, doi:10.1038/290481a0 | P (via secondary source) | The k-point (w/L)^(k−2) law |
| A. Desolneux, L. Moisan, J.-M. Morel, "Meaningful alignments", Int. J. Comput. Vis. 40(1) (2000) 7–23 | V | NFA = number of tests × probability; expectation without independence |
| A. Desolneux, L. Moisan, J.-M. Morel, *From Gestalt Theory to Image Analysis*, Springer, ISBN 978-0-387-74378-3 | P (year 2008 and series not seen) | Optional |
| C. Bucher, "Asymptotic sampling for high-dimensional reliability analysis", Probab. Eng. Mech. 24(4) (2009) 504–510, doi:10.1016/j.probengmech.2009.03.002 | V | The closest MC analogue, to contrast with |
| M. T. Sichani, S. R. K. Nielsen, C. Bucher, "Applications of asymptotic sampling on high dimensional structural dynamic problems", Struct. Saf. 33(4–5) (2011) 305–316, doi:10.1016/j.strusafe.2011.05.002 | V | Asymptotic Sampling's sensitivity to support points |
| A. Naess, B. J. Leira, O. Batsevych, "System reliability analysis by enhanced Monte Carlo simulation", Struct. Saf. 31(5) (2009) 349–355, doi:10.1016/j.strusafe.2009.02.004 | V (via Crossref query; one note guessed a different DOI, so use this one) | Second analogue, with a fitted tail |
| K. Breitung, "Asymptotic approximations for multinormal integrals", J. Eng. Mech. 110(3) (1984) 357–366, doi:10.1061/(ASCE)0733-9399(1984)110:3(357) | P (end page 366 vs 367) | Optional; the opposite (large-β) asymptotic regime |
| S. K. Au, J. L. Beck, "Estimation of small failure probabilities in high dimensions by subset simulation", Probab. Eng. Mech. 16(4) (2001) 263–277 | P (DOI not seen) | Generic alternative rare-event method |
| R. Y. Rubinstein, D. P. Kroese, *Simulation and the Monte Carlo Method*, 3rd ed., Wiley, 2016, ISBN 978-1-118-63216-1 | V (DOI U) | Rare-event MC textbook |
| J. A. Bucklew, *Introduction to Rare Event Simulation*, Springer, 2004, ISBN 978-0-387-20078-1 | P | Optional |
| S. Asmussen, P. W. Glynn, *Stochastic Simulation: Algorithms and Analysis*, Springer, 2007 | P | Optional |
| S. Coles, *An Introduction to Statistical Modeling of Extreme Values*, Springer, London, 2001, doi:10.1007/978-1-4471-3675-0 | V | Peaks-over-threshold with a fixed shape |
| W. V. Li, Q.-M. Shao, "Gaussian processes: inequalities, small ball probabilities and applications", Handbook of Statistics 19 (2001) 533–597, doi:10.1016/S0169-7161(01)19019-X | V | Small-ball terminology (Gaussian and infinite-dimensional, so background only) |
| χ² small-x law: DLMF §8.7 or Abramowitz & Stegun 6.5.29 | U | The incomplete-gamma series |
| H. Hotelling, "Tubes and spheres in n-spaces, and a class of statistical problems", Amer. J. Math. 61(2) (1939) 440–460 | P (pages from memory) | The volume of a tube around a curve: the geometric K (added 10 Oct 2026) |
| H. Weyl, "On the volume of tubes", Amer. J. Math. 61(2) (1939) 461–472 | P (pages from memory) | The general tube formula; leading term = area of the manifold × sphere volume |
| S. Johansen, I. M. Johnstone, "Hotelling's theorem on the volume of tubes: some illustrations in simultaneous inference and data analysis", Ann. Statist. 18(2) (1990) 652–684 | V | Statistical use of the tube formula for nonlinear least squares |
| S. Kuriki, A. Takemura, "The volume-of-tubes formula: computational methods and statistical applications", arXiv:math/0511502 | V | Review |

### Muon Collider

| Citation | Status | Cite for |
|---|---|---|
| P. Andreetto et al., "MUSIC: a multi-purpose detector concept for physics at the 10 TeV muon collider", Eur. Phys. J. C 86 (2026) 554, doi:10.1140/epjc/s10052-026-15654-8, arXiv:2511.23273 | V | The likely source of the example geometry (confirm) |
| C. Bell et al., "MAIA: A new detector concept for a 10 TeV muon collider", arXiv:2502.00181 | V (no journal yet) | Detector context; BIB density |
| C. Accettura et al., "Towards a muon collider", Eur. Phys. J. C 83 (2023) 864, doi:10.1140/epjc/s10052-023-11889-x; Erratum ibid. 84 (2024) 36 | V | General reference |
| C. Accettura et al. (IMCC), "Interim report for the International Muon Collider Collaboration", CERN Yellow Rep. Monogr. 2/2024, doi:10.23731/CYRM-2024-002, arXiv:2407.12450 | V | Baseline design and BIB |
| C. Accettura et al. (IMCC), "The Muon Collider", arXiv:2504.21417 | V (unpublished as checked) | ESPP input (optional) |
| K. M. Black et al., "Muon Collider Forum report", JINST 19 (2024) T02015, doi:10.1088/1748-0221/19/02/T02015 | V | Snowmass context |
| M. Casarsa, D. Lucchesi, L. Sestini, "Experimentation at a muon collider", Ann. Rev. Nucl. Part. Sci. 74 (2024) 233–261, doi:10.1146/annurev-nucl-102622-011319 | V | BIB review |
| N. Bartosik et al., "Detector and Physics Performance at a Muon Collider", JINST 15 (2020) P05001, doi:10.1088/1748-0221/15/05/P05001 | V | First full simulation with BIB |
| F. Collamati et al., "Advanced assessment of beam-induced background at a muon collider", JINST 16 (2021) P11009, doi:10.1088/1748-0221/16/11/P11009 | V | FLUKA BIB |
| V. Di Benedetto et al., "A Study of Muon Collider Background Rejection Criteria in Silicon Vertex and Tracker Detectors", JINST 13 (2018) P09004, doi:10.1088/1748-0221/13/09/P09004 | P (journal fields inferred from the DOI) | Timing and double-layer cuts |

### Track fitting and reconstruction

| Citation | Status | Cite for |
|---|---|---|
| R. Frühwirth, "Application of Kalman filtering to track and vertex fitting", Nucl. Instrum. Meth. A 262 (1987) 444–450, doi:10.1016/0168-9002(87)90887-4 | V | The Kalman fit |
| P. Billoir, "Track fitting with multiple scattering: a new method", Nucl. Instrum. Meth. 225 (1984) 352–366, doi:10.1016/0167-5087(84)90274-6 | V | Write it as "Nucl. Instrum. Meth." (before the A/B split) |
| P. Billoir, "Progressive track recognition with a Kalman-like fitting procedure", Comput. Phys. Commun. 57 (1989) 390–394, doi:10.1016/0010-4655(89)90249-X | V | Combinatorial track following |
| P. Billoir, S. Qian, "Simultaneous pattern recognition and track fitting by the Kalman filtering method", Nucl. Instrum. Meth. A 294 (1990) 219–228, doi:10.1016/0168-9002(90)91835-Y | V | Same |
| V. Karimäki, "Effective circle fitting for particle trajectories", Nucl. Instrum. Meth. A 305 (1991) 187–191, doi:10.1016/0168-9002(91)90533-V | V | Helix/circle fit |
| N. I. Chernov, G. A. Ososkov, "Effective algorithms of circle fitting", Comput. Phys. Commun. 33 (1984) 329–333, doi:10.1016/0010-4655(84)90137-1 | V (title uses "of") | Circle fit |
| J. F. Crawford, "A non-iterative method for fitting circular arcs to measured points", Nucl. Instrum. Meth. 211 (1983) 223–225, doi:10.1016/0167-5087(83)90575-6 | V | Circle fit |
| K. Levenberg, Q. Appl. Math. 2 (1944) 164–168, doi:10.1090/qam/10666; D. W. Marquardt, J. Soc. Ind. Appl. Math. 11 (1963) 431–441, doi:10.1137/0111030 | V | Nonlinear least squares (LM) |
| R. L. Gluckstern, "Uncertainties in track momentum and direction, due to multiple scattering and measurement errors", Nucl. Instrum. Meth. 24 (1963) 381–389, doi:10.1016/0029-554X(63)90347-1 | V | Resolution formulae |
| CMS Collaboration, "Description and performance of track and primary-vertex reconstruction with the CMS tracker", JINST 9 (2014) P10009, doi:10.1088/1748-0221/9/10/P10009 | V | Fake-rate definitions at the LHC |
| ATLAS Collaboration, "Performance of the ATLAS track reconstruction algorithms in dense environments in LHC Run 2", Eur. Phys. J. C 77 (2017) 673, doi:10.1140/epjc/s10052-017-5225-7 | V | Same |
| X. Ai et al., "A Common Tracking Software Project", Comput. Softw. Big Sci. 6 (2022) 8, doi:10.1007/s41781-021-00078-8 | V | ACTS (the Muon Collider reconstruction) |
| P. V. C. Hough, US Patent 3,069,654 (1962); R. O. Duda, P. E. Hart, Commun. ACM 15 (1972) 11–15, doi:10.1145/361237.361242 | V | Hough-like methods (background for Casarsa et al.) |
| S. Amrouche et al., "The Tracking Machine Learning Challenge: Throughput Phase", Comput. Softw. Big Sci. 7 (2023) 1, doi:10.1007/s41781-023-00094-w | V | Optional; fake metric definitions |

### Reviews and textbooks

| Citation | Status | Cite for |
|---|---|---|
| R. Mankel, "Pattern recognition and event reconstruction in particle physics experiments", Rep. Prog. Phys. 67 (2004) 553, doi:10.1088/0034-4885/67/4/R03, arXiv:physics/0402039 | V (end page 622 U) | Combinatorial explosion |
| A. Strandlie, R. Frühwirth, "Track and vertex reconstruction: From classical to adaptive methods", Rev. Mod. Phys. 82 (2010) 1419–1458, doi:10.1103/RevModPhys.82.1419 | V | Standard review |
| R. Frühwirth, A. Strandlie, *Pattern Recognition, Tracking and Vertex Reconstruction in Particle Detectors*, Springer Cham, 2021 (open access), doi:10.1007/978-3-030-65771-0 | V (© 2021, not 2020) | Standard textbook |
| R. Frühwirth, M. Regler, R. K. Bock, H. Grote, D. Notz, *Data Analysis Techniques for High-Energy Physics*, 2nd ed., Cambridge University Press, 2000, ISBN 9780521635486 | P (**do not use doi:10.1017/CBO9780511534836, which belongs to another book**; 10.1017/CBO9781139142199 is U) | χ² and least squares |
| S. Navas et al. (PDG), "Review of Particle Physics", Phys. Rev. D 110 (2024) 030001, doi:10.1103/PhysRevD.110.030001 | V (chapter authors U; check for a 2026 edition) | Statistics chapter (χ² conventions) |

## IOP and JINST require a named-tool disclosure in the Acknowledgements

JINST has no AI policy of its own. It "applies the principles of publication ethics endorsed by IoPP" ([JINST about](https://jinst.sissa.it/jinst/help/helpLoader.jsp?pgType=about)). IOP Publishing's page "Generative AI Tools" (metadata last modified 2026-09-28) requires authors who used generative AI to "disclose this in a separate AI disclosure statement as part of the Acknowledgement section". Where the use is relevant to methodology, it must "additionally [be] described in detail in the Methods". The statement must "identify the tool or model used, describe how it was used" and "explain the extent to which the authors verified, validated and critically reviewed the AI-assisted outputs". Routine spelling, grammar and language editing does not need a separate statement. AI tools "cannot be legally accountable for published work", and "all named authors are fully responsible for all material presented in their manuscript", so AI cannot be listed as an author ([IOP Publishing Support](https://publishingsupport.iopscience.iop.org/questions/generative-ai-tools/)). The page also covers two related points. AI-generated explanatory figures must be clearly distinguishable from data. In responding to referees, AI may be used only to polish language, and reviewer reports must not be uploaded to AI tools.

IOP's template is: *"During the preparation of this manuscript, the authors used [Tool Name, Version] to support [describe use]."* For this paper, a statement could read:

> "During the preparation of this manuscript, the authors used Claude (Anthropic, [model/version]) to support the literature search and the verification of bibliographic references, and [drafting/editing of text, and/or development of analysis code — state what applies]. All references were checked against the original sources, and all derivations, code and numerical results were independently verified by the authors, who take full responsibility for the content."

If AI helped write the Monte Carlo or fitting code, also describe that briefly in the Methods section.

## Conclusion

The defensible claim is not "a new formula for fakes". It is that a heuristic that track-trigger designers, including the authors, have used informally for decades can be made exact and calibrated, at the 10⁻²⁶ level, with crude Monte Carlo. That framing turns each likely referee objection into a citation the paper already makes: alignments for the exponent, the a contrario NFA for the expectation, Asymptotic Sampling and Enhanced Monte Carlo for the calibration, and your own Sec. 8.5 for the scaling. The exact λ² identity and the fixed exponent are what lift the method above its relatives.

The interpretation as an algorithm-independent lower limit is the paper's most quotable result and its least protected one. It depends on a careful statement of which algorithms it bounds, which should be settled before submission.


## Update, 10 October 2026

The track-manifold form of K (paper 1, draft 4 of the outline) has a direct mathematical ancestor in the volume-of-tubes formula of Hotelling and Weyl (1939), widely used in statistics for significance tests on nonlinear regression models (Johansen & Johnstone 1990; Kuriki & Takemura). No use of it to compute combinatorial fake rates in tracking was found. It should be cited as the source of the leading-order result; what is new is its application to random hit combinations, with a uniform rather than Gaussian measure.
