# BIB n-tuple format and first plots

## How to run the analysis

Everything is run from the Analysis folder in Dropbox, which holds the
three editable settings files:

```
cd ~/Dropbox/Documents/MuonColliderSimulation/Analysis
# edit __cuts_config.txt (cuts), __smearing_config.txt (resolutions)
# and/or __input_files_config.txt (which data and geometry files to use)
./run_all
```

Input files live in two folders of `MuonColliderSimulation`: the ROOT
files in `Data/`, the detector-geometry XML files in `Geometry/`.
`__input_files_config.txt` just names them (BIB plus, minus, ipp, the
combined BIB file, the signal sample, and the 4 geometry files). The
combined BIB file is derived from plus + minus (+ ipp): before running,
`./run_all` checks that it contains exactly those files (same entries,
run/event numbers and hit counts, in order) and builds or rebuilds it
automatically when it doesn't - e.g. after new plus/minus/ipp files
were put in `Data/`.

`./run_all` first checks the settings files (misspelled, unknown or
missing settings, values that aren't numbers or true/false, negative
widths/sigmas, duplicates - nothing runs until they are fixed), the
input ROOT files and the Python packages. It then runs steps 1-4 and
archives the whole run in `runs/<date>_<time>/`: every plot and table,
the two settings files exactly as used (copied at the start, so editing
them during a run has no effect on it), the complete log
(`run_log.txt`), the code version (`code_version.txt`) and a short
`summary.txt` (settings + BIB rejection / signal efficiency per
subsystem). Only when every step succeeds are the `step*` folders next
to the settings files replaced by that run's results; `latest_run.txt`
says which run they show. A failed or interrupted run is kept as
`runs/<date>_<time>_FAILED` (or `_INTERRUPTED`) and leaves the `step*`
folders untouched. To compare runs: `cat runs/*/summary.txt`.

The 14 key plots of each run are also gathered in one folder,
`_highlights/` (in the Analysis folder for the latest run, and in each
`runs/<date>_<time>/`): the BIB density before/after cuts, the
track-finding efficiency vs pT, and the z-intercept, pT and corrected-
time distributions with the signal overlaid - both without cuts (step
3, `*_with_signal.png`) and N-1 (step 4, `*_n1.png`), full range and
zoomed. The list is `HIGHLIGHTS` at the top of `run_all_steps.sh`.

The code folder (`~/code/MuonCollider`) keeps only templates of the
settings files (`templates/`); the copies in the Analysis folder are the
ones used, and each run keeps a copy of all three. `./run_all` finds the
code in `~/code/MuonCollider` (or set `MUONCOLLIDER_CODE=/path`), reads
the input files from `Data/` and `Geometry/` in the folder above the
Analysis folder (or set `SIM_DIR=/path`), and calls `run_all_steps.sh`
(see `./run_all --help`). Each run's `code_version.txt` also records the
size and date of every input file used.

## Source files
`ntu_bib_plus_1evt.root` and `ntu_bib_minus_1evt.root`, tree `HTAtree` —
1 entry each, one seed/primary particle (`part_pdg = -14`, muon
antineutrino) whose secondaries deposited ~8.05M hits in the tracker
(8,050,751 for "plus", 8,057,637 for "minus" — likely charge-conjugate
BIB samples; their incidence-angle and z-axis-intercept signs are
consistently flipped between the two files, another confirmation
they're mirror samples).

**Merged file**: `ntu_bib_2evt.root` (in the same folder) combines both
into a single 2-entry file — one entry per original file, all branches
carried over unchanged, nothing recomputed — via `merge_bib_files.py`:
```
python3 merge_bib_files.py <plus.root> <minus.root> <output.root>
```
`bib_common.load_hits()` now flattens hits across **all** entries of a
tree by default (`entry=None`), so it reads this combined file
transparently (16,108,388 total hits = exact sum of plus+minus) and
still gives identical results on the original single-entry plus/minus
files (no other code needed to change). Going forward, point any script
at `ntu_bib_2evt.root` to get combined plus+minus results in one run
instead of running twice and adding by hand. The original plus/minus
files and their separate results are kept too (still useful as a
charge-conjugate cross-check, per the "Notable feature" note above).

The merged file's total hit content (plus+minus summed) is taken to
represent the BIB background expected in a **single MC collision**
(bunch crossing) — this matters for how signal is compared to it below
(see "Signal sample and overlay").

## Hit branches (per-hit parallel std::vectors)
`hit_index`, `hit_mcp` (unused, historical), `hit_id0` (cellID, see below),
position `hit_x/y/z` (mm), local coords `hit_u/v` (mm) with resolution
`hit_du/dv` (mm), local axes `hit_ux/uy/uz`, `hit_vx/vy/vz`, momentum
`hit_px/py/pz` (GeV), time `hit_t` (ns).

## hit_id0 cellID encoding
`"system:5,side:-2,layer:6,module:11,sensor:8"` (32-bit uint), decoded as:
```
system = id0 & 0x1F
side   = (id0 >> 5) & 0x3
layer  = (id0 >> 7) & 0x3F
module = (id0 >> 13) & 0x7FF
sensor = (id0 >> 24) & 0xFF
```
- system: 1=VXD barrel, 2=VXD endcap, 3=IT barrel, 4=IT endcap,
  5=OT barrel, 6=OT endcap
- side: 0 = barrel; 1 / 3 = endcap +z / -z
- layer: 0-based layer/disk index within the subsystem

"Detector region" for this analysis = a distinct (system, side, layer)
combination — the maximum granularity available. 41 such regions exist
in this file (5 VXD barrel layers, 4+4 VXD endcap disks, 3 IT barrel
layers, 7+7 IT endcap disks, 3 OT barrel layers, 4+4 OT endcap disks).
Full per-region hit counts and (r, z) extents: see `region_summary.csv`
(generated by `make_basic_plots.py`).

## Subsystem hit fractions (this file)
VXD barrel 31%, IT barrel 25%, OT barrel 22%, OT endcap 8%,
IT endcap 8%, VXD endcap 7%.

## Plots
All density plots use a **log y-axis** — hit density spans 5+ orders
of magnitude across subsystems (VXD barrel ~9 hits/mm² down to OT
endcap ~0.01 hits/mm²), so a linear scale hides all but the densest
bars.

The endcap disk plot (`density_per_disk_endcaps.png`) orders disks by
their **signed** z (mean_z, ascending) — a single continuous axis from
the most negative -z disk through to the most positive +z disk — not
by |z| grouped by side. Same convention should be used for any future
per-disk plot.

## Hit density: mean and peak
Density is reported in hits/mm² rather than raw hit counts, and two
different metrics are computed per region:

- **Mean density** = total hits in the region / true sensitive area of
  the region, from the detector geometry (see below).
- **Peak density** = a robust 2D "hot spot" estimate: hit positions are
  unrolled to a flat local 2D coordinate (barrel: arc-length = mean_r
  × phi, and z; endcap: x, y directly), binned into a 2D grid, and the
  99th percentile of the nonzero bin densities is reported (the raw
  max is also kept for reference, as `peak_density_max_hits_per_mm2`).

**Bin size is chosen per (system,side,layer) region, not fixed
globally** (`bib_common.add_peak_density`, called with
`bin_size_mm=None`): mean hit density varies by 5+ orders of magnitude
across regions, so one fixed bin size is either too fine for sparse
regions or too coarse for dense ones. The bin size is derived from
each region's own mean density to target ~20 hits per occupied bin
(`PEAK_TARGET_HITS_PER_BIN` in `make_basic_plots.py`), clipped to
[1mm, 200mm]. Two diagnostic fields are recorded per region in
`region_summary.csv`: `peak_bin_size_mm` (the bin size actually used)
and `peak_mean_hits_per_bin` (the achieved statistics) — check these
if a region's peak value still looks noisy.

Peak density is always ≥ mean density everywhere; that consistency
check caught two real bugs in turn, both from before adaptive binning
was in place:
1. An early mean-density area estimate (inferred from which
   hits/modules were struck) badly underestimated the true installed
   area for the sparser IT/OT barrel layers — fixed by switching to
   true geometry-based area (see next section).
2. With a single fixed 1mm × 1mm bin, sparse regions (IT/OT endcap
   especially) had close to 1 hit per occupied bin — with only one
   simulated event, the "peak" bin count was then just small integers
   (1, 2, 3 hits/mm²) dominated by single-hit Poisson noise, not real
   spatial structure, so many disks/layers showed suspiciously
   identical peak values. A uniform 10mm bin helped densely-populated
   regions but still left the sparsest OT endcap disks at only
   ~1.3 hits/occupied-bin. Fixed by moving to the adaptive per-region
   bin size described above, which reaches ~20-26 hits/bin everywhere,
   including OT endcap.

## Sensitive area: from the true detector geometry, not hits
Region area is computed from the dd4hep/lcgeo compact XML geometry
(key4hep/k4geo, `MuColl/MuSIC/compact/MuSIC_v2/MuSIC_v2.xml` and its
includes `Vertex_o2_v06_01.xml`, `InnerTracker_o2_v07_01.xml`,
`OuterTracker_o2_v07_01.xml`), not inferred from which modules were
hit in a given event — the latter undercounts area whenever an event
doesn't strike every installed module (a big effect for the sparser
IT/OT barrel layers: mean density there dropped by roughly 40x after
switching to the geometry-based area). All 4 XML files live in
`~/Dropbox/Documents/MuonColliderSimulation/Geometry/` (named in
`__input_files_config.txt`) and are parsed by `geometry.py` (`build_area_lookup`), which resolves the XML's
symbolic constants and computes exact area for VXD barrel/endcap and
IT/OT barrel/endcap (41/41 regions matched). Local copies of the
geometry files are used in preference to fetching from GitHub.

## Track incidence angle (longitudinal / transverse), relative to sensor normal
For each hit, the track momentum direction (`hit_px/py/pz`) is compared
to the local sensor-plane normal to characterize how obliquely the
particle struck the detector element. Implemented as
`bib_common.add_incidence_angles(hits)`, diagnosed/plotted by
`incidence_angle_plots.py`. These diagnostic plots are treated as a
standard part of the results (not just a one-off check), so they're
generated and archived alongside the step 1 density plots on every run.

**Geometry**, per hit at position (x,y,z), r = sqrt(x²+y²):
- Local cylindrical frame: ρ̂ = (x,y,0)/r (radially outward),
  φ̂ = (-y,x,0)/r (azimuthal), ẑ = (0,0,1).
- Sensor normal n̂ = normalize(**u** × **v**), from the module's own
  local axes (`hit_ux/uy/uz`, `hit_vx/vy/vz` — already per-hit in the
  ntuple, so this captures the *true* module orientation, including
  any tilt, not an assumed barrel/endcap geometry), oriented outward
  (n̂ · hit_position > 0).
- Track direction p̂ = (px,py,pz)/|p|.

**Longitudinal component** (`theta_long_deg`) — signed angle between p̂
and n̂ as seen *within the meridian plane* (spanned by ρ̂ and ẑ, i.e.
the plane containing the Z axis and the hit — this is the plane the
user asked for):
```
theta_long = atan2(p_rho*n_z - p_z*n_rho, p_rho*n_rho + p_z*n_z)
```
**Transverse component** (`theta_trans_deg`) — signed angle
representing p̂'s azimuthal (X-Y plane) tilt relative to the normal:
```
theta_trans = atan2(p . phi_hat, p_hat . n_hat)
```
Also stored: `theta_full_deg` = arccos(p̂·n̂) (the plain 3D incidence
angle, a cross-check reference) and `n_tilt_deg` = arcsin(|n̂·φ̂|) (how
far the module's own normal sits out of the meridian plane).

**Verified on data**: endcap disks are exactly flat (n_tilt_deg = 0 to
machine precision), so there theta_long/theta_trans reconstruct
theta_full essentially exactly (tan² identity holds to ~1e-7 relative).
Barrel layers have a genuine tilt (real detector feature, likely a
ladder overlap/stereo angle): VXD barrel n_tilt_deg mean 6.1° (up to
16°), IT barrel 1.6°, OT barrel 0.4° — decreasing with radius. This
means for VXD barrel especially, theta_long/theta_trans don't exactly
recombine to theta_full (~11% median relative difference in the tan²
identity) — expected and not a bug, just a reminder that the two
components are each well-defined and meaningful per hit, but the
"add back up to the full angle" check is only approximate where the
module is tilted.

**Notable feature in this sample**: theta_long is *bimodal* per
subsystem (a dip near 0°, symmetric peaks around ±30–60°) — hits
landing at exactly normal incidence are disfavored relative to oblique
hits. theta_trans, by contrast, peaks sharply at 0° with symmetric
falling tails, as expected. Plausible physical explanation: BIB
particles are dominated by secondaries streaming roughly along the
beam direction rather than radially outward from the IP, so barrel
layers get hit obliquely more often than head-on. Not yet confirmed
against expectations by the user.

## Z-axis intercept of the meridian-plane track (companion to theta_long)
For each hit, `z_axis_intercept_mm` (also in `bib_common.add_incidence_angles`)
extrapolates the hit's local (ρ, z) trajectory — a straight line through
the hit position with slope (p_rho, p_z), i.e. exactly the meridian-plane
projection that theta_long is computed from — back to where it crosses
the Z axis (ρ=0):
```
z_axis_intercept = z - r * (p_z / p_rho)
```
A purely kinematic quantity (no dependence on the sensor normal), = NaN
for tracks running along z within the meridian plane (p_rho ≈ 0, an
undefined intercept — left as NaN rather than ±infinity). Plotted per
subsystem (`z_axis_intercept_per_subsystem.png`), clipped to ±3000mm
(`Z0_RANGE_MM` in `incidence_angle_plots.py`); hits outside that range
or with undefined p_rho are excluded from the histogram, not binned at
the edges — the excluded fraction is reported per subsystem (typically
<1%, but up to ~5-6% for IT/OT endcap, which have more near-grazing
hits). **Result**: sharply peaked near z≈0 with realistic tails in both
directions — consistent with most BIB secondaries' local trajectories
pointing back toward the IP region.

**Zoomed version** (`z_axis_intercept_per_subsystem_zoom.png`): same
quantity, same 121 bins, but restricted to ±100mm (`Z0_ZOOM_RANGE_MM`)
instead of ±3000mm — i.e. much finer resolution (~1.65mm/bin vs.
~50mm/bin) right around z=0, where the full-range plot only has a
couple of bins to work with. Excluded fraction (outside ±100mm, or
undefined) is printed per subsystem — much larger than the full-range
exclusion, since this window is much narrower (roughly 51% for VXD
barrel up to ~89% for IT/OT endcap). **Result, unexpected**: rather
than being peaked at z=0 as the coarse full-range view suggested, the
zoomed view shows a *dip* right at z=0, rising toward both edges of
the ±100mm window, in every subsystem — a real feature the full-range
binning was too coarse to resolve, not yet explained. (The signal
overlay below shows the signal sample fills in exactly this dip — see
"Signal sample and overlay" section.)

## Transverse curvature 1/R (companion to theta_trans)
For each hit, `inv_radius_per_mm` (also in `bib_common.add_incidence_angles`)
is the curvature of the circle in the transverse (X-Y) plane that passes
through both the origin (beam axis) and the hit, and is tangent there to
the track's transverse-momentum direction — i.e. the curvature a particle
produced at the origin would need in order to reach this hit going in this
transverse direction. Derived from the standard tangent-chord relation for
a circle: the chord from the origin to the hit has length r and subtends
angle 2ψ at the center, where ψ is the angle between the transverse
momentum direction and ρ̂ (the chord direction), giving chord = 2R·sin(ψ):
```
psi = atan2(p_phi, p_rho)          # also stored as psi_transverse_deg
inv_radius_per_mm = 2 * sin(psi) / r
```
Unlike `z_axis_intercept_mm`, this is bounded for every hit (every hit has
r > 0, so no divide-by-zero/undefined case), with |1/R| ≤ 2/r — no
exclusion or clipping needed. **Verified on data**: well-behaved, peaked
near 0 with realistic tails, spread scaling down with radius as expected
(VXD barrel p01/p99 ≈ ∓64/m, down to OT barrel/endcap ≈ ∓1-2/m); plus/minus
files mirror each other in sign, consistent with the charge-conjugate
sample check.

**Plot axis relabeled in GeV/c** (`inv_radius_per_subsystem.png`): the
histogram itself is unchanged (same bins, same ±80 1/m range, still
literally 1/R in 1/m internally), but the x-axis ticks are now printed as
transverse momentum `pT [GeV/c] = 0.3 * B[T] * R[m]`, assuming B = 5T
(`B_FIELD_T` in `incidence_angle_plots.py`). Because pT is a *reciprocal*
of 1/R, this is a nonlinear relabeling, not a rescaling: tick **positions**
stay evenly spaced in 1/R (`pt_ticks_for_axis()`, step = 20 1/m) — so they
don't crowd together — and only the printed numbers become nonlinear pT
values, converging toward ±∞ at the center (1/R = 0, a perfectly straight,
infinite-pT track) and shrinking toward ~0.019 GeV/c at the ±80 1/m edges.
Sign of pT follows the sign of the curvature. B is only a labeling
assumption here — the stored `inv_radius_per_mm`/`psi_transverse_deg`
values themselves don't depend on it, so relabeling for a different field
strength just changes `B_FIELD_T`, no re-run of the underlying calculation
needed.

**Zoomed version** (`inv_radius_per_subsystem_zoom.png`): crops the same
curvature axis as the full-range plot above (still centered and binned in
1/R, **not** re-binned directly in pT) down to the narrow window
`|1/R| <= GEV_PER_INV_M / PT_ZOOM_RANGE_GEV` — i.e. exactly where the
relabeled pT tick would read ±`PT_ZOOM_RANGE_GEV` (±1 GeV/c) at the
edges — then rebins that window finely (121 bins) and relabels ticks with
`pt_ticks_for_axis()`, same as the full plot. This preserves the same
center-high/edge-low pT sense as the full-range plot (important: an
earlier version of this plot instead re-binned directly in pT, which put
pT=0 at the center and made pT increase outward — backwards relative to
the full-range plot; that version has been replaced). It reveals the
internal structure of what was a single central spike in the full-range
plot. Excluded fraction (`|1/R|` outside this window, i.e. the
sharply-curved, low-pT tail below `PT_ZOOM_RANGE_GEV`) is printed per
subsystem: ~89% for VXD barrel down to <1% for OT barrel — confirming VXD
is dominated by soft/curved particles and OT by hard/straight ones.

Diagnostic outputs: per-subsystem histograms of theta_long, theta_trans,
theta_full, z_axis_intercept_mm (full range + zoomed to ±100mm),
inv_radius_per_mm (full range relabeled in pT GeV/c, plus a zoomed,
same-convention version cropped to the |pT| >= 1 GeV/c window), plus a 2D
(theta_long, theta_trans) histogram, and `incidence_angle_summary.csv`
(mean/std/percentiles per subsystem, plus z_axis_intercept_mm's finite
fraction/percentiles and inv_radius_per_m's mean/std/percentiles — the
CSV itself still stores 1/R in 1/m, not pT, and doesn't yet carry the
zoomed-window exclusion fractions, which are console-only for now).
Run on the original plus/minus files (`step2_incidence_angles/`,
`step2_incidence_angles_minus/`) and on the merged file
(`step2_incidence_angles_combined/`). Not yet used to define any density
cuts.

## Time-of-flight correction (companion diagnostic to step 2)
For each hit, `bib_common.add_time_of_flight(hits)` computes the time a
muon would take to travel from the IP (origin) to the hit position along
its actual curved trajectory, and subtracts that expected time-of-flight
from the measured hit time `t`, giving `t_corrected_ns`. Requires
`add_incidence_angles(hits)` to have been called first (reuses
`psi_transverse_deg`/`inv_radius_per_mm`). For a genuine on-time signal
muon from the IP, `t_corrected_ns` should cluster tightly around zero
regardless of how far out the hit is — verified on data, see below —
unlike the raw hit time, which simply grows with distance travelled.

**Geometry**, two separate pieces:
1. **Transverse arc length** from the origin to the hit, along the same
   circle already used for `inv_radius_per_mm` (passes through the
   origin and the hit, tangent there to the transverse momentum
   direction, subtending angle 2ψ at its center):
   ```
   s_transverse = r * psi / sin(psi)     # psi in radians
   ```
   Numerically stable as psi→0 (straight track), where it correctly
   limits to s_transverse → r.
2. **Total (3D) path length**: a helix in a solenoidal (z-directed)
   field has constant pitch (pz and transverse-momentum magnitude are
   each separately conserved along the true path), so
   `s_total = s_transverse / sin(theta_polar)`, where
   `sin(theta_polar) = pT_raw/p_raw` uses the hit's own raw local
   momentum direction (the real, known pitch angle).

**Momentum/speed estimate**: per the request, the transverse momentum
magnitude used for the speed (beta) calculation comes from the
curvature (`inv_radius_per_mm` → pT = 0.3·B·R, same as the pT axis on
`inv_radius_per_subsystem.png`), not directly from the raw hit
momentum — this is what a real curvature-based track-pT measurement
would give — combined with the same raw-direction polar angle above:
```
pT_curv = 0.3 * B_FIELD_T / |inv_radius_per_m|
p_estimate = pT_curv / sin(theta_polar)
beta = p_estimate / sqrt(p_estimate^2 + m_muon^2)
tof_expected_ns = s_total / (beta * c)          # c = 299.792458 mm/ns
t_corrected_ns = hit_t - tof_expected_ns
```
A fixed particle mass (muon, PDG value 0.1056583755 GeV, by default) is
assumed for every hit regardless of what actually produced it — that's
the point: this tests which hits are consistent with an on-time muon
from the IP, whatever they actually are.

**Verified on data** (`time_of_flight_plots.py`, run per-subsystem on
both the signal file and the merged BIB file):
- **Signal**: `t_corrected_ns` peaks sharply at exactly 0 in every
  subsystem (median ≈ 0.000-0.002 ns everywhere; p05/p95 within about
  ±0.1-0.6 ns in the inner subsystems), confirming the correction works
  as intended.
- **BIB**: does *not* peak at zero — broad, offset distribution, median
  +0.3 to +1.1 ns depending on subsystem, no sharp central spike. So
  `t_corrected_ns` is a genuine BIB-vs-signal discriminant in its own
  right, on top of the angle/curvature variables from step 2.
- **Outlier tail (both samples, expected and diagnosed, not a bug)**: a
  minority of hits (~5.5% for signal, ~27% for BIB) show very large
  |t_corrected_ns| (up to hundreds of ns). These concentrate almost
  entirely among hits with very low raw p_T (median ~2 MeV, vs. ~2.7 GeV
  for the bulk in the signal sample) — i.e. low-momentum secondary/
  knock-on hits, not the primary track, for which the "originated at
  the IP" assumption behind the whole correction is simply wrong. Large
  excursions there are an expected consequence of that assumption
  breaking down for non-primary hits, not a defect in the calculation.
- **Asymmetric shape (discussed with the user)**: the corrected-time
  distribution isn't just shifted positive on average, it's genuinely
  skewed - suppressed just below 0, then a shoulder/bump around
  +0.3-0.5 ns, then a long positive tail, in most subsystems. Likely
  explanation: t_corrected_ns is a *minimum-delay* construction (the
  fastest a muon born at the IP at t=0 could reach the hit along the
  idealized helix), so nothing can arrive earlier than that floor, but
  off-IP origin, extra scattering, or a longer real path than the
  idealized helix can only push a hit later, never earlier - a
  one-sided physical constraint that naturally produces a
  suppressed-negative, long-positive-tail shape rather than a
  symmetric one. Not yet decomposed into sub-populations (e.g. the
  dip-to-shoulder feature might separate hits that scatter once near
  the IP region from hits originating further upstream along the
  beamline) - a possible follow-up.

Outputs, in `Analysis/step2_time_of_flight[_minus|_combined|_signal]/`:
`hit_time_raw_per_subsystem.png` (measured t, unchanged from before),
`time_corrected_per_subsystem.png` (full ±20ns range),
`time_corrected_per_subsystem_zoom.png` (±2ns zoom, shows the peak
clearly), and `time_of_flight_summary.csv` (n_hits, raw/corrected
mean/std/percentiles, and the >1ns outlier fraction, per subsystem).

## Signal sample and overlay (step 3)
`ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root` (same folder,
same `HTAtree` format/branches as the BIB files): a muon-gun physics
signal sample, 100,000 independent single-muon events (`part_pdg` = ±13),
theta in [10°,170°], phi in [0°,360°], and a ±1.5mm vertex-z spread
(`dz1p5`) — 1,331,666 hits total, ~13.3 hits/event. **The filename is
misleading**: despite "pt1p5GeV", the sample is not generated at a
fixed pT. Verified on the generated `part_px`/`part_py`: pT is flat in
**1/pT**, with a floor at pT = 1.5 GeV/c (`pt.min() = 1.5000042`) and a
long tail out to very high pT (`pt.max() = 120341.9`, mean 17.4, std
511 GeV/c in the 100k-event sample) — confirmed by binning evenly in
1/pT and finding equal track counts per bin (~6475-6846 out of 100k
across 15 bins spanning the full range). So "1.5GeV" names the pT
floor of the generator, not a fixed beam momentum; this matters for how
the curvature/pT plots below are read, and is the basis for the
1/pT-binned tracking-efficiency measurement in step 4. Since
`load_hits()` flattens across all entries by default, this large
multi-event file loads exactly the same way as the merged BIB file, no
special-casing needed.

**Physical picture**: each simulated MC collision is expected to contain
the full BIB background (the merged plus+minus file, taken as one
collision's worth) *plus* the hits from a single signal muon (one event
from the muon-gun file). So BIB and signal are two contributions to the
*same* collision, not two separate, independently-normalized samples —
this is what fixes the normalization used everywhere below.

`signal_overlay_plots.py` produces the spatial r-z hit map (BIB hits
colored by subsystem, signal hits overlaid in black) plus a numeric
comparison CSV:
```
python3 signal_overlay_plots.py <bib.root> <signal.root> [output_dir] [geometry_dir]
```
run as `signal_overlay_plots.py ntu_bib_2evt.root
ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root
Analysis/step3_signal_overlay`. (An earlier version of this script also
drew three "hit density per subsystem/layer/disk" bar charts with a
signal bar alongside BIB mean/peak — dropped: BIB density already spans
~5-9 orders of magnitude on its own, and the per-event signal density
added another few, stretching the plots without adding information the
CSV doesn't already carry; the plain BIB-only density plots already
exist from step 1.)

Outputs, in `Analysis/step3_signal_overlay/`:
- `rz_map_with_signal.png` — the r-z hit map (as in step 1) with BIB hits
  colored by subsystem and signal hits overlaid in black (subsampled,
  60k points). Because the signal sample spans nearly the full polar
  angle range (10°-170°), it densely traces the geometric outline of
  every layer/disk in the tracker; BIB hits show up as color concentrated
  toward small r/z (near the beam) beneath that outline.
- `signal_vs_bib_summary.csv` — per-subsystem BIB n_hits/mean/peak
  density (hits/mm², whole BIB sample) vs. signal n_hits/n_events/
  mean-per-event density (hits/mm² from one signal muon) — kept as a
  numeric reference table even though no bar chart of it is produced
  any more.

**Signal overlay on the incidence-angle/curvature/time plots**
(`signal_overlay_angle_plots.py`), producing 6 "with signal" companions
to the step 2 plots:
```
python3 signal_overlay_angle_plots.py <bib.root> <signal.root> [output_dir]
```
run as `signal_overlay_angle_plots.py ntu_bib_2evt.root
ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root
Analysis/step3_signal_overlay` (same output folder as
`signal_overlay_plots.py`, so the r-z map/CSV and these 6 plots live
together; under `./run_all` both scripts simply write into that one
folder).

**Normalization: hits / collision / bin.** Per the physical picture
above, both curves are put on the same, addable footing:
- BIB drawn as raw hit counts per bin (unweighted) — the merged
  plus+minus file already represents one collision's full background.
- signal drawn as hit counts per bin weighted by `1/n_sig_events` — the
  average contribution of the one signal muon per collision.
(An earlier version of this script instead normalized each histogram to
unit area, i.e. a pure shape comparison discarding both samples' actual
scale — replaced, since what's wanted is how the two actually add up
per collision, not just where each one's shape peaks.) Produces, one
panel per subsystem, BIB (blue) and signal (green) step outlines
overlaid, y-axis "hits / collision / bin", log scale:
- `z_axis_intercept_per_subsystem_with_signal.png` (full ±3000mm range)
- `z_axis_intercept_per_subsystem_zoom_with_signal.png` (±100mm zoom)
- `inv_radius_per_subsystem_with_signal.png` (full range, pT-relabeled)
- `inv_radius_per_subsystem_zoom_with_signal.png` (|pT|>=1 GeV/c zoom)
- `time_corrected_per_subsystem_with_signal.png` (full ±20ns range of
  `t_corrected_ns`, same convention as `time_of_flight_plots.py`)
- `time_corrected_per_subsystem_zoom_with_signal.png` (±2ns zoom)

**Result (angle/curvature)**: BIB dominates signal by roughly 3-5 orders
of magnitude per collision across most of both quantities' range (e.g.
in the pT zoom, BIB sits around 10³-10⁵ hits/collision/bin vs. signal's
10⁻³-10⁻¹) — the honest per-collision picture, in contrast to the
shape-only version which made the two look comparable by discarding
scale. Signal is much more sharply concentrated than BIB in both
quantities, as expected for real muons from the
nominal IP vs. diffuse BIB secondaries: signal's small per-collision
contribution is a narrow spike at high pT in the full-range curvature
plot and a broad plateau across the whole ±1 GeV/c zoom window
(consistent with the sample's flat-in-1/pT generation, which populates
every pT scale roughly equally rather than concentrating at one value), while BIB spreads much more
broadly toward low pT/high curvature. For z-axis intercept, signal's
tiny per-collision peak sits tightly around 0 (the true IP vertex,
±1.5mm generator spread) and lands almost exactly where BIB itself
shows its unexplained dip at z=0 in the zoomed view — consistent with
that BIB dip being a genuine transport feature (most BIB secondaries'
local trajectories don't extrapolate back through z=0 itself) rather
than an artifact, though still not fully explained.

**Result (time-of-flight-corrected time)**: in the same
hits/collision/bin convention, BIB (blue) sits broadly across the whole
±2ns zoom window at roughly 10³-10⁵ hits/collision/bin with no sharp
feature at zero, while signal (green) is a tiny (~10⁻⁴-1
hits/collision/bin) but very sharp spike exactly at `t_corrected_ns=0`
in every subsystem — the same qualitative BIB-vs-signal discrimination
seen in the standalone per-file comparison (see "Time-of-flight
correction" above), now on the same addable per-collision footing as
the other overlay plots. Confirms `t_corrected_ns` is a genuine
candidate cut variable in the same sense as the angle/curvature
variables, not just a shape difference.

## Selection cuts (step 4)
Three quantities were picked as cut variables, each of which should sit
near a signal-consistent "expected" value for a genuine on-time IP muon
and is spread much more broadly for BIB (see their sections above):
`t_corrected_ns` (peaks at 0 for signal), `z_axis_intercept_mm` (peaks at
0, the true IP, for signal), and the curvature-derived transverse
momentum (large/straight for a moderate-to-high-momentum signal muon,
small for the sharply-curved low-momentum BIB secondaries that dominate
the background).

**Editable config, not hardcoded values**: cut values live in
`__cuts_config.txt` (plain INI, stdlib `configparser`), one section per
cut, each with a `halfwidth` (and `enabled` flag) so a hypothesis can be
changed and re-run with no code edits:
```ini
[t_corrected_ns]
halfwidth = 0.3         # ns
enabled = true
zoom_halfwidth = 2.0     # ns - display range of the N-1 zoomed plot only

[z_axis_intercept_mm]
halfwidth = 15.0         # mm
enabled = true
zoom_halfwidth = 100.0   # mm - display range of the N-1 zoomed plot only

[momentum_gev]
halfwidth = 5.0          # GeV/c (see note below on how this is applied)
enabled = true
zoom_halfwidth = 1.0     # GeV/c - display range of the N-1 zoomed plot only

[track]
min_hits_found = 5           # min. surviving hits for a track to count as "found" (step 4 part 2)
exclude_vertex_hits = false  # if true, vertex-detector hits don't count toward min_hits_found
```
Loaded/applied by `bib_common.load_cuts(path)` / `bib_common.apply_cuts(hits, cuts)`.
`t_corrected_ns` and `z_axis_intercept_mm` are simple symmetric windows,
accept `|value| <= halfwidth`. `momentum_gev` is **not** a symmetric
window on pT itself - pT = 0.3·B/|1/R| diverges as R → Infinity (a
perfectly straight, best-reconstructed track), and worse, pT is
*discontinuous* at 1/R = 0 (it jumps between +Infinity and -Infinity
depending which side of R = Infinity you approach from), so no finite
window in pT-space can ever include that case - a real gap the first
version of this cut had, caught before running it at scale. Instead the
cut is applied directly on the curvature `inv_radius_per_mm` itself
(bounded and continuous everywhere): accept
`|1/R| <= (0.3*B_FIELD_T/1000) / halfwidth_gev`, i.e. "reconstructed
|pT| >= halfwidth_gev", with R = Infinity always accepted since 1/R = 0
sits at the exact center of that window. `apply_cuts()` returns both the
combined accept mask and each cut's own mask (for a cutflow breakdown).

**`zoom_halfwidth`** (per cut) is purely a plot display range - the
half-width of that variable's zoomed histogram in `n1_cut_plots.py`
(`bib_common.ZOOM_HALFWIDTH_DEFAULTS` supplies the value above as the
fallback when a cut's section doesn't set it, so older configs without
this key still work unchanged). It has no effect on the cut itself; it
exists so the zoomed view can be widened to keep the cut-threshold line
visible after loosening that cut's `halfwidth`, without a code change.

**`exclude_vertex_hits`** (`[track]` section): if true, hits in the
vertex detector (VXD barrel/endcap) are excluded from the per-track
surviving-hit count that `min_hits_found` is compared against in
`track_efficiency.py`, even when those hits pass the three selection
cuts - i.e. "found" then means `min_hits_found`+ surviving hits *outside*
the vertex detector. Default false (vertex hits count the same as any
other subsystem's).

**Running step 4.** Step 4 (`apply_cuts.py`, `track_efficiency.py`,
`n1_cut_plots.py`) runs as part of `./run_all` together with steps 1-3
(see "How to run the analysis" at the top); every run is archived with
the exact settings that produced it in `Analysis/runs/<date>_<time>/`.
Steps 1-3 are re-run too, since the measurement smearing affects them
and the step-3 overlay plots draw the current cut thresholds.

**`apply_cuts.py`** (step 4, part 1) runs the configured cuts on both BIB
and signal and reports the comparison with vs. without cuts - the point
of this step:
```
python3 apply_cuts.py <bib.root> <signal.root> [cuts_config] [output_dir] [geometry_dir]
```
run as `apply_cuts.py ntu_bib_2evt.root
ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root __cuts_config.txt
Analysis/step4_cuts`. Outputs, in `Analysis/step4_cuts/`:
- `cutflow_bib.csv` / `cutflow_signal.csv` - per subsystem (and an "ALL"
  row): hit counts passing each individual cut on its own, and all three
  combined, plus the combined pass fraction.
- `density_before_after_cuts.csv` - per subsystem: **BIB only** n_hits/
  mean density/peak density, before and after cuts, plus BIB rejection
  fraction (using the true geometry area, same as step 1/3). Signal hit
  density was dropped from this table and plot per the user's
  instruction - it isn't an interesting quantity on its own (signal
  efficiency is instead measured properly by track, see
  `track_efficiency.py` below); the cutflow CSV above still reports raw
  signal hit counts/pass fractions.
- `density_before_after_cuts.png` - single panel (BIB mean density with
  vs. without cuts), log-scale bar chart per subsystem.

**Result, starting hypothesis** (±0.3ns, ±15mm, |pT|>=5 GeV/c, all three
combined): BIB is suppressed by 99.3-99.98% depending on subsystem
(99.73% overall), while signal *hits* (not yet tracks, see below)
surviving all three cuts run 25-29% depending on subsystem (27.1%
overall) - i.e. this starting cut set removes roughly 2-4 orders of
magnitude of BIB per subsystem while keeping about a quarter of genuine
signal hits. The momentum cut is the single biggest driver of BIB
rejection in every subsystem (BIB is dominated by low-momentum,
sharply-curved secondaries); the time and z-intercept cuts contribute
less on their own but tighten the combined selection further. Not yet
iterated on - this is the starting point the user specified;
`__cuts_config.txt` is designed to make trying other hypotheses
(tighter/looser windows, disabling one cut) a config edit and a re-run,
not a code change.

**`track_efficiency.py`** (step 4, part 2) measures the quantity that
actually matters for physics: not just what fraction of signal *hits*
survive the cuts, but what fraction of signal *tracks* are still
reconstructible afterward. Per the user's definition: for each
simulated muon-gun track (one event = one generated muon), count how
many of its own hits survive the three cuts; the track counts as
"found" if that count is >= `min_hits_found` (`__cuts_config.txt`,
`[track]` section, default 5). If `[track]`'s `exclude_vertex_hits` is
true, surviving hits in the vertex detector (VXD barrel/endcap) don't
count toward that total, so "found" then requires `min_hits_found`+
surviving hits *outside* the vertex detector - useful for asking
whether a track is still findable from its outer-tracker hits alone.
Default false (vertex hits count the same as any other subsystem's).
```
python3 track_efficiency.py <signal1.root> [signal2.root ...] [--cuts __cuts_config.txt] [--out output_dir] [--bins N]
```
run as `track_efficiency.py
ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root --out
Analysis/step4_track_efficiency`. Because the sample is generated flat
in 1/pT rather than at a fixed pT (see "Signal sample and overlay"
above), tracks are grouped into bins evenly spaced in **1/pT** (not
pT) - this keeps roughly equal statistics per bin all the way into the
high-pT tail, whereas equal-width bins in pT itself would leave the
tail nearly empty. Efficiency = (# tracks found)/(# tracks examined)
is computed per bin, each with a binomial uncertainty. Default binning
is **150 bins** (`N_BINS_DEFAULT`, `--bins` to override) - increased
from an initial 15 to better resolve the sharp turn-on described below.
Outputs, in `Analysis/step4_track_efficiency/`:
- `track_efficiency_vs_pt.csv` - per bin: 1/pT range/center, pT
  range/mean, n tracks examined/found, efficiency, uncertainty.
- `track_efficiency_vs_pt.png` - efficiency (%) vs. 1/pT bin center,
  with vertical (binomial) and horizontal (bin half-width) error bars.
  **Axis convention**: linear in **1/pT**, not pT and not logarithmic -
  matching the variable the binning is actually evenly spaced in (same
  convention as `inv_radius_per_subsystem.png` in step 2). Tick
  *positions* stay evenly spaced in 1/pT; only the printed tick
  *labels* are relabeled as the reciprocal pT value in GeV/c
  (`pt_ticks_for_inv_pt_axis()`), so the axis reads directly in
  physical pT while the underlying spacing (and the visible point
  density) still reflects the flat-in-1/pT generation. Tick
  *positions* are chosen so the printed pT labels land on round-ish
  numbers (5, 10, 20, 50, 100, 200 GeV/c, ...) rather than the
  arbitrary values evenly-spaced-in-1/pT ticks would print; because
  round numbers spaced evenly in pT pile up on top of each other near
  1/pT=0 (the high-pT end) on a linear-in-1/pT axis, candidates are
  walked outward from the low-pT/turn-on end and a candidate is kept
  only once it's far enough (in 1/pT) from the last kept tick - so the
  flat, uninformative high-pT region above the last kept tick is left
  unlabeled rather than crowded with overlapping numbers. Short,
  unlabeled minor ticks then fill in *each gap between consecutive
  kept major ticks* at a round step sized to divide that gap into
  roughly 5 parts (6, 7, 8, 9 between 5 and 10 GeV/c; 12, 14, 16, 18
  between 10 and 20 GeV/c; 60, 70, 80, 90 between 50 and 100 GeV/c;
  etc., via `_nice_minor_step`) - built gap-by-gap like this, rather
  than from one global "1-9 x 10^k" candidate set, so every major
  interval gets minor ticks (a plain tens/hundreds sequence would
  otherwise skip intervals like 10-20 entirely) and none spill out
  past the last kept major tick. Y-axis has
  headroom above 100% (`ylim(0, 108)`) so points sitting at/near full
  efficiency aren't clipped against the top edge. X-axis is cropped
  just past the last bin (scanning from high pT down) with efficiency
  > 1% - the long flat near-zero tail extending down to the pT=1.5
  GeV/c generator floor carries no information and was wasting most of
  the plot's horizontal space, so it's dropped from the view (still in
  the CSV).

**Result**: a sharp turn-on right at the momentum-cut threshold. With
the finer (150-bin) binning: ~0% efficiency for pT well below
threshold, a jump to ~22% in the bin just below 5 GeV/c, ~92% in the
bin just above, and ~100% for pT above that - i.e. the transition is
now resolved into a few intermediate points rather than the single
coarse [4.5, 5.6] GeV/c transition bin seen with the original 15-bin
version. This is an expected consequence of the momentum cut's current
halfwidth (5 GeV/c, applied as a hard `|pT| >= 5 GeV/c` requirement,
see "Selection cuts" above) - a track below threshold has essentially
none of its hits surviving that one cut, so `min_hits_found` is almost
never reached, while a track above threshold has essentially all of
them survive. The efficiency curve is still close to a step function
rather than a gradual turn-on even at this finer resolution, which is
worth keeping in mind when choosing whether to loosen the momentum
halfwidth in a future iteration (see "Next step").

**Memory note**: `apply_cuts.py` is the first script to combine
`add_incidence_angles` + `add_time_of_flight` (many intermediate per-hit
arrays) with `region_table`/`add_peak_density` (which needs the full hit
arrays again) in one process, on the full 16M-hit BIB sample - this
combination was enough to get OOM-killed on the Mac on the first attempt.
Fixed by dropping per-hit fields no longer needed as soon as they're
consumed (`slim_hits()`/inline `del`s in `apply_cuts.py`, plus explicit
`gc.collect()` calls) - worth keeping in mind for any future script that
similarly stacks several of these per-hit-array-heavy steps together in
one run.

**`n1_cut_plots.py`** (step 4, part 3) produces "N-1" cut-validation
plots — the standard technique of plotting a cut variable's own
distribution after applying every *other* cut but not the one being
plotted, so the cut threshold can be judged against an unbiased view of
where signal and BIB actually separate, rather than a view already
shaped by that same cut. With 3 cuts, N-1 means 2 of the 3 applied at a
time: the z-axis-intercept plots apply the time and momentum cuts (not
z); the momentum/curvature plots apply the time and z cuts (not
momentum); the time plots apply the z and momentum cuts (not time). A
disabled cut in `__cuts_config.txt` simply drops out of its own N-1
combinations (`bib_common.apply_cuts()` already returns an all-True
mask for a disabled cut), so this still does the right thing with
fewer than 3 cuts turned on.
```
python3 n1_cut_plots.py <bib.root> <signal.root> [cuts_config] [output_dir]
```
run as `n1_cut_plots.py ntu_bib_2evt.root
ntu_muongun_pt1p5GeV_theta10-170_phi0-360_dz1p5_100k.root __cuts_config.txt
Analysis/step4_n1_cuts`. Reuses `signal_overlay_angle_plots.py`'s 6-plot
layout and hits/collision/bin BIB-vs-signal convention (z-axis intercept
full+zoom, curvature/pT full+zoom, corrected-time full+zoom - see
"Signal overlay on the incidence-angle/curvature/time plots" above), so
the N-1 plots read the same way as their unfiltered step-3 counterparts,
just with the other two cuts applied. Each plot also draws the current
(enabled) cut's own threshold as a vertical dashed line, so the
threshold's position relative to the signal/BIB separation can be
judged directly. Each variable's zoomed plot uses that cut's own
`zoom_halfwidth` from `__cuts_config.txt` (see "Selection cuts" above) as
its display range, rather than a fixed value, so the zoomed view can be
widened alongside a loosened `halfwidth` without the threshold line
falling outside it. Outputs, in `Analysis/step4_n1_cuts/`:
`z_axis_intercept_per_subsystem_n1.png`/`_zoom_n1.png`,
`inv_radius_per_subsystem_n1.png`/`_zoom_n1.png`,
`time_corrected_per_subsystem_n1.png`/`_zoom_n1.png`.

## Tooling
ROOT is not installed in the analysis sandbox; the analysis is done with
Python (uproot/awkward/numpy/matplotlib), which reads the same ROOT
ntuples and gives identical results. Code lives in
`~/code/MuonCollider/` on the user's Mac:
- `bib_common.py` — `decode_id0`, `load_hits` (flattens across all
  entries of a tree by default — see "Merged file" above; also computes
  a per-hit `event_id` field identifying which event/track each hit
  came from, used by `track_efficiency.py` to group hits back into
  tracks), `region_table`,
  `add_peak_density` (adaptive per-region bin sizing),
  `add_incidence_angles` (longitudinal/transverse track-incidence
  angles, plus z_axis_intercept_mm and the transverse curvature
  inv_radius_per_mm/psi_transverse_deg), `add_time_of_flight`
  (time-of-flight-corrected hit time, described above; requires
  `add_incidence_angles` to have been called first; also adds the
  signed curvature-derived `pT_curv_gev`), `load_cuts`/`apply_cuts`
  (selection cuts, described above), `load_track_params` (reads the
  `[track]` section of `__cuts_config.txt`, currently just
  `min_hits_found`), `mask_hits` (filter every per-hit
  array field of a hits dict by a boolean mask, keeping metadata
  fields unchanged - used to re-run region_table/add_peak_density on
  only the hits passing cuts), `subsystem_density_table`,
  `write_logfile`, `prepare_output_dir` (when a script is run by hand,
  archives old results into `OldResults/` instead of overwriting;
  skipped under `./run_all`, which gives every run its own folder),
  `default_cuts_config` (the cuts file a script uses when run by hand
  without one: `./__cuts_config.txt` if present, else the template), and
  the smearing functions `load_smearing_config`,
  `apply_position_time_smearing`, `apply_angle_smearing`,
  `describe_smearing`.
- `geometry.py` — parses the true detector geometry XML to get exact
  sensitive area per region (`build_area_lookup`,
  `annotate_rows_with_geometry_area`).
- `merge_bib_files.py` — merges single-entry BIB files (plus, minus,
  ipp) into one combined file, one entry per input (see "Merged file"
  above); its `merge()` is what `./run_all` uses to (re)build the
  combined file, via `input_files.py`.
- `input_files.py` — reads `__input_files_config.txt` (paths of all input
  files), and checks/rebuilds the combined BIB file
  (`ensure_combined`).
- `make_basic_plots.py` — step 1 main script; produces per-subsystem/
  layer/disk density plots (mean vs. peak, log y-axis, disks ordered
  by signed z), an r-z hit map, `region_summary.csv`,
  `subsystem_summary.csv`, and `run_log.txt`. Peak-density target
  statistics and percentile are set at the top of the file
  (`PEAK_TARGET_HITS_PER_BIN`, `PEAK_PERCENTILE`).
- `incidence_angle_plots.py` — step 2 (part 1, diagnostic) script;
  produces the incidence-angle, z-axis-intercept (full range + zoom),
  and transverse-curvature (full range + zoom) histograms and summary
  CSV described above; run alongside step 1 as standard output. Does
  not yet apply any cuts or touch hit density. `B_FIELD_T` controls the
  pT axis label/binning; `Z0_ZOOM_RANGE_MM` and `PT_ZOOM_RANGE_GEV`
  control the two zoomed-window widths (the pT zoom's actual curvature
  cutoff is derived from `PT_ZOOM_RANGE_GEV` via `GEV_PER_INV_M`, so
  the zoom stays in the same 1/R-centered convention as the full plot).
  Also exports `pt_ticks_for_axis`, `panel_grid`, and its constants for
  reuse by `signal_overlay_angle_plots.py`.
- `time_of_flight_plots.py` — step 2 (part 3) script; produces the raw
  and TOF-corrected hit-time histograms and summary CSV described in
  "Time-of-flight correction" above. Also exports `TC_RANGE_NS` and
  `TC_ZOOM_RANGE_NS` for reuse by `signal_overlay_angle_plots.py`.
- `signal_overlay_plots.py` — step 3 script; produces the r-z hit map
  with signal overlaid and the BIB-vs-signal density summary CSV,
  described above (no longer produces density bar charts with signal).
- `signal_overlay_angle_plots.py` — step 3 companion script; overlays
  the signal sample on the 6 incidence-angle/curvature/time-of-flight
  diagnostic plots (z-axis intercept full+zoom, transverse curvature/pT
  full+zoom, TOF-corrected time full+zoom), in hits/collision/bin,
  described above.
- `templates/__cuts_config.txt`, `templates/__smearing_config.txt` —
  templates of the two settings files; the live, editable copies are in
  the Analysis folder (see "How to run the analysis" at the top).
- `check_configs.py` — checks the two settings files before a run and
  prints a readable summary of them; `./run_all` refuses to start if it
  finds a problem.
- `apply_cuts.py` — step 4 (part 1) script; applies the configured cuts
  to BIB and signal and reports the with-cuts-vs-without-cuts
  comparison (BIB density only - signal density dropped per the user's
  instruction), described above.
- `track_efficiency.py` — step 4 (part 2) script; measures
  track-finding efficiency vs. generated pT, binned evenly in 1/pT to
  match the signal sample's flat-in-1/pT generation, described above.
- `n1_cut_plots.py` — step 4 (part 3) script; produces N-1 cut-
  validation plots (each cut variable plotted with the other cuts
  applied, not its own), reusing `signal_overlay_angle_plots.py`'s
  layout/conventions, described above.
- `run_all_steps.sh` — the full pipeline (steps 1-4) behind the
  `./run_all` launcher in the Analysis folder (see "How to run the
  analysis" at the top). Replaces the former `run_step4.sh`.

Step 1, step 2 (both the angle and time-of-flight scripts) are run for
the original plus/minus files AND the merged file on every run, plus
once on the signal file for the time-of-flight validation; step 3 (both
scripts) and step 4 use the merged BIB file plus the signal file. Each
run's results go to
`~/Dropbox/Documents/MuonColliderSimulation/Analysis/runs/<date>_<time>/`
(`step1_basic_plots/`, `step1_basic_plots_minus/`,
`step1_basic_plots_combined/`, `step2_incidence_angles/`,
`step2_incidence_angles_minus/`, `step2_incidence_angles_combined/`,
`step2_time_of_flight/`, `step2_time_of_flight_minus/`,
`step2_time_of_flight_combined/`, `step2_time_of_flight_signal/`,
`step3_signal_overlay/` (shared by both step-3 scripts), `step4_cuts/`,
`step4_track_efficiency/`, `step4_n1_cuts/`, plus the two settings
files, `run_log.txt`, `code_version.txt` and `summary.txt`); the same
`step*` folders (and `_highlights/`) directly in `Analysis/` always hold
the latest successful run. `runs/` also keeps the older step-4-only archives made before this
workflow existed (suffix `_step4-only`, each with a `NOTE.txt`).

## Next step (in progress)
The step-4 cut framework is in place, the starting hypothesis has been
run (see "Selection cuts" above: ~99.7% BIB rejection overall), and
track-finding efficiency vs. pT is now measured directly
(`track_efficiency.py`: ~0% below ~4.5 GeV/c, ~51% at [4.5,5.6] GeV/c,
~100% above ~5.6 GeV/c). That result - a near-step-function turn-on
right at the 5 GeV/c momentum-cut threshold - is itself a candidate
next topic: it's a direct consequence of how hard the current momentum
cut is, and loosening its halfwidth would likely trade some BIB
rejection for a smoother, more physically informative efficiency curve
(worth discussing with the user before changing it). More generally,
next: iterate on the cut values in `__cuts_config.txt` (tighter/looser
windows, or disabling individual cuts) - `./run_all` regenerates
everything that depends on them in one command - to see how the
BIB-rejection/track-efficiency trade-off moves, and decide with the
user whether/how to combine the three cuts differently (e.g. an
optimization rather than three independently-chosen windows). Also
open: the dip at z=0 in the zoomed z_axis_intercept plots, and the
asymmetric shape of the BIB corrected-time distribution (see
"Time-of-flight correction" above) are each plausibly explained but
not yet decomposed into their underlying sub-populations. The N-1
cut-validation plots (`n1_cut_plots.py`) are now available to help
judge each cut threshold against an unbiased view of the signal/BIB
separation as this iteration proceeds.
