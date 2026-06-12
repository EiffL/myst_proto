---
title: Configuration-Space BAO Distances from DESI DR1
subtitle: A reproduction across eight tracer–redshift bins, pre- and post-reconstruction
authors:
  - name: DESI Collaboration
exports:
  - format: pdf
---

## Abstract

We reproduce the DESI DR1 configuration-space baryon-acoustic-oscillation (BAO)
measurement across all eight tracer–redshift bins — BGS, LRG1, LRG2, LRG3,
ELG1, ELG2, the combined LRG3+ELG1 bin, and QSO, together spanning
$0.1 < z < 2.1$ — by fitting a damped-BAO template to the two-point correlation
function $\xi(s)$ before and after density-field reconstruction. Working in the
$\alpha_\mathrm{iso}$–$\alpha_\mathrm{AP}$ dilation basis (with
$\alpha_\mathrm{iso}$ alone for the 1D-fit tracers BGS, ELG1, and QSO), the
baseline configuration recovers a well-defined acoustic feature in every
post-reconstruction tracer, and reconstruction sharply tightens the isotropic
dilation: for ELG1, $\sigma(\alpha_\mathrm{iso})$ falls from
{astra:value}`bao_alpha_values tracer=elg1 recon=Pre col=alpha1_std` to
{astra:value}`bao_alpha_values tracer=elg1 recon=Post col=alpha1_std`. The
strongest single constraint, the combined LRG3+ELG1 bin, reaches
$\alpha_\mathrm{iso} =$
{astra:value}`bao_alpha_values tracer=lrg3_elg1 recon=Post col=alpha1_mean err=alpha1_std`.
Propagated to distances, we obtain $D_M/r_d$, $D_H/r_d$, and $D_V/r_d$ at the
percent level per tracer — for LRG3+ELG1, $D_V/r_d =$
{astra:value}`bao_distance_table tracer=lrg3_elg1 col=DV_over_rd pm` at
$z_\mathrm{eff} =$ {astra:value}`bao_distance_table tracer=lrg3_elg1 col=z_eff` —
in good agreement with the Planck 2018 $\Lambda$CDM standard ruler.

## 1. Introduction

In the hot plasma of the early Universe, sound waves launched by primordial
overdensities propagated outward until recombination released the photons and
the waves stalled. The distance each wave had travelled — the sound horizon at
the drag epoch, $r_d$ — was thereby frozen into the clustering of matter as a
preferred comoving separation, visible today as a localised peak in the galaxy
two-point correlation function $\xi(s)$ near $100\,h^{-1}$ Mpc. Because $r_d$
is calibrated to a quarter of a percent by the cosmic microwave background,
this baryon-acoustic-oscillation (BAO) feature is the most robust standard
ruler in large-scale structure: measuring its apparent size across and along
the line of sight yields the comoving distance $D_M(z)/r_d$ and the Hubble
distance $D_H(z)/r_d$, mapping the expansion history of the Universe.

The ruler is not pristine. Over cosmic time, bulk flows and non-linear
structure growth displace galaxies from their original positions, which
{astra:prior-insight}`nonlinear_peak_smearing|blurs the acoustic peak` and
degrades the achievable distance precision by a
{astra:prior-insight}`precision_loss_factor_three|factor of roughly three`.
Density-field reconstruction was introduced to undo this damage: by estimating
the large-scale displacement field from the observed density and moving
galaxies back along it, reconstruction
{astra:prior-insight}`recon_sharpens_bao_peak|re-sharpens the acoustic peak`
and {astra:prior-insight}`recon_factor_two_demonstrated|recovers roughly a factor of two in distance precision`.

This article reproduces the DESI DR1 BAO measurement in configuration space,
fitting every tracer both before and after reconstruction so the gain is
measured rather than assumed; Fourier-space $P(k)$ is out of scope. The
headline result is twofold: the acoustic feature is
{astra:finding}`bao_detected_post_recon` in every tracer after reconstruction,
and reconstruction {astra:finding}`recon_reduces_alpha_iso_sigma` — delivering
{astra:finding}`subpercent_alpha_iso_precision` for the best-measured LRG bins.

## 2. Data

The measurement consumes three classes of input. The first is the raw DESI DR1
LSS clustering catalogs — galaxy positions plus eighteen random realisations
per Galactic cap, which define the survey geometry and selection. The second
is the tabulated DESI fiducial cosmology
({astra:prior-insight}`planck2018_headline_parameters|Planck 2018 ΛCDM`),
which converts redshifts to comoving distances; every apparent BAO scale is
measured relative to this fiducial expectation. The third is the set of
RascalC semi-analytic covariance matrices for $\xi(s)$ — one pre- and one
post-reconstruction file per tracer, taken from the published DESI release.

The covariances are ingested as published rather than recomputed, and this is
the load-bearing data decision: because each RascalC matrix is calibrated
against the fiducial DESI reconstruction and clustering configuration, adopting
it fixes the $\xi(s)$ binning to the published linear grid and, as discussed
below, pins several upstream pipeline choices to their fiducial settings.

## 3. Methods

The pipeline runs in three stages: {astra:analysis}`reconstruction` produces
shifted catalogs, {astra:analysis}`clustering` measures correlation functions
from them, and a template-fitting stage turns each correlation function into
posterior constraints on the BAO scale — one MCMC chain per (tracer,
reconstruction state), sixteen chains in total.

**Reconstruction.** The linear displacement field is estimated from the
Gaussian-smoothed galaxy density and applied symmetrically to galaxies and
randoms (the {astra:decision}`reconstruction.convention` choice, RecSym),
moving structure approximately back to its initial position and thereby
sharpening the acoustic peak. The single load-bearing knob is the
{astra:decision}`smoothing_radius`: it sets the scale $\Sigma_\mathrm{sm}$ of
the smoothing applied before the displacement solve, and because the same
scale re-enters the BAO damping template downstream, the reconstruction stage
inherits it from the root of the analysis rather than choosing independently —
the two stages cannot drift apart.

**Clustering.** We measure the Landy–Szalay $\xi(s,\mu)$ multipoles on the raw
catalogs (pre-reconstruction) and on the shifted catalogs
(post-reconstruction), per tracer redshift slice. Because the covariance is
pinned to the published RascalC grid (the
{astra:decision}`covariance_source` decision), the $s$-binning, $\mu$-binning,
and estimator are locked to that grid, and the non-fiducial smoothing options
are declared incompatible with it — the baseline universe only validates at
fiducial smoothing. The one remaining free choice here,
{astra:decision}`clustering.imaging_weights`, is null-tested by the same
covariance at sub-$\sigma$ significance and left unlocked.

**Template fitting.** Each chain fits a damped-BAO template with `desilike` +
`emcee`. The template is built in the fiducial cosmology and the fit measures
how far the acoustic feature in the data is dilated away from it:
$\alpha_\mathrm{iso}$ rescales the feature isotropically and traces
$D_V/r_d$, while $\alpha_\mathrm{AP}$ warps it anisotropically and traces the
Alcock–Paczyński ratio $D_M/D_H$. The five 2D tracers (the LRG bins, ELG2, and
LRG3+ELG1) fit both parameters from the monopole and quadrupole; the three
sparser 1D tracers (BGS, ELG1, QSO) fit $\alpha_\mathrm{iso}$ alone from the
monopole. The smooth, BAO-free part of each multipole is absorbed by the
{astra:decision}`broadband` model — the spline form is fiducial and has been
shown to {astra:prior-insight}`spline_broadband_fiducial|give α consistent with the polynomial form`,
with the residual absorbed into the
modelling-error budget. The non-linear smearing of the peak is modelled by
damping parameters controlled jointly by the
{astra:decision}`damping_prior` and {astra:decision}`damping_centers`,
anchored to the result that
{astra:prior-insight}`damping_prior_width_justification|mis-centred damping priors bias α`.
Three further template-shape choices —
{astra:decision}`fog_placement`, {astra:decision}`dilate_model_component`, and
{astra:decision}`dewiggling_method` — are pinned to fiducial and exposed to
document the modelling-systematic budget they underlie; the fiducial dilation
acts on the wiggle component only, the expected behaviour given that
reconstruction {astra:prior-insight}`recon_reduces_bao_damping|reduces the non-linear BAO damping`.

## 4. Results

### 4.1 Detection and peak sharpening

The most direct view of the measurement is the acoustic feature itself.
[](#output-bao_fit_plot) isolates it by subtracting the smooth part of the
best-fit model from each measured multipole: in every tracer the
post-reconstruction peak is visibly narrower and better matched by the
template than its pre-reconstruction counterpart — the peak
{astra:finding}`bao_peak_sharpens_post_recon`.

:::{astra:output} bao_fit_plot
:::

To quantify whether the feature is detected at all, each post-reconstruction
correlation function is refit with the BAO wiggles removed, and
[](#output-bao_detection_plot) profiles the $\chi^2$ difference between the
two models as a function of $\alpha_\mathrm{iso}$. Every tracer develops a
well-defined minimum near the fiducial scale — weakest for the sparse 1D
tracers, strongest for the combined LRG3+ELG1 bin.

:::{astra:output} bao_detection_plot
:::

### 4.2 The dilation parameters

Condensing the chains, the post-reconstruction isotropic dilation is consistent
with unity across the suite — the acoustic scale in DESI DR1 sits where the
fiducial cosmology predicts. The strongest 2D bins reach the sub-percent
regime: $\alpha_\mathrm{iso} =$
{astra:value}`bao_alpha_values tracer=lrg3_elg1 recon=Post col=alpha1_mean err=alpha1_std`
for the combined LRG3+ELG1 bin and
{astra:value}`bao_alpha_values tracer=lrg3 recon=Post col=alpha1_mean err=alpha1_std`
for LRG3. The gain from reconstruction is largest where the
pre-reconstruction feature is most degraded — for ELG1,
$\sigma(\alpha_\mathrm{iso})$ contracts from
{astra:value}`bao_alpha_values tracer=elg1 recon=Pre col=alpha1_std` to
{astra:value}`bao_alpha_values tracer=elg1 recon=Post col=alpha1_std`, and for
LRG2 from {astra:value}`bao_alpha_values tracer=lrg2 recon=Pre col=alpha1_std`
to {astra:value}`bao_alpha_values tracer=lrg2 recon=Post col=alpha1_std`. The
full set of fits is collected in [](#output-bao_alpha_values):

:::{astra:output} bao_alpha_values
:::

The fits are statistically well-behaved — $\chi^2/\mathrm{dof}$
{astra:finding}`bao_fit_chi2_near_dof`. The combined LRG3+ELG1 bin returns
$\chi^2 =$ {astra:value}`bao_alpha_values tracer=lrg3_elg1 recon=Post col=chi2`
for {astra:value}`bao_alpha_values tracer=lrg3_elg1 recon=Post col=dof` degrees
of freedom; QSO is the least well-behaved, at $\chi^2 =$
{astra:value}`bao_alpha_values tracer=qso recon=Post col=chi2` for
{astra:value}`bao_alpha_values tracer=qso recon=Post col=dof`.

### 4.3 Cosmological distances

The dilation parameters convert directly into distances: each chain carries
$D_M/r_d$, $D_H/r_d$, and $D_V/r_d$ as derived parameters, so the constraints
in [](#output-bao_distance_table) inherit the full non-Gaussian shape of the
posteriors rather than a linearised propagation. For the combined LRG3+ELG1
bin at $z_\mathrm{eff} =$
{astra:value}`bao_distance_table tracer=lrg3_elg1 col=z_eff` we measure
$D_M/r_d =$
{astra:value}`bao_distance_table tracer=lrg3_elg1 col=DM_over_rd pm`,
$D_H/r_d =$ {astra:value}`bao_distance_table tracer=lrg3_elg1 col=DH_over_rd pm`,
and $D_V/r_d =$
{astra:value}`bao_distance_table tracer=lrg3_elg1 col=DV_over_rd pm`. The 1D
tracers contribute $D_V/r_d$ alone:
{astra:value}`bao_distance_table tracer=bgs col=DV_over_rd pm` at
$z_\mathrm{eff} =$ {astra:value}`bao_distance_table tracer=bgs col=z_eff` (BGS),
{astra:value}`bao_distance_table tracer=elg1 col=DV_over_rd pm` at
{astra:value}`bao_distance_table tracer=elg1 col=z_eff` (ELG1), and
{astra:value}`bao_distance_table tracer=qso col=DV_over_rd pm` at
{astra:value}`bao_distance_table tracer=qso col=z_eff` (QSO).

:::{astra:output} bao_distance_table
:::

Placed on a BAO Hubble diagram against 6dFGS, WiggleZ, SDSS DR16, and DES Y6
([](#output-hubble_diagram_plot)), the DESI DR1 distances trace the Planck 2018
$\Lambda$CDM prediction across the full redshift range, with the standard-ruler
anchor set by
{astra:prior-insight}`planck2018_rdrag_precision|the Planck 2018 sound horizon`.

:::{astra:output} hubble_diagram_plot
:::

## 5. Systematics and robustness

The reported uncertainties are not statistical-only: a combined modelling +
HOD + fiducial-cosmology systematic budget is folded into the quoted errors at
the aggregation step, toggled by the
{astra:decision}`systematic_error_treatment` decision as a table rebuild
rather than a refit — the chains themselves are untouched, and the
statistical-only errors are kept alongside for comparison. The budget is
anchored to the DESI companion modelling result that
{astra:prior-insight}`combined_systematic_budget|sets the combined systematic budget`,
dominated by the wiggle/no-wiggle split, FoG placement, and
template-dilation choices catalogued in the decision register.
Fiducial-cosmology dependence is swept through
{astra:decision}`template_cosmology` (a set of AbacusSummit grids), and the
remaining {astra:decision}`fit_range`, {astra:decision}`ells`, and
{astra:decision}`fitting_method` decisions exist for sensitivity tests; none
moves the baseline result beyond its quoted budget.

## 6. Conclusions

A configuration-space BAO analysis of DESI DR1 detects the acoustic feature in
all eight tracer bins post-reconstruction, reaches sub-percent isotropic
precision in the strongest LRG bins, and yields a self-consistent set of
$D_M/r_d$, $D_H/r_d$, $D_V/r_d$ distances at the percent level — reproducing
the DESI DR1 BAO distance ladder. Every numerical claim above is traceable to
a registered ASTRA finding, the decision that configured it, and the
materialised output product it summarises.
