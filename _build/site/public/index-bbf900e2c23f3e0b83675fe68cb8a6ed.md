---
title: Configuration-Space BAO Distances from DESI DR1
subtitle: A reproduction across eight tracer–redshift bins, pre- and post-reconstruction
authors:
  - name: DESI Collaboration
exports:
  - format: pdf
---

> Every italic-underlined term below is a live ASTRA object: hover a
> {astra:finding}`bao_detected_post_recon` finding or a
> {astra:decision}`covariance_source` decision for its full record, and the
> figures carry their own provenance. The prose is the author's; the evidence
> is pulled from `astra.yaml` and the materialised result products.

## Abstract

We reproduce the DESI DR1 configuration-space baryon-acoustic-oscillation (BAO)
measurement across all tracer–redshift bins — BGS, LRG1, LRG2, LRG3, ELG1,
ELG2, the combined LRG3+ELG1 bin, and QSO — by fitting a damped-BAO template to
the two-point correlation function $\xi(s)$ before and after density-field
reconstruction. Working in the $\alpha_\mathrm{iso}$–$\alpha_\mathrm{AP}$ basis
(with $\alpha_\mathrm{iso}$ alone for the 1D tracers BGS, ELG1, QSO), the
baseline configuration recovers a well-defined acoustic feature in every
post-reconstruction tracer and sharply tightens the isotropic dilation: for
ELG1, $\sigma(\alpha_\mathrm{iso})$ falls from
{astra:value}`bao_alpha_values tracer=elg1 recon=Pre col=alpha1_std` to
{astra:value}`bao_alpha_values tracer=elg1 recon=Post col=alpha1_std`
post-reconstruction. The strongest single constraint, the combined LRG3+ELG1
bin, reaches $\alpha_\mathrm{iso} =$
{astra:value}`bao_alpha_values tracer=lrg3_elg1 recon=Post col=alpha1_mean err=alpha1_std`.
Propagated to distances, we obtain $D_M/r_d$, $D_H/r_d$, and $D_V/r_d$ at the
percent level per tracer — for LRG3+ELG1, $D_V/r_d =$
{astra:value}`bao_distance_table tracer=lrg3_elg1 col=DV_over_rd pm` at
$z_\mathrm{eff} =$ {astra:value}`bao_distance_table tracer=lrg3_elg1 col=z_eff` —
consistent with the Planck 2018 $\Lambda$CDM standard ruler.

## 1. Introduction

The BAO feature imprinted in the galaxy correlation function at the
sound-horizon scale $r_d$ is the most robust standard ruler in large-scale
structure. Its utility is limited by non-linear evolution,
which smears the acoustic peak and degrades the achievable distance precision by
a factor of a few — a {astra:prior-insight}`precision_loss_factor_three|few-fold precision loss` that
density-field reconstruction was introduced to recover. We measure the BAO scale
in DESI DR1 in configuration space, both before and after reconstruction, and
propagate the result to cosmological distances. Fourier-space $P(k)$ is out of
scope.

The headline result is twofold: the acoustic feature is
{astra:finding}`bao_detected_post_recon` in every tracer after reconstruction,
and reconstruction {astra:finding}`recon_reduces_alpha_iso_sigma` — delivering
{astra:finding}`subpercent_alpha_iso_precision` for the best-measured LRG bins.

## 2. Data

The measurement consumes three input classes: the raw DESI DR1 LSS clustering
catalogs (galaxies plus eighteen random realisations per cap), the tabulated
fiducial cosmology supplying $z \to r(z)$, and the RascalC analytic covariance
matrices — one pre- and one post-reconstruction file per tracer. The covariances
are ingested as published rather than recomputed, which fixes the $\xi(s)$
binning to the published linear grid and, as
discussed below, pins several upstream choices to their fiducial settings. The
full input register is tabulated in the [analysis record](#data-products).

## 3. Methods

The pipeline runs in two upstream stages — {astra:analysis}`reconstruction`
then {astra:analysis}`clustering` — feeding a template-fitting stage that
produces one MCMC chain per (tracer, reconstruction state).

**Reconstruction.** The linear displacement field is estimated from the
Gaussian-smoothed density and applied symmetrically to galaxies and randoms
(the {astra:decision}`convention` choice, RecSym), which sharpens the acoustic
peak. The single load-bearing knob is the {astra:decision}`smoothing_radius`:
it sets the $\Sigma_\mathrm{sm}$ entering the displacement solve, and because the
same scale re-enters the BAO damping template downstream, the reconstruction
stage inherits it from the root rather than choosing independently.

**Clustering.** We measure the Landy–Szalay $\xi(s,\mu)$ multipoles on the raw
and shifted catalogs. Because the covariance is pinned to the published RascalC
grid (the {astra:decision}`covariance_source` decision), the binning and
estimator are locked to that grid and the non-fiducial smoothing options are
declared incompatible with it — so the baseline universe validates only at
fiducial smoothing. The one remaining free choice,
{astra:decision}`imaging_weights`, is null-tested by the same covariance at
sub-$\sigma$ significance and left unlocked.

**Template fitting.** Each chain fits the fiducial damped-BAO template with
`desilike` + `emcee`. The 2D tracers carry $\alpha_\mathrm{iso}$ and
$\alpha_\mathrm{AP}$ (with $D_M/r_d$, $D_H/r_d$ derived in-chain); the 1D tracers
fit $\alpha_\mathrm{iso}$ only. The dominant modelling dial is the
{astra:decision}`broadband` model — the spline form is fiducial and has been
shown to {astra:prior-insight}`spline_broadband_fiducial|give α consistent with the polynomial form`, with the residual
absorbed into the modelling-error budget. BAO damping is controlled jointly by
the {astra:decision}`damping_prior` and {astra:decision}`damping_centers`,
anchored to the result that {astra:prior-insight}`damping_prior_width_justification|mis-centred damping priors bias α`.
Three further template-shape choices — {astra:decision}`fog_placement`,
{astra:decision}`dilate_model_component`, and {astra:decision}`dewiggling_method`
— are pinned to fiducial and exposed to document the modelling-systematic
budget they underlie; the fiducial dilation acts on the wiggle component, the
expected behaviour given that reconstruction {astra:prior-insight}`recon_reduces_bao_damping|reduces the non-linear BAO damping`.

## 4. Results

### 4.1 Detection and peak sharpening

The isolated BAO feature, fit per tracer pre and post-reconstruction, is shown
in [](#output-bao_fit_plot). The peak {astra:finding}`bao_peak_sharpens_post_recon`,
and the corresponding $\Delta\chi^2(\alpha_\mathrm{iso})$ profiles against a
no-BAO reference ([](#output-bao_detection_plot)) develop a well-defined minimum
near the fiducial scale in every post-reconstruction tracer.

:::{astra:output} bao_fit_plot
:::

:::{astra:output} bao_detection_plot
:::

### 4.2 The dilation parameters

Condensing the chains, the post-reconstruction isotropic dilation is consistent
with unity across the suite. The strongest 2D bins reach the sub-percent regime:
$\alpha_\mathrm{iso} =$
{astra:value}`bao_alpha_values tracer=lrg3_elg1 recon=Post col=alpha1_mean err=alpha1_std`
for the combined LRG3+ELG1 bin and
{astra:value}`bao_alpha_values tracer=lrg3 recon=Post col=alpha1_mean err=alpha1_std`
for LRG3. The gain from reconstruction is largest where non-linear evolution is
worst — for ELG1, $\sigma(\alpha_\mathrm{iso})$ contracts from
{astra:value}`bao_alpha_values tracer=elg1 recon=Pre col=alpha1_std` to
{astra:value}`bao_alpha_values tracer=elg1 recon=Post col=alpha1_std`, and for
LRG2 from {astra:value}`bao_alpha_values tracer=lrg2 recon=Pre col=alpha1_std`
to {astra:value}`bao_alpha_values tracer=lrg2 recon=Post col=alpha1_std`. The
full $\alpha_\mathrm{iso}/\alpha_\mathrm{AP}$ table, with the $\chi^2$ and
degrees of freedom per fit, is:

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

Propagating the post-reconstruction chains to the distance basis yields the
constraints in [](#output-bao_distance_table). For the combined LRG3+ELG1 bin at
$z_\mathrm{eff} =$ {astra:value}`bao_distance_table tracer=lrg3_elg1 col=z_eff`
we measure $D_M/r_d =$
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
$\Lambda$CDM prediction, which sets the standard-ruler anchor through
{astra:prior-insight}`planck2018_rdrag_precision|the Planck 2018 sound horizon`.

:::{astra:output} hubble_diagram_plot
:::

## 5. Systematics and robustness

The reported uncertainties fold a combined modelling + HOD + fiducial-cosmology
systematic budget into the chains at the aggregation step, toggled by the
{astra:decision}`systematic_error_treatment` decision as a table rebuild rather
than a refit. The budget is anchored to the DESI companion modelling result that
{astra:prior-insight}`combined_systematic_budget|sets the combined systematic budget`, dominated by the wiggle/no-wiggle
split, FoG placement, and template-dilation choices catalogued in the decision
register. Fiducial-cosmology dependence is swept through
{astra:decision}`template_cosmology` (a set of AbacusSummit grids), and the
remaining {astra:decision}`fit_range`, {astra:decision}`ells`, and
{astra:decision}`fitting_method` decisions exist for sensitivity tests; none
moves the baseline result beyond its quoted budget.

## 6. Conclusions

A configuration-space BAO analysis of DESI DR1 detects the acoustic feature in
all eight tracer bins post-reconstruction, reaches sub-percent isotropic
precision in the strongest LRG bins, and yields a self-consistent set of
$D_M/r_d$, $D_H/r_d$, $D_V/r_d$ distances at the percent level — reproducing
the DESI DR1 BAO distance ladder. Every numerical claim above is traceable to a
registered ASTRA finding, the decision that configured it, and the materialised
output product it summarises, catalogued below.

---

(analysis-record)=
## Appendix — Analysis record

The article's claims hyperlink into this register of ASTRA objects: the findings
that back each result, the decisions that configured the pipeline, and the
input/output product tables. Prior-literature insights need no visible register —
each inline `{astra:prior-insight}` reference and each decision option carries its
claim and source on hover, and the plugin emits the hidden hover targets those
references resolve against automatically.

### A.1 Findings

:::{astra:finding} bao_peak_sharpens_post_recon
:compact:
:::
:::{astra:finding} bao_detected_post_recon
:compact:
:::
:::{astra:finding} recon_reduces_alpha_iso_sigma
:compact:
:::
:::{astra:finding} subpercent_alpha_iso_precision
:compact:
:::
:::{astra:finding} bao_fit_chi2_near_dof
:compact:
:::

### A.2 Decision register

:::{astra:decision} covariance_source
:::
:::{astra:decision} smoothing_radius
:::
:::{astra:decision} smoothing_radius_qso
:::
:::{astra:decision} broadband
:::
:::{astra:decision} damping_prior
:::
:::{astra:decision} damping_centers
:::
:::{astra:decision} dewiggling_method
:::
:::{astra:decision} fog_placement
:::
:::{astra:decision} dilate_model_component
:::
:::{astra:decision} template_cosmology
:::
:::{astra:decision} systematic_error_treatment
:::
:::{astra:decision} fit_range
:::
:::{astra:decision} ells
:::
:::{astra:decision} fitting_method
:::

(data-products)=
### A.4 Data products

Inputs consumed by the analysis:

:::{astra:inputs}
:::

Outputs registered by the analysis:

:::{astra:outputs}
:::

### A.5 Sub-analyses

:::{astra:subanalysis} reconstruction
:::
:::{astra:subanalysis} clustering
:::
