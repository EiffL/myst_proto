---
title: Clustering
subtitle: Two-point correlation function ξ(s, μ) for DESI DR1 tracers
---

> Sub-analysis page. Components are addressed with the `clustering/` path prefix.

## Overview

Clustering measures the Landy–Szalay ξ(s, μ) multipoles (ℓ = 0, 2, 4) on a
linear 4 Mpc/h s-binning from 0 to 200 Mpc/h, NGC + SGC combined. Pre-recon ξ
is computed from the raw catalogs; post-recon ξ uses the RecSym estimator with
the shifted catalogs from [reconstruction](/reconstruction), sliced into tracer
z-bins.

The binning, μ-binning, estimator, and ℓ range are all locked to the published
RascalC covariance grid, which is why this stage carries only a single decision
rather than a set of binning knobs.

## Decision

The one free choice is whether to apply imaging-systematics weights. The same
RascalC covariance null-tests this toggle at 0.27 σ, so it is left unlocked:

:::{astra} clustering.decisions.imaging_weights
:::

## Outputs

The headline diagnostic is the stacked ξ multipoles figure across tracers:

:::{astra} clustering.outputs.xi_multipoles_plot
:::

The full per-tracer registry — one `xi_pre_recon_<tracer>` and one
`xi_post_recon_<tracer>` per tracer — is:

:::{astra} clustering.outputs
:::
