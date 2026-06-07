---
title: Reconstruction
subtitle: Density-field reconstruction for DESI DR1 parent samples
---

> Sub-analysis page. The `astra:*` components below are scoped into the
> `reconstruction` sub-analysis with the `reconstruction.` path prefix.

## Overview

Density-field reconstruction runs first in the pipeline. The linear
displacement field Ψ is estimated from the Gaussian-smoothed raw density and
applied symmetrically to galaxies and randoms (RecSym), sharpening the BAO peak
before clustering measures ξ(s). One pass per parent sample (BGS, LRG, ELG, QSO)
emits parent-level shifted catalogs; per-z-bin slicing is deferred to the
[clustering](/clustering) stage.

## Decisions

The reconstruction algorithm is IterativeFFT (fiducial), with MultiGrid exposed
as a sensitivity alternative:

:::{astra:decision} reconstruction.algorithm
:::

The displacement convention is RecSym — both data and randoms are shifted, which
is what keeps the post-reconstruction ξ estimator unbiased:

:::{astra:decision} reconstruction.convention
:::

The Gaussian smoothing scale is not decided here: it is inherited from the root
via `from: ../smoothing_radius`, so the Σ_sm used in the displacement solve and
the Σ_sm in the downstream BAO damping template stay identical by construction.

## Outputs

Reconstruction emits a shifted galaxy + random catalog per parent sample,
alongside diagnostic mean-displacement metrics:

:::{astra:outputs} reconstruction
:::
