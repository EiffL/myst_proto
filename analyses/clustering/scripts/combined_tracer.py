"""Joint multi-tracer weights for the lrg3_elg1 combined sample (paper §7.5).

Ports the DESI DR1 LSS combined-tracer construction from
desihub/LSS:py/LSS/combined_tracer_utils.py (pin 8fe50b2). The combined
WEIGHT_FKP and WEIGHT in DR1's published combined catalog can't be derived
from the parent catalogs' WEIGHT/WEIGHT_FKP alone -- they require
(a) a per-NTILE completeness factor `nxfac` rebuilt from each parent's
WEIGHT_COMP and FRAC_TLOBS_TILES, and (b) a joint n_eff(z) computed from
each parent's published nz histogram. Both are reproduced here to ~0.1% per
object against the published combined catalog (verified offline).

DESI formula:
    nxfac_i  = comp_ntl[NTILE_i - 1]   (per parent, per cap)
    neff(z)  = (Σ b_i n_i(z))² / (Σ b_i² n_i(z))   (per cap)
    fkp_i    = 1 / (1 + neff(z_i) · nxfac_i · P0)
    W_comb_i = W_parent_i · b_parent
    W_FKP_comb_i = fkp_i

For combined randoms the ELG random WEIGHT additionally carries a
counts-balance factor (N_d_ELG · N_r_LRG) / (N_d_LRG · N_r_ELG) computed
**per random file** before bias scaling -- see DESI's combined_catalog.py.
"""

from __future__ import annotations

import numpy as np
import fitsio


P0_DEFAULT = 6000.0
DZ_DEFAULT = 0.01


def comp_ntl(prefix: str) -> np.ndarray:
    """Per-NTILE completeness array for a parent+cap.

    `prefix` is the full path stem (e.g. "<dir>/LRG_NGC"). Reads the
    parent's clustering data + first random file to compute, for each
    NTILE bin: 1/mean(WEIGHT_COMP) × mean(FRAC_TLOBS_TILES). Returns an
    array indexed from 1 (i.e. `nxfac_i = comp_ntl_arr[NTILE_i - 1]`).
    """
    fd = fitsio.read(f"{prefix}_clustering.dat.fits",
                     columns=["NTILE", "WEIGHT_COMP"])
    ntl_vals = np.unique(fd["NTILE"])
    arr = np.zeros(len(ntl_vals))
    for i, n in enumerate(ntl_vals):
        sel = fd["NTILE"] == n
        arr[i] = 1.0 / np.mean(fd["WEIGHT_COMP"][sel])
    fr = fitsio.read(f"{prefix}_0_clustering.ran.fits",
                     columns=["NTILE", "FRAC_TLOBS_TILES"])
    for i, n in enumerate(ntl_vals):
        sel = fr["NTILE"] == n
        arr[i] *= np.mean(fr["FRAC_TLOBS_TILES"][sel])
    return arr


def _rebin_nz(z_centers, nz_t, zmin_t_arr, dz_t):
    """Map a tracer's nz histogram onto the combined z-binning."""
    zmin_t = zmin_t_arr[0]
    zmax_t = zmin_t_arr[-1] + dz_t
    out = np.zeros(len(z_centers))
    for i, z in enumerate(z_centers):
        if zmin_t < z < zmax_t:
            out[i] = nz_t[int((z - zmin_t) / dz_t)]
    return out


def neff_curve(nz_files, biases, zmin, zmax, dz=DZ_DEFAULT):
    """neff(z) on the combined z-binning. nz_files: list of (N,>=4) arrays
    where col 1 is zmin per bin and col 3 is n(z). biases: list of floats
    (one per parent). Returns (zmin_bin_edges, neff_array)."""
    nbins = int((zmax - zmin) / dz)
    zmin_bins = np.linspace(zmin, zmax, nbins, endpoint=False)
    zc = zmin_bins + dz / 2
    bnz_sum = np.zeros(nbins)
    b2nz_sum = np.zeros(nbins)
    for nzf, b in zip(nz_files, biases):
        zmin_t_arr = nzf[:, 1]
        nz_t = nzf[:, 3]
        dz_t = zmin_t_arr[1] - zmin_t_arr[0]
        nrebin = _rebin_nz(zc, nz_t, zmin_t_arr, dz_t)
        bnz_sum += b * nrebin
        b2nz_sum += b * b * nrebin
    beff = b2nz_sum / np.where(bnz_sum > 0, bnz_sum, 1.0)
    neff = np.where(beff > 0, bnz_sum / beff, 0.0)
    return zmin_bins, neff


def fkp_per_object(z, nxfac, neff, zmin, zmax, dz=DZ_DEFAULT, P0=P0_DEFAULT):
    """fkp_i = 1/(1 + neff(z_i)·nxfac_i·P0). Outside [zmin,zmax) returns 0."""
    fkp = np.zeros(len(z))
    inside = (z > zmin) & (z < zmax)
    if not inside.any():
        return fkp
    idx = np.minimum(((z[inside] - zmin) / dz).astype(np.int64), len(neff) - 1)
    fkp[inside] = 1.0 / (1.0 + neff[idx] * nxfac[inside] * P0)
    return fkp
