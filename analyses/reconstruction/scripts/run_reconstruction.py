"""Standard density-field reconstruction for DESI DR1 parent samples (paper §4.2).

Paper §4.2: "The catalogs were reconstructed across the entire redshift range
of each tracer simultaneously" -- so recon runs once per parent sample
(LRG: 0.4-1.1, ELG: 0.8-1.6), not per z-bin. Per-z-bin slicing happens
downstream in the clustering sub-analysis using the saved gal_z / ran_z.

Loads a parent's NGC+SGC clustering catalogs (data + randoms) over the
parent's full z-range, runs reconstruction (algorithm/convention/smoothing/
bias/iterations driven by the universe's decisions), saves the shifted
catalogs. Run via `prism run post_recon_catalog_<parent>`.

Parent-specific settings (catalog prefix, z-range, bias, growth, smoothing)
come from the root `scripts/tracers.py` PARENTS registry.
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from mpi4py import MPI
from scipy.interpolate import interp1d


# Expose the root `scripts/` dir so we can import the parent registry.
_ROOT_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
sys.path.insert(0, str(_ROOT_SCRIPTS))
from tracers import PARENTS, get_parent  # noqa: E402


# ---- DESI box geometry (matches desipipe compute_reconstruction) ------------
CELLSIZE = 4.0              # DESI desipipe hardcoded; nmesh derived by pyrecon
BOX_PAD = 1.2               # 20% FFT padding
# Random file count comes from `parent.nran` in the registry (per-tracer DESI
# values: BGS=1, LRG=8, ELG=10, QSO=4). LSS provides up to 18 random files.

# Decision -> implementation mappings
SMOOTHING = {"sm10": 10.0, "sm15": 15.0, "sm20": 20.0, "sm30": 30.0}
N_ITER = {"n1": 1, "n3": 3, "n5": 5}

CAT_DIR = "/global/cfs/cdirs/desi/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5"
COSMO_FILE = "../../data/desi_fiducial_cosmology.dat"  # relative to sub-analysis cwd


def parse_args():
    p = argparse.ArgumentParser()
    # The script writes TWO outputs (post_recon_catalog_<parent> + mean_displacement_<parent>)
    # in one invocation. {output} from the recipe engine is one of them; we derive the
    # shared universe dir as parent of parent and place both alongside.
    p.add_argument("--output", default=None,
                   help="Destination directory (rendered from {output}); universe dir is its parent.")
    p.add_argument("--universe", default=None,
                   help="Universe name (defaults to parent dir of --output).")
    p.add_argument("--parent", required=True, choices=sorted(PARENTS))
    # Decision flags are optional with defaults matching the baseline universe so that
    # direct CLI invocation (`python run_reconstruction.py --parent X --universe Y`) works.
    p.add_argument("--algorithm", choices=["iterative_fft", "multigrid"],
                   default="iterative_fft")
    p.add_argument("--convention", choices=["recsym", "reciso"], default="recsym")
    p.add_argument("--smoothing_radius", choices=["sm10", "sm15", "sm20"], default="sm15",
                   help="Smoothing for BGS/LRG/ELG parents (paper Table 4: sm15 fiducial).")
    p.add_argument("--smoothing_radius_qso", choices=["sm20", "sm30"], default="sm30",
                   help="Smoothing for QSO parent only (paper Table 4: sm30 fiducial).")
    p.add_argument("--bias_input", choices=["fixed", "estimated_from_data"], default="fixed")
    p.add_argument("--n_iterations", choices=list(N_ITER), default="n3")
    args = p.parse_args()
    if args.output is None and args.universe is None:
        p.error("must pass either --output or --universe")
    return args


def log(msg, t0, rank):
    if rank == 0:
        print(f"[{time.time() - t0:6.1f}s] {msg}", flush=True)


def load_z_to_r():
    """Tabulated z -> r(z) for the DESI fiducial cosmology (matches TabulatedDESI)."""
    z, _, r = np.loadtxt(COSMO_FILE, comments="#", unpack=True)
    return interp1d(z, r, kind="cubic")


def radec_z_to_xyz(ra_deg, dec_deg, z, r_of_z):
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    r = r_of_z(z)
    return np.column_stack([
        r * np.cos(dec) * np.cos(ra),
        r * np.cos(dec) * np.sin(ra),
        r * np.sin(dec),
    ]).astype(np.float64)


def load_catalog(kind, parent, r_of_z, n_files=None):
    if n_files is None:
        n_files = parent.nran
    """Load NGC+SGC catalog for this parent (kind='dat' galaxies or 'ran' randoms).

    Filters to the parent's full z-range so downstream z-bin slicing is cheap
    (gal_z / ran_z are saved alongside positions).
    """
    pos_list, w_list, z_list = [], [], []
    for cap in ["NGC", "SGC"]:
        files = (
            [f"{CAT_DIR}/{parent.catalog_prefix}_{cap}_clustering.dat.fits"]
            if kind == "dat"
            else [f"{CAT_DIR}/{parent.catalog_prefix}_{cap}_{i}_clustering.ran.fits"
                  for i in range(n_files)]
        )
        for fn in files:
            with fits.open(fn) as hdu:
                d = hdu[1].data
                sel = (d["Z"] >= parent.z_min) & (d["Z"] <= parent.z_max)
                pos_list.append(radec_z_to_xyz(d["RA"][sel], d["DEC"][sel], d["Z"][sel], r_of_z))
                w_list.append((d["WEIGHT"][sel] * d["WEIGHT_FKP"][sel]).astype(np.float64))
                z_list.append(d["Z"][sel].astype(np.float64))
    return (np.concatenate(pos_list), np.concatenate(w_list), np.concatenate(z_list))


def desi_growth_rate(z, omega_m=0.3153):
    """f(z) ~= Omega_m(z)^0.55 for DESI fiducial cosmology (Omega_m=0.3153).

    Accurate to ~0.1% vs full class growth_rate for LCDM at z<1 -- using this
    instead of cosmoprimo.fiducial.DESI() to avoid the pyclass dependency.
    """
    e2 = omega_m * (1 + z) ** 3 + (1 - omega_m)
    return (omega_m * (1 + z) ** 3 / e2) ** 0.55


def compute_z_eff(z, weights):
    """Weighted mean redshift, matching desipipe convention."""
    return float(np.sum(z * weights) / np.sum(weights))


def fit_bias_from_pk(gal_pos, ran_pos, gal_w, ran_w, z_eff, mpicomm):
    """Estimate b1 from large-scale FKP P0(k) (k < 0.1 h/Mpc)."""
    from cosmoprimo import Cosmology
    from pypower import CatalogFFTPower

    # Subsample randoms for speed
    if len(ran_pos) > 10_000_000:
        rng = np.random.RandomState(42)
        idx = rng.choice(len(ran_pos), 10_000_000, replace=False)
        ran_pos, ran_w = ran_pos[idx], ran_w[idx]

    kedges = np.arange(0.0, 0.4, 0.005)
    # Size the P(k) box to contain all particles, matching pyrecon's boxpad=1.2
    pk_boxsize = 2 * BOX_PAD * max(np.abs(gal_pos).max(), np.abs(ran_pos).max())
    res = CatalogFFTPower(
        data_positions1=gal_pos, data_weights1=gal_w,
        randoms_positions1=ran_pos, randoms_weights1=ran_w,
        edges=kedges, ells=(0,), los="firstpoint",
        cellsize=CELLSIZE, boxsize=pk_boxsize, boxcenter=0.0,
        resampler="tsc", interlacing=2, position_type="pos", mpicomm=mpicomm,
    )
    k = res.poles.k
    P0 = res.poles(ell=0, complex=False)
    # Eisenstein-Hu transfer — no pyclass dependency; bias fit is insensitive
    # to the small (~%-level) EH-vs-full-TF offset at k < 0.1 h/Mpc.
    cosmo = Cosmology(Omega_cdm=0.2589, Omega_b=0.0486, h=0.6736, n_s=0.9649,
                      sigma8=0.8147, engine="eisenstein_hu")
    Pk_th = cosmo.get_fourier().pk_interpolator(of="delta_cb")(k, z=z_eff)

    valid = (k > 0.02) & (k < 0.1) & np.isfinite(P0) & (Pk_th > 0)
    b_eff_sq = float(np.mean(P0[valid] / Pk_th[valid]))
    f = desi_growth_rate(z_eff)
    # Kaiser RSD: b_eff^2 = b1^2 + (2/3) b1 f + f^2/5  ->  solve for b1
    disc = max(b_eff_sq - 4 * f**2 / 45, 0.01)
    return max(-f / 3 + np.sqrt(disc), 0.1)


def main():
    args = parse_args()
    parent = get_parent(args.parent)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    t0 = time.time()

    if args.output:
        out_dir_arg = Path(args.output)
        universe = args.universe or out_dir_arg.parent.name
        universe_dir = str(out_dir_arg.parent)
    else:
        universe = args.universe
        universe_dir = f"results/{universe}"
    args.universe = universe

    # Paper Table 4: QSO reconstruction uses sm30 while BGS/LRG/ELG use sm15.
    smoothing_tag = args.smoothing_radius_qso if parent.id == "qso_full" else args.smoothing_radius
    smoothing = SMOOTHING[smoothing_tag]
    n_iter = N_ITER[args.n_iterations]

    log(f"Parent: {parent.id}  ({parent.catalog_prefix}, "
        f"z in [{parent.z_min}, {parent.z_max}])", t0, rank)
    log(f"Decisions: algo={args.algorithm} conv={args.convention} "
        f"smooth={smoothing} bias={args.bias_input} niter={n_iter}", t0, rank)

    # ---- Load catalogs (rank 0 only; pyrecon scatters via mpiroot=0) --------
    if rank == 0:
        log("Loading DESI DR1 catalogs...", t0, rank)
        r_of_z = load_z_to_r()
        gal_pos, gal_w, gal_z = load_catalog("dat", parent, r_of_z)
        ran_pos, ran_w, ran_z = load_catalog("ran", parent, r_of_z)
        z_eff = compute_z_eff(gal_z, gal_w)
        f_growth = parent.growth_f
        log(f"  {len(gal_pos):,} galaxies, {len(ran_pos)/1e6:.1f}M randoms "
            f"(ratio={len(ran_pos)/len(gal_pos):.0f}x)", t0, rank)
        log(f"  z_eff = {z_eff:.4f} (weighted), "
            f"f = {f_growth:.4f} (paper Table 4, whole-sample)", t0, rank)
    else:
        gal_pos = gal_w = gal_z = ran_pos = ran_w = ran_z = None
        z_eff = f_growth = None
    z_eff = comm.bcast(z_eff, root=0)
    f_growth = comm.bcast(f_growth, root=0)

    # ---- Bias ---------------------------------------------------------------
    if args.bias_input == "fixed":
        bias = parent.bias
        log(f"Using fixed bias b1 = {bias:.3f}", t0, rank)
    else:
        log("Fitting bias from large-scale P0(k)...", t0, rank)
        bias = fit_bias_from_pk(gal_pos, ran_pos, gal_w, ran_w, z_eff,
                                mpicomm=MPI.COMM_SELF if rank == 0 else None)
        bias = comm.bcast(bias, root=0)
        log(f"Fitted bias b1 = {bias:.3f}", t0, rank)

    # ---- Reconstruction -----------------------------------------------------
    if args.algorithm == "iterative_fft":
        from pyrecon import IterativeFFTReconstruction as Recon
        algo_kw = {"niterations": n_iter}
    else:
        from pyrecon import MultiGridReconstruction as Recon
        algo_kw = {}  # MultiGrid is non-iterative

    log(f"Running {args.algorithm} (cellsize={CELLSIZE}, boxpad={BOX_PAD}, "
        f"los=local)...", t0, rank)
    recon = Recon(
        f=f_growth, bias=bias, los="local",
        data_positions=gal_pos, data_weights=gal_w,
        randoms_positions=ran_pos, randoms_weights=ran_w,
        cellsize=CELLSIZE, boxpad=BOX_PAD,
        wrap=False, smoothing_radius=smoothing,
        fft_plan="estimate", mpiroot=0, mpicomm=comm, dtype="f8",
        **algo_kw,
    )

    # Convention: RecSym shifts both data and randoms by Psi_d = D + f(D.n)n.
    # RecIso shifts data the same way but randoms only by the real-space disp D.
    rand_field = "disp+rsd" if args.convention == "recsym" else "disp"

    log("Reading shifted positions...", t0, rank)
    gal_pos_shift = recon.read_shifted_positions(gal_pos, field="disp+rsd", mpiroot=0)
    ran_pos_shift = recon.read_shifted_positions(ran_pos, field=rand_field, mpiroot=0)

    # ---- Save (rank 0) ------------------------------------------------------
    if rank == 0:
        import json
        mean_disp = float(np.linalg.norm(gal_pos_shift - gal_pos, axis=1).mean())
        log(f"Mean displacement: {mean_disp:.3f} Mpc/h", t0, rank)

        # Convention output paths: results/<universe>/<output_id>/<filename>
        cat_oid = f"post_recon_catalog_{parent.id}"
        metric_oid = f"mean_displacement_{parent.id}"
        cat_dir = f"{universe_dir}/{cat_oid}"
        metric_dir = f"{universe_dir}/{metric_oid}"
        os.makedirs(cat_dir, exist_ok=True)
        os.makedirs(metric_dir, exist_ok=True)

        cat_path = f"{cat_dir}/{cat_oid}.npz"
        np.savez_compressed(
            cat_path,
            gal_pos=gal_pos, gal_pos_shift=gal_pos_shift,
            gal_weights=gal_w, gal_z=gal_z,
            ran_pos=ran_pos, ran_pos_shift=ran_pos_shift,
            ran_weights=ran_w, ran_z=ran_z,
            bias=bias, f_growth=f_growth, mean_displacement=mean_disp,
            cellsize=CELLSIZE, boxpad=BOX_PAD,
            algorithm=args.algorithm, convention=args.convention,
            smoothing_radius=smoothing, n_iterations=n_iter,
            parent=parent.id, catalog_prefix=parent.catalog_prefix, region="GCcomb",
            z_min=parent.z_min, z_max=parent.z_max, z_eff=z_eff,
        )
        log(f"Saved {cat_path}", t0, rank)
        sys.stdout.flush()

        metric_path = f"{metric_dir}/{metric_oid}.json"
        with open(metric_path, "w") as fh:
            json.dump({"value": mean_disp, "units": "Mpc/h"}, fh)
        log(f"Saved {metric_path}", t0, rank)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
