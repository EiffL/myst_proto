"""Two-point correlation xi(s, mu) for DESI DR1 tracers (paper §4.1).

Computes pre- and post-reconstruction xi via pycorr's Landy-Szalay estimator,
NGC/SGC separately, then combined (GCcomb). Pre-recon: data x randoms from raw
catalogs filtered to the tracer z-range. Post-recon: RecSym (shifted data +
shifted randoms + unshifted randoms) from the parent-sample reconstruction
NPZ, sliced to the tracer z-range via the stored gal_z / ran_z arrays.

Composite tracers (lrg3_elg1, paper §7.5) are handled via the
`combined_tracer` helper, which ports DESI's joint FKP construction
(per-NTILE nxfac × neff(z) joint FKP, plus bias-multiplied WEIGHT) so we
build the combined weights ourselves rather than cloning the VAC.

**Per-random-file split:** pycorr's TwoPointCorrelationFunction objects sum
their pair counts under `+`, so we compute one TPCF per random file and sum
the results rather than concatenating all randoms into one giant array. This
keeps RR cost ~linear in N_RANDOM_FILES instead of quadratic, which is the
trick David Valcin pointed out -- bumping randoms from 4 to 10 actually gets
*faster* this way while pulling closer to the noise-free RR limit.

Settings (fixed to match desipipe / DESI RascalC covariance):
  mode='smu', s-edges=np.arange(0, 201, 4), mu-edges=linspace(-1, 1, 201),
  los='firstpoint', estimator='landyszalay' (pycorr default).

Tracer-specific configuration (parents, z-range) comes from the root
`scripts/tracers.py` registry. Writes (one directory per declared output_id):
  results/<universe>/xi_pre_recon_<tracer>/xi_pre_recon_<tracer>.npy
  results/<universe>/xi_post_recon_<tracer>/xi_post_recon_<tracer>.npy
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from mpi4py import MPI
from pycorr import TwoPointCorrelationFunction
from scipy.interpolate import interp1d


_ROOT_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
sys.path.insert(0, str(_ROOT_SCRIPTS))
from tracers import TRACERS, PARENTS, get as get_tracer, get_parent  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import combined_tracer as ct  # noqa: E402


# Random file counts come from the registry: tracer.nran (pre-recon) and
# parent.nran (recon). For post-recon we additionally clamp to whatever the
# existing NPZ was actually built with — see _detect_recon_nran below.
SEDGES = np.arange(0, 201, 4)
MUEDGES = np.linspace(-1, 1, 201)

CAT_DIR = "/global/cfs/cdirs/desi/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5"
COSMO_FILE = "../../data/desi_fiducial_cosmology.dat"
RECON_NPZ = "../reconstruction/results/{universe}/post_recon_catalog_{parent}.npz"


def parse_args():
    p = argparse.ArgumentParser()
    # The recipe engine renders {output} as a single output dir
    # (e.g. results/<universe>/xi_pre_recon_<tracer>/). compute_xi writes BOTH
    # pre and post xi for the tracer in one invocation, so we derive the
    # shared parent universe dir from --output (parents[1]) and place pre/post
    # under sibling output_id subdirs there. Keeps `lc run xi_pre_recon_X` and
    # `lc run xi_post_recon_X` interchangeable while honouring the engine's
    # convention that the script writes into {output}.
    p.add_argument("--output", default=None,
                   help="Destination directory (rendered from {output}); the universe "
                        "dir is derived from it (parent of parent).")
    p.add_argument("--universe", default=None,
                   help="Universe name (defaults to parent dir of --output).")
    p.add_argument("--tracer", required=True, choices=sorted(TRACERS))
    p.add_argument("--imaging_weights", choices=["on", "off"], default="on",
                   help="`on`: WEIGHT * WEIGHT_FKP (DESI fiducial). "
                        "`off`: (WEIGHT/WEIGHT_SYS) * WEIGHT_FKP -- drops "
                        "imaging-systematics correction (paper §5.3 null test).")
    args = p.parse_args()
    if args.output is None and args.universe is None:
        p.error("must pass either --output or --universe")
    return args


def _imaging_W(d, imaging_weights):
    """WEIGHT (or WEIGHT/WEIGHT_SYS), no FKP."""
    w = d["WEIGHT"].astype(np.float64)
    if imaging_weights == "off":
        w = w / d["WEIGHT_SYS"].astype(np.float64)
    return w


def _w_full(d, imaging_weights):
    """Per-DESI parent weight: imaging × WEIGHT_FKP. Used for single-tracer."""
    return _imaging_W(d, imaging_weights) * d["WEIGHT_FKP"].astype(np.float64)


def log(msg, t0, rank):
    if rank == 0:
        print(f"[{time.time() - t0:6.1f}s] {msg}", flush=True)


def load_z_to_r():
    z, _, r = np.loadtxt(COSMO_FILE, comments="#", unpack=True)
    return interp1d(z, r, kind="cubic")


def radec_z_to_xyz(ra_deg, dec_deg, z, r_of_z):
    ra = np.deg2rad(ra_deg); dec = np.deg2rad(dec_deg)
    r = r_of_z(z)
    return np.column_stack([
        r * np.cos(dec) * np.cos(ra),
        r * np.cos(dec) * np.sin(ra),
        r * np.sin(dec),
    ]).astype(np.float64)


def is_ngc(pos):
    ra = np.degrees(np.arctan2(pos[:, 1], pos[:, 0])) % 360
    return (ra > 80) & (ra < 300)


# ---------------------------------------------------------------------------
# Single-tracer pre-recon loaders
# ---------------------------------------------------------------------------

def load_data_by_cap(tracer, r_of_z, imaging_weights):
    """Pre-recon data per cap: {cap: (pos, w)}."""
    parent = get_parent(tracer.parent[0])
    out = {}
    for cap in ("NGC", "SGC"):
        with fits.open(f"{CAT_DIR}/{parent.catalog_prefix}_{cap}_clustering.dat.fits") as hdu:
            d = hdu[1].data
            sel = (d["Z"] >= tracer.z_min) & (d["Z"] <= tracer.z_max)
            out[cap] = (radec_z_to_xyz(d["RA"][sel], d["DEC"][sel], d["Z"][sel], r_of_z),
                        _w_full(d[sel], imaging_weights))
    return out


def load_random_by_cap(tracer, fi, r_of_z, imaging_weights):
    """Pre-recon randoms for one random file index `fi`: {cap: (pos, w)}."""
    parent = get_parent(tracer.parent[0])
    out = {}
    for cap in ("NGC", "SGC"):
        with fits.open(f"{CAT_DIR}/{parent.catalog_prefix}_{cap}_{fi}_clustering.ran.fits") as hdu:
            d = hdu[1].data
            sel = (d["Z"] >= tracer.z_min) & (d["Z"] <= tracer.z_max)
            out[cap] = (radec_z_to_xyz(d["RA"][sel], d["DEC"][sel], d["Z"][sel], r_of_z),
                        _w_full(d[sel], imaging_weights))
    return out


# ---------------------------------------------------------------------------
# Combined-tracer pre-recon loaders (joint FKP, bias-multiplied WEIGHT;
# ELG randoms get per-(file,cap) counts-balance factor).
# ---------------------------------------------------------------------------

def _combined_neff_per_cap(tracer, cap):
    parents = [get_parent(p) for p in tracer.parent]
    biases = [p.bias for p in parents]
    nz_files = [np.loadtxt(f"{CAT_DIR}/{p.catalog_prefix}_{cap}_nz.txt")
                for p in parents]
    _, neff = ct.neff_curve(nz_files, biases, tracer.z_min, tracer.z_max)
    return neff


def _combined_fkp(d, sel, neff, comp_arr, tracer):
    z = np.asarray(d["Z"][sel]).astype(np.float64)
    ntile = np.asarray(d["NTILE"][sel]).astype(np.int64)
    return ct.fkp_per_object(z, comp_arr[ntile - 1], neff,
                             tracer.z_min, tracer.z_max)


def load_combined_data(tracer, r_of_z, imaging_weights):
    """Pre-recon combined data. Returns ({cap: (pos, w)}, n_d) where n_d is
    {(parent_id, cap): sum(W·fkp)} -- needed for the random counts-balance factor."""
    parents = [get_parent(pid) for pid in tracer.parent]
    out = {}
    n_d = {}
    for cap in ("NGC", "SGC"):
        neff = _combined_neff_per_cap(tracer, cap)
        positions, weights = [], []
        for p in parents:
            comp_arr = ct.comp_ntl(f"{CAT_DIR}/{p.catalog_prefix}_{cap}")
            with fits.open(f"{CAT_DIR}/{p.catalog_prefix}_{cap}_clustering.dat.fits") as hdu:
                d = hdu[1].data
                sel = (d["Z"] >= tracer.z_min) & (d["Z"] <= tracer.z_max)
                fkp = _combined_fkp(d, sel, neff, comp_arr, tracer)
                W = _imaging_W(d[sel], imaging_weights)
                positions.append(radec_z_to_xyz(d["RA"][sel], d["DEC"][sel], d["Z"][sel], r_of_z))
                weights.append(W * fkp * p.bias)
                n_d[(p.id, cap)] = float(np.sum(W * fkp))
        out[cap] = (np.concatenate(positions), np.concatenate(weights))
    return out, n_d


def load_combined_random(tracer, fi, r_of_z, imaging_weights, n_d):
    """Pre-recon combined random for one file index. {cap: (pos, w)}."""
    parents = [get_parent(pid) for pid in tracer.parent]
    out = {}
    for cap in ("NGC", "SGC"):
        neff = _combined_neff_per_cap(tracer, cap)
        loaded = []  # (parent, pos, w_with_bias, n_r)
        for p in parents:
            comp_arr = ct.comp_ntl(f"{CAT_DIR}/{p.catalog_prefix}_{cap}")
            with fits.open(f"{CAT_DIR}/{p.catalog_prefix}_{cap}_{fi}_clustering.ran.fits") as hdu:
                d = hdu[1].data
                sel = (d["Z"] >= tracer.z_min) & (d["Z"] <= tracer.z_max)
                fkp = _combined_fkp(d, sel, neff, comp_arr, tracer)
                W = _imaging_W(d[sel], imaging_weights)
                pos = radec_z_to_xyz(d["RA"][sel], d["DEC"][sel], d["Z"][sel], r_of_z)
                loaded.append((p, pos, W * fkp * p.bias, float(np.sum(W * fkp))))
        # Apply (N_d_E·N_r_L)/(N_d_L·N_r_E) factor to non-first-parent's randoms.
        base_p, _, _, base_nr = loaded[0]
        positions, weights = [], []
        for pi, (p, pos, w, nr) in enumerate(loaded):
            if pi >= 1:
                w = w * (n_d[(p.id, cap)] * base_nr / (n_d[(base_p.id, cap)] * nr))
            positions.append(pos)
            weights.append(w)
        out[cap] = (np.concatenate(positions), np.concatenate(weights))
    return out


# ---------------------------------------------------------------------------
# Post-recon loaders. NPZ rows are pre-aligned with parent FITS in iteration
# order [NGC dat | SGC dat] for galaxies and [NGC_0..N-1 | SGC_0..N-1] for
# randoms. We re-read the parent FITS once to recover NTILE (combined-tracer)
# and per-row file index (for per-file random splitting).
# ---------------------------------------------------------------------------

def _detect_recon_nran(npz, parent):
    """Infer how many random files went into the recon NPZ by counting how
    many parent random files match its total ran row count. Lets us read
    NPZs built with old N=4 even though parent.nran might now be 8 or 10."""
    target = len(npz["ran_z"])
    cumul = 0
    for fi in range(20):  # safety cap; LSS has up to 18
        try:
            with fits.open(f"{CAT_DIR}/{parent.catalog_prefix}_NGC_{fi}_clustering.ran.fits") as hdu:
                d = hdu[1].data
                sel = (d["Z"] >= parent.z_min) & (d["Z"] <= parent.z_max)
                cumul += int(sel.sum())
            with fits.open(f"{CAT_DIR}/{parent.catalog_prefix}_SGC_{fi}_clustering.ran.fits") as hdu:
                d = hdu[1].data
                sel = (d["Z"] >= parent.z_min) & (d["Z"] <= parent.z_max)
                cumul += int(sel.sum())
        except FileNotFoundError:
            break
        if cumul == target:
            return fi + 1
    raise RuntimeError(
        f"Could not infer recon nran for {parent.id}: NPZ has {target} ran rows, "
        f"reading {fi+1} random files yields {cumul}."
    )


def _reload_parent(parent, kind, cols, n_random):
    """Concatenated columns + per-row file index, aligned to NPZ rows. Uses
    exactly `n_random` random files per cap."""
    out = {c: [] for c in cols}
    fidx = []
    for cap in ("NGC", "SGC"):
        files = (
            [(f"{CAT_DIR}/{parent.catalog_prefix}_{cap}_clustering.dat.fits", -1)]
            if kind == "dat"
            else [(f"{CAT_DIR}/{parent.catalog_prefix}_{cap}_{i}_clustering.ran.fits", i)
                  for i in range(n_random)]
        )
        for fn, fi in files:
            with fits.open(fn) as hdu:
                d = hdu[1].data
                sel = (d["Z"] >= parent.z_min) & (d["Z"] <= parent.z_max)
                for c in cols:
                    out[c].append(np.asarray(d[c][sel]))
                fidx.append(np.full(int(sel.sum()), fi, dtype=np.int64))
    return ({c: np.concatenate(v) for c, v in out.items()},
            np.concatenate(fidx))


class PostReconSingle:
    """Pre-loads the NPZ + parent FITS columns once, exposes per-(cap, fi)
    random accessors so xi can be computed file-by-file."""
    def __init__(self, universe, tracer, imaging_weights):
        self.tracer = tracer
        self.iw = imaging_weights
        self.parent = get_parent(tracer.parent[0])
        npz = np.load(RECON_NPZ.format(universe=universe, parent=self.parent.id))
        self.npz = npz
        self.nran = _detect_recon_nran(npz, self.parent)
        # Galaxy weights: NPZ has them with imaging_weights="on"; reload for "off".
        if imaging_weights == "on":
            self.gw_full = npz["gal_weights"]
            self.rw_full = npz["ran_weights"]
        else:
            gcols, _ = _reload_parent(self.parent, "dat",
                                      ["WEIGHT", "WEIGHT_SYS", "WEIGHT_FKP"], self.nran)
            assert len(gcols["WEIGHT"]) == len(npz["gal_z"])
            self.gw_full = (gcols["WEIGHT"]/gcols["WEIGHT_SYS"]) * gcols["WEIGHT_FKP"]
            rcols, _ = _reload_parent(self.parent, "ran",
                                      ["WEIGHT", "WEIGHT_SYS", "WEIGHT_FKP"], self.nran)
            assert len(rcols["WEIGHT"]) == len(npz["ran_z"])
            self.rw_full = (rcols["WEIGHT"]/rcols["WEIGHT_SYS"]) * rcols["WEIGHT_FKP"]
        # Per-NPZ-row random file index (0..nran-1).
        _, self.r_fidx = _reload_parent(self.parent, "ran", ["Z"], self.nran)
        # Cap masks + tracer-z mask.
        self.is_g_ngc = is_ngc(npz["gal_pos"])
        self.is_r_ngc = is_ngc(npz["ran_pos"])
        self.gz_mask = (npz["gal_z"] >= tracer.z_min) & (npz["gal_z"] <= tracer.z_max)
        self.rz_mask = (npz["ran_z"] >= tracer.z_min) & (npz["ran_z"] <= tracer.z_max)

    def galaxies(self):
        out = {}
        for cap, ngc in [("NGC", self.is_g_ngc), ("SGC", ~self.is_g_ngc)]:
            m = self.gz_mask & ngc
            out[cap] = (self.npz["gal_pos_shift"][m], self.gw_full[m])
        return out

    def random(self, fi):
        """Returns {cap: (rs, rw, ru)} for random file `fi`."""
        m_file = (self.r_fidx == fi)
        out = {}
        for cap, ngc in [("NGC", self.is_r_ngc), ("SGC", ~self.is_r_ngc)]:
            m = self.rz_mask & ngc & m_file
            out[cap] = (self.npz["ran_pos_shift"][m],
                        self.rw_full[m],
                        self.npz["ran_pos"][m])
        return out


class PostReconCombined:
    """Combined-tracer post-recon. For each parent: NPZ + reload FITS
    columns (NTILE/Z/WEIGHT/WEIGHT_SYS) once, plus precomputed neff and
    per-NTILE comp_ntl per cap. Per-(cap, fi) accessor combines the two
    parents with bias-multiplied weights and joint FKP, plus the
    counts-balance factor on non-first-parent's randoms."""
    def __init__(self, universe, tracer, imaging_weights):
        self.tracer = tracer
        self.iw = imaging_weights
        self.parents = [get_parent(pid) for pid in tracer.parent]
        self.neff = {cap: _combined_neff_per_cap(tracer, cap) for cap in ("NGC", "SGC")}
        self.comp = {(p.id, cap): ct.comp_ntl(f"{CAT_DIR}/{p.catalog_prefix}_{cap}")
                     for p in self.parents for cap in ("NGC", "SGC")}
        self.npz = {p.id: np.load(RECON_NPZ.format(universe=universe, parent=p.id))
                    for p in self.parents}
        self.nran = {p.id: _detect_recon_nran(self.npz[p.id], p) for p in self.parents}
        self.gcols = {}
        self.rcols = {}
        self.r_fidx = {}
        for p in self.parents:
            gc, _ = _reload_parent(p, "dat", ["Z", "NTILE", "WEIGHT", "WEIGHT_SYS"],
                                    self.nran[p.id])
            assert len(gc["Z"]) == len(self.npz[p.id]["gal_z"])
            self.gcols[p.id] = gc
            rc, fidx = _reload_parent(p, "ran", ["Z", "NTILE", "WEIGHT", "WEIGHT_SYS"],
                                       self.nran[p.id])
            assert len(rc["Z"]) == len(self.npz[p.id]["ran_z"])
            self.rcols[p.id] = rc
            self.r_fidx[p.id] = fidx
        # Per-parent-per-cap N_d (over data, in tracer z-range) and per-(file, cap) N_r.
        self.n_d = {}            # (pid, cap) -> N_d
        self.n_r = {}            # (pid, cap, fi) -> N_r
        for p in self.parents:
            npz = self.npz[p.id]
            gc = self.gcols[p.id]
            tz_g = (npz["gal_z"] >= tracer.z_min) & (npz["gal_z"] <= tracer.z_max)
            is_g_ngc = is_ngc(npz["gal_pos"])
            for cap, ngc in [("NGC", is_g_ngc), ("SGC", ~is_g_ngc)]:
                m = tz_g & ngc
                ntile = gc["NTILE"][m].astype(np.int64)
                z_arr = gc["Z"][m].astype(np.float64)
                W_arr = gc["WEIGHT"][m].astype(np.float64)
                if imaging_weights == "off":
                    W_arr = W_arr / gc["WEIGHT_SYS"][m].astype(np.float64)
                comp_arr = self.comp[(p.id, cap)]
                fkp = ct.fkp_per_object(z_arr, comp_arr[ntile - 1],
                                        self.neff[cap], tracer.z_min, tracer.z_max)
                self.n_d[(p.id, cap)] = float(np.sum(W_arr * fkp))
            rc = self.rcols[p.id]
            fidx = self.r_fidx[p.id]
            tz_r = (npz["ran_z"] >= tracer.z_min) & (npz["ran_z"] <= tracer.z_max)
            is_r_ngc = is_ngc(npz["ran_pos"])
            for fi in range(self.nran[p.id]):
                for cap, ngc in [("NGC", is_r_ngc), ("SGC", ~is_r_ngc)]:
                    m = tz_r & ngc & (fidx == fi)
                    if not m.any():
                        self.n_r[(p.id, cap, fi)] = 0.0
                        continue
                    ntile = rc["NTILE"][m].astype(np.int64)
                    z_arr = rc["Z"][m].astype(np.float64)
                    W_arr = rc["WEIGHT"][m].astype(np.float64)
                    if imaging_weights == "off":
                        W_arr = W_arr / rc["WEIGHT_SYS"][m].astype(np.float64)
                    comp_arr = self.comp[(p.id, cap)]
                    fkp = ct.fkp_per_object(z_arr, comp_arr[ntile - 1],
                                            self.neff[cap], tracer.z_min, tracer.z_max)
                    self.n_r[(p.id, cap, fi)] = float(np.sum(W_arr * fkp))

    def _combined_w_for(self, p, mask, cap, kind):
        cols = self.gcols[p.id] if kind == "dat" else self.rcols[p.id]
        ntile = cols["NTILE"][mask].astype(np.int64)
        z_arr = cols["Z"][mask].astype(np.float64)
        W_arr = cols["WEIGHT"][mask].astype(np.float64)
        if self.iw == "off":
            W_arr = W_arr / cols["WEIGHT_SYS"][mask].astype(np.float64)
        comp_arr = self.comp[(p.id, cap)]
        fkp = ct.fkp_per_object(z_arr, comp_arr[ntile - 1],
                                self.neff[cap], self.tracer.z_min, self.tracer.z_max)
        return W_arr * fkp * p.bias

    def galaxies(self):
        out = {"NGC": [[], []], "SGC": [[], []]}
        for p in self.parents:
            npz = self.npz[p.id]
            tz = (npz["gal_z"] >= self.tracer.z_min) & (npz["gal_z"] <= self.tracer.z_max)
            is_g_ngc = is_ngc(npz["gal_pos"])
            for cap, ngc in [("NGC", is_g_ngc), ("SGC", ~is_g_ngc)]:
                m = tz & ngc
                out[cap][0].append(npz["gal_pos_shift"][m])
                out[cap][1].append(self._combined_w_for(p, m, cap, "dat"))
        return {cap: (np.concatenate(p), np.concatenate(w))
                for cap, (p, w) in out.items()}

    def random(self, fi):
        out = {"NGC": [[], [], []], "SGC": [[], [], []]}
        base = self.parents[0]
        for pi, p in enumerate(self.parents):
            npz = self.npz[p.id]
            fidx = self.r_fidx[p.id]
            tz = (npz["ran_z"] >= self.tracer.z_min) & (npz["ran_z"] <= self.tracer.z_max)
            is_r_ngc = is_ngc(npz["ran_pos"])
            for cap, ngc in [("NGC", is_r_ngc), ("SGC", ~is_r_ngc)]:
                m = tz & ngc & (fidx == fi)
                w = self._combined_w_for(p, m, cap, "ran")
                if pi >= 1:
                    nr_p  = self.n_r[(p.id, cap, fi)]
                    nr_b  = self.n_r[(base.id, cap, fi)]
                    nd_p  = self.n_d[(p.id, cap)]
                    nd_b  = self.n_d[(base.id, cap)]
                    if nr_p > 0 and nd_b > 0:
                        w = w * (nd_p * nr_b / (nd_b * nr_p))
                out[cap][0].append(npz["ran_pos_shift"][m])
                out[cap][1].append(w)
                out[cap][2].append(npz["ran_pos"][m])
        return {cap: (np.concatenate(rs), np.concatenate(rw), np.concatenate(ru))
                for cap, (rs, rw, ru) in out.items()}


# ---------------------------------------------------------------------------
# pycorr wrappers
# ---------------------------------------------------------------------------

def _xi_pre(gp, gw, rp, rw, comm, nthreads):
    return TwoPointCorrelationFunction(
        mode="smu", edges=(SEDGES, MUEDGES),
        data_positions1=gp, data_weights1=gw,
        randoms_positions1=rp, randoms_weights1=rw,
        los="firstpoint", position_type="pos",
        mpicomm=comm, mpiroot=0, nthreads=nthreads,
    )


def _xi_post(gp, gw, rs, rw, ru, comm, nthreads):
    return TwoPointCorrelationFunction(
        mode="smu", edges=(SEDGES, MUEDGES),
        data_positions1=gp, data_weights1=gw,
        shifted_positions1=rs, shifted_weights1=rw,
        randoms_positions1=ru, randoms_weights1=rw,
        los="firstpoint", position_type="pos",
        mpicomm=comm, mpiroot=0, nthreads=nthreads,
    )


def main():
    args = parse_args()
    tracer = get_tracer(args.tracer)
    is_combined = len(tracer.parent) > 1
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nthreads = int(os.environ.get("OMP_NUM_THREADS", 32))
    t0 = time.time()

    if args.output:
        # {output} → results/<universe>/<output_id>/. Universe dir = parents[1].
        out_dir_arg = Path(args.output)
        universe = args.universe or out_dir_arg.parent.name
        universe_dir = str(out_dir_arg.parent)
    else:
        universe = args.universe
        universe_dir = f"results/{universe}"
    args.universe = universe  # for downstream RECON_NPZ.format(...) call
    pre_id = f"xi_pre_recon_{tracer.id}"
    post_id = f"xi_post_recon_{tracer.id}"
    pre_dir = f"{universe_dir}/{pre_id}"
    post_dir = f"{universe_dir}/{post_id}"
    if rank == 0:
        os.makedirs(pre_dir, exist_ok=True)
        os.makedirs(post_dir, exist_ok=True)
        log(f"Tracer: {tracer.id}  (parents={list(tracer.parent)}, "
            f"z in [{tracer.z_min}, {tracer.z_max}], pre-recon nran={tracer.nran})",
            t0, rank)

    # ------------------------------- Pre-recon -------------------------------
    if rank == 0:
        log("Loading pre-recon data...", t0, rank)
        r_of_z = load_z_to_r()
        if is_combined:
            gal_caps, n_d = load_combined_data(tracer, r_of_z, args.imaging_weights)
        else:
            gal_caps = load_data_by_cap(tracer, r_of_z, args.imaging_weights)
            n_d = None
        for cap in ("NGC", "SGC"):
            log(f"  data {cap}: {len(gal_caps[cap][0]):,} gal", t0, rank)
    else:
        gal_caps = None; n_d = None; r_of_z = None

    cap_results = []
    for cap in ("NGC", "SGC"):
        xi_cap = None
        for fi in range(tracer.nran):
            if rank == 0:
                if is_combined:
                    ran_caps = load_combined_random(tracer, fi, r_of_z, args.imaging_weights, n_d)
                else:
                    ran_caps = load_random_by_cap(tracer, fi, r_of_z, args.imaging_weights)
                gp, gw = gal_caps[cap]
                rp, rw = ran_caps[cap]
                if fi == 0:
                    log(f"  {cap} ran file 0: {len(rp):,}", t0, rank)
            else:
                gp = gw = rp = rw = None
            log(f"Pre-recon xi: {cap} ran-file {fi+1}/{tracer.nran}...", t0, rank)
            xi_i = _xi_pre(gp, gw, rp, rw, comm, nthreads)
            xi_cap = xi_i if xi_cap is None else xi_cap + xi_i
        cap_results.append(xi_cap)

    if rank == 0:
        xi_gccomb = sum(cap_results)
        pre_path = f"{pre_dir}/{pre_id}.npy"
        xi_gccomb.save(pre_path)
        log(f"Saved {pre_path}", t0, rank)
        del gal_caps, cap_results, xi_gccomb
    comm.Barrier()

    # ------------------------------- Post-recon ------------------------------
    if rank == 0:
        log("Loading post-recon catalogs...", t0, rank)
        if is_combined:
            ctx = PostReconCombined(args.universe, tracer, args.imaging_weights)
            post_nran = min(ctx.nran.values())     # combined post-recon limited by min over parents
        else:
            ctx = PostReconSingle(args.universe, tracer, args.imaging_weights)
            post_nran = ctx.nran
        log(f"  post-recon nran={post_nran} (from existing NPZ)", t0, rank)
        gal_caps = ctx.galaxies()
        for cap in ("NGC", "SGC"):
            log(f"  gal {cap}: {len(gal_caps[cap][0]):,}", t0, rank)
    else:
        ctx = None; gal_caps = None; post_nran = None
    post_nran = comm.bcast(post_nran, root=0)

    cap_results = []
    for cap in ("NGC", "SGC"):
        xi_cap = None
        for fi in range(post_nran):
            if rank == 0:
                ran_caps = ctx.random(fi)
                gp, gw = gal_caps[cap]
                rs, rw, ru = ran_caps[cap]
                if fi == 0:
                    log(f"  {cap} ran file 0: {len(rs):,}", t0, rank)
            else:
                gp = gw = rs = rw = ru = None
            log(f"Post-recon xi: {cap} ran-file {fi+1}/{post_nran}...", t0, rank)
            xi_i = _xi_post(gp, gw, rs, rw, ru, comm, nthreads)
            xi_cap = xi_i if xi_cap is None else xi_cap + xi_i
        cap_results.append(xi_cap)

    if rank == 0:
        xi_gccomb = sum(cap_results)
        post_path = f"{post_dir}/{post_id}.npy"
        xi_gccomb.save(post_path)
        log(f"Saved {post_path}", t0, rank)

    comm.Barrier()


if __name__ == "__main__":
    main()
