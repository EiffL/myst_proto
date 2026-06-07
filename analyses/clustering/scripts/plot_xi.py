"""Plot pre- and post-recon xi_ell(s) multipoles as s^2 * xi(s), per tracer.

Loads the xi outputs for every tracer in the registry that has both
xi_pre_recon_<tracer>.npy and xi_post_recon_<tracer>.npy present, plots
ell=0,2,4 for both pre- and post-recon over s in [30, 180] Mpc/h, one row
per tracer.
"""

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from pycorr import TwoPointCorrelationFunction


# Expose the root `scripts/` dir so we can import the tracer registry.
_ROOT_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
sys.path.insert(0, str(_ROOT_SCRIPTS))
from bao_utils import TRACER_COLORS  # noqa: E402
from tracers import TRACERS, get as get_tracer, rascalc_path  # noqa: E402


ELLS = (0, 2, 4)
S_LIM = (30, 180)                    # focus on the BAO-relevant range
COV_S_EDGES = np.linspace(20.0, 200.0, 46)
COV_ELLS = (0, 2, 4)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--output", default=None,
                   help="Destination directory (rendered from {output} by the recipe engine).")
    p.add_argument("--universe", default=None,
                   help="Universe name (defaults to parent dir of --output).")
    args = p.parse_args()
    if args.output is None and args.universe is None:
        p.error("must pass either --output or --universe")
    return args


def load_xi(path):
    r = TwoPointCorrelationFunction.load(path)
    s, xi = r(ells=ELLS, return_sep=True, return_std=False)
    return s, dict(zip(ELLS, xi))


def load_cov_sigmas(cov_path):
    """Diagonal 1-sigma uncertainties per ell on the cov's native s-grid."""
    cov = np.genfromtxt(cov_path)
    nbins = len(COV_S_EDGES) - 1
    assert cov.shape == (nbins * len(COV_ELLS),) * 2, f"unexpected cov shape {cov.shape}"
    s_mid = 0.5 * (COV_S_EDGES[1:] + COV_S_EDGES[:-1])
    sigma = {
        ell: np.sqrt(np.diag(cov)[i * nbins:(i + 1) * nbins])
        for i, ell in enumerate(COV_ELLS)
    }
    return s_mid, sigma


def main():
    args = parse_args()
    if args.output:
        cli_out_dir = Path(args.output)
        universe = args.universe or cli_out_dir.parent.name
    else:
        universe = args.universe
        cli_out_dir = None
    args.universe = universe
    out_dir = f"results/{universe}"

    tracers = []
    for tracer_id in sorted(TRACERS):
        tracer = get_tracer(tracer_id)
        pre_id = f"xi_pre_recon_{tracer.id}"
        post_id = f"xi_post_recon_{tracer.id}"
        pre_path = f"{out_dir}/{pre_id}/{pre_id}.npy"
        post_path = f"{out_dir}/{post_id}/{post_id}.npy"
        if not (os.path.exists(pre_path) and os.path.exists(post_path)):
            print(f"  skipped {tracer.id}: missing xi output")
            continue
        tracers.append((tracer, pre_path, post_path))

    if not tracers:
        raise SystemExit(f"No xi outputs found under {out_dir}.")

    nrows = len(tracers)
    fig, axes = plt.subplots(nrows, 3, figsize=(15, 4.5 * nrows),
                             sharex=True, squeeze=False)

    for row, (tracer, pre_path, post_path) in enumerate(tracers):
        s_pre, xi_pre = load_xi(pre_path)
        s_post, xi_post = load_xi(post_path)
        # Post-recon cov gives a useful error band around the post-recon curve
        s_cov, sigma_cov = load_cov_sigmas(rascalc_path(tracer, "post"))
        color = TRACER_COLORS.get(tracer.id, "0.3")

        for col, ell in enumerate(ELLS):
            ax = axes[row, col]
            mask_pre = (s_pre >= S_LIM[0]) & (s_pre <= S_LIM[1])
            mask_post = (s_post >= S_LIM[0]) & (s_post <= S_LIM[1])
            if ell in sigma_cov:
                xi_post_interp = np.interp(s_cov, s_post, xi_post[ell])
                band_lo = s_cov**2 * (xi_post_interp - sigma_cov[ell])
                band_hi = s_cov**2 * (xi_post_interp + sigma_cov[ell])
                m_cov = (s_cov >= S_LIM[0]) & (s_cov <= S_LIM[1])
                ax.fill_between(s_cov[m_cov], band_lo[m_cov], band_hi[m_cov],
                                color=color, alpha=0.2, lw=0)
            ax.plot(s_pre[mask_pre], (s_pre**2 * xi_pre[ell])[mask_pre],
                    color=color, lw=1.4, ls=":")
            ax.plot(s_post[mask_post], (s_post**2 * xi_post[ell])[mask_post],
                    color=color, lw=1.4)
            ax.axvline(100.0, color="gray", lw=0.6, ls=":", alpha=0.6)
            ax.axhline(0.0, color="gray", lw=0.6, ls="-", alpha=0.4)
            if row == nrows - 1:
                ax.set_xlabel(r"$s$  [Mpc/$h$]")
            ax.set_ylabel(rf"$s^2\, \xi_{ell}(s)$  [(Mpc/$h$)$^2$]")
            ax.set_title(rf"{tracer.id.upper()}  $\ell = {ell}$")
            ax.set_xlim(*S_LIM)
            ax.grid(True, alpha=0.2)
            if row == 0 and ell == 0:
                ax.legend(handles=[
                    Line2D([0], [0], color="black", lw=1.4, ls=":", label="pre-recon"),
                    Line2D([0], [0], color="black", lw=1.4, label="post-recon (RecSym)"),
                    Patch(facecolor="black", alpha=0.2, label=r"RascalC $\pm1\sigma$"),
                ], fontsize=9, loc="best")

    fig.suptitle(f"DESI DR1 xi(s) multipoles  —  universe: {args.universe}",
                 fontsize=12)
    fig.tight_layout()

    plot_dir = str(cli_out_dir) if cli_out_dir is not None else f"{out_dir}/xi_multipoles_plot"
    os.makedirs(plot_dir, exist_ok=True)
    fig_path = f"{plot_dir}/xi_multipoles_plot.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"Saved {fig_path}")


if __name__ == "__main__":
    main()
