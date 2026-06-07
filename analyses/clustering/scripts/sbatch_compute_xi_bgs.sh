#!/bin/bash
# xi for BGS (uses existing bgs_full recon NPZ).
#SBATCH --job-name=xi_bgs
#SBATCH --account=m4031
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=logs/xi_bgs_%j.out
#SBATCH --error=logs/xi_bgs_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY="$SLURM_SUBMIT_DIR/.venv/bin/python"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cd analyses/clustering
srun -n "$SLURM_NTASKS" -c "$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    "$PY" -u scripts/compute_xi.py --universe baseline --tracer bgs
