#!/bin/bash
# xi for combined tracer lrg3_elg1 (paper §7.5). Uses the joint-FKP DIY
# in compute_xi.py + combined_tracer.py — no VAC clone.
#SBATCH --job-name=xi_lrg3_elg1
#SBATCH --account=m4031
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=logs/xi_lrg3_elg1_%j.out
#SBATCH --error=logs/xi_lrg3_elg1_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY="$SLURM_SUBMIT_DIR/.venv/bin/python"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cd analyses/clustering
srun -n "$SLURM_NTASKS" -c "$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    "$PY" -u scripts/compute_xi.py --universe baseline --tracer lrg3_elg1
