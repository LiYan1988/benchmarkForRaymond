#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-24:0:0
#SBATCH --output=30dmd_%a.stdout
#SBATCH --error=30dmd_%a.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --array=5-20
#SBATCH --mem=64G
module load matlab/R2015a
module load gurobi/6.5.1

echo $TMPDIR

matlab -nodesktop -r "dmd30_${SLURM_ARRAY_TASK_ID}_bpsk;" -logfile matlab_output_${SLURM_ARRAY_TASK_ID}

