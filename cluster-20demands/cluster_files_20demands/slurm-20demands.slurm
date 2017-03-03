#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-5:0:0
#SBATCH --output=20demands_%a.stdout
#SBATCH --error=20demands_%a.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --array=1-2

module load matlab/R2015a
module load gurobi/6.5.1

echo $TMPDIR

matlab -nodesktop -r "matlab_20demands_${SLURM_ARRAY_TASK_ID};" -logfile matlab_output_${SLURM_ARRAY_TASK_ID}

