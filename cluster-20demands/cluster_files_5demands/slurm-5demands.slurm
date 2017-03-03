#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=0-10:0:0
#SBATCH --output=5demands_%a.stdout
#SBATCH --error=5demands_%a.stderr
#SBATCH --partition=economy
#SBATCH --account=maite_group
#SBATCH --array=1-30

module load matlab/R2015a
module load gurobi/6.5.1

echo $TMPDIR

matlab -nodesktop -r "test_bm_batch1_5demands_${SLURM_ARRAY_TASK_ID};" -logfile matlab_output_${SLURM_ARRAY_TASK_ID}

