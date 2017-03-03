#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=00:100:00
#SBATCH --output=output_test2
#SBATCH --error=error_test2
#SBATCH --partition=serial
#SBATCH --account=maite_group

module load matlab/R2015a
module load gurobi/6.5.1

matlab -nodesktop -r "test_bm_batch1_5demands2;" -logfile matlab_output2
