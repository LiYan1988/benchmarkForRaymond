#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=00:100:00
#SBATCH --output=output_test1
#SBATCH --error=error_test1
#SBATCH --partition=serial
#SBATCH --account=maite_group

module load matlab/R2016a
module load gurobi/6.5.1

matlab -nodesktop -r "addpath(genpath('/scratch/ly6j/YALMIP')); parpool(10); test_bm_batch1_5demands; exit" -logfile matlab_output1
