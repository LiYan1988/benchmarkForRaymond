clc;
clear;
close all;

nSimulations = 2;
nDemands = 5;
modulation = 'bpsk';
folderName = sprintf('cluster_files_%ddemands3_%s', nDemands, modulation);
nNodes = 1;
cpuPerTask = 4;
ntasks_per_node = 1;
nDays = 0;
nHours = 10;
nMinutes = 0;
nSeconds = 0;
mem = 64;
partition = 'economy';
group = 'maite_group';

if ~exist(folderName, 'dir')
    mkdir(folderName)
end
listing = dir('.');

% copy all files into cluster_files
for i=1:length(listing)
    if ~(strcmp(listing(i).name, '.') || strcmp(listing(i).name, '..')...
            || strcmp(listing(i).name, 'template_slurm.slurm')...
            || strcmp(listing(i).name, 'template_matlab.m')...
            || listing(i).isdir...
            || strcmp(listing(i).name, 'template_slurm_array.slurm'))
        copyfile(listing(i).name, folderName)
    end
end

template_matlab = regexp( fileread('template_matlab.m'), '\n', 'split');
for i=1:nSimulations
    template_matlab{11} = sprintf('nDemands = %d;', nDemands);
    template_matlab{12} = sprintf('idx = %d;', i);
    if strcmp(modulation, 'bpsk')
        template_matlab{20} = sprintf('result(idx) = fcn_bm(freqMax, demandPair, demandPairMatrix, 15, %d);', 1);
    elseif strcmp(modulation, 'qpsk')
        template_matlab{20} = sprintf('result(idx) = fcn_bm(freqMax, demandPair, demandPairMatrix, 15, %d);', 2);
    end
    fid = fopen(sprintf('%s/matlab_%ddemands_%d_%s.m', folderName, nDemands, i, modulation), 'w');
    fprintf(fid, '%s\n', template_matlab{:});
    fclose(fid);
end

template_slurm = regexp( fileread('template_slurm_array.slurm'), '\n', 'split');
template_slurm{1} = sprintf('#!/bin/bash');
template_slurm{2} = sprintf('#SBATCH --nodes=%d', nNodes);
% template_slurm{3} = sprintf('#SBATCH --ntasks-per-node=%d', ntasks_per_node);
template_slurm{3} = sprintf('#SBATCH --cpus-per-task=%d', cpuPerTask);
template_slurm{4} = sprintf('#SBATCH --time=%d-%d:%d:%d', nDays, nHours, nMinutes, nSeconds);
template_slurm{5} = sprintf('#SBATCH --output=%ddemands_%%a.stdout', nDemands);
template_slurm{6} = sprintf('#SBATCH --error=%ddemands_%%a.stderr', nDemands);
template_slurm{7} = sprintf('#SBATCH --partition=%s', partition);
template_slurm{8} = sprintf('#SBATCH --account=%s', group);
template_slurm{9} = sprintf('#SBATCH --array=%d-%d', 1, nSimulations);
template_slurm{10} = sprintf('#SBATCH --mem=%dG', mem);
template_slurm{11} = sprintf('module load matlab/R2015a');
template_slurm{12} = sprintf('module load gurobi/6.5.1');
template_slurm{13} = sprintf('');
template_slurm{14} = sprintf('echo $TMPDIR');
template_slurm{15} = sprintf('');
template_slurm{16} = sprintf('matlab -nodesktop -r "matlab_%ddemands_${SLURM_ARRAY_TASK_ID};" -logfile matlab_output_${SLURM_ARRAY_TASK_ID}',nDemands);

fid = fopen(sprintf('%s/slurm-%ddemands.slurm', folderName, nDemands), 'w');
fprintf(fid, '%s\n', template_slurm{:});
fclose(fid);