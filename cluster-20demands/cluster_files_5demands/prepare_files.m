clc;
clear;
close all;

nSimulations = 3;
nDemands = 5;
folderName = sprintf('cluster_files_%ddemands', nDemands);
nNodes = 1;
nTasks = 10;
nDays = 0;
nHours = 10;
nMinutes = 0;
nSeconds = 0;
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
    fid = fopen(sprintf('%s/matlab-%ddemands-%d.m', folderName, nDemands, i), 'w');
    fprintf(fid, '%s\n', template_matlab{:});
    fclose(fid);
end

template_slurm = regexp( fileread('template_slurm_array.slurm'), '\n', 'split');
template_slurm{2} = sprintf('#SBATCH --nodes=%d', nNodes);
template_slurm{3} = sprintf('#SBATCH --ntasks=%d', nTasks);
template_slurm{4} = sprintf('#SBATCH --time=%d-%d:%d:%d', nDays, nHours, nMinutes, nSeconds);
template_slurm{5} = sprintf('#SBATCH --output=%ddemands_%%a.stdout', nDemands);
template_slurm{6} = sprintf('#SBATCH --error=%ddemands_%%a.stderr', nDemands);
template_slurm{7} = sprintf('#SBATCH --partition=%s', partition);

fid = fopen(sprintf('%s/slurm-%ddemands.slurm', folderName, nDemands), 'w');
fprintf(fid, '%s\n', template_slurm{:});
fclose(fid);