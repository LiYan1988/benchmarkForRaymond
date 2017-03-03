% simulation benchmark for Raymond, 5 demands, 20 simulations
addpath(genpath('/scratch/ly6j/YALMIP'));
addpath(genpath('/share/apps/gurobi/6.5.1/matlab'));

parpool(10);

%% Generate traffic matrix
freqMax = 8000;

%% Run optimization - 1. preprocessing
nDemands = 5;
idx = 1;

demandName = sprintf('../demands/demands_14nodes_matlab_%d.mat', idx);
load(demandName);
demandPair = m(1:nDemands, :);
demandPairMatrix = demandPairMatrix(1:nDemands, :);
clear m;
tic;
result(idx) = fcn_bm(freqMax, demandPair, demandPairMatrix, 15, 1);
runtime = toc;
fprintf('\nTraffic %d is finished using %.2f seconds.\n', idx, runtime);
if ~exist('results', 'dir')
    mkdir('results')
end
resultName = sprintf('results/benchmark_N%d_%d.mat', nDemands, idx);
save(resultName)

exit;

