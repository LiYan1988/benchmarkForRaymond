% simulation benchmark for Raymond, 5 demands, 20 simulations
addpath(genpath('/scratch/ly6j/YALMIP'));
addpath(genpath('/share/apps/gurobi/6.5.1/matlab'));

parpool(10); 

%% Generate traffic matrix
freqMax = 800;

%% Run optimization - 1. preprocessing
nDemands = 5;

for j = 1 : 1
    demandName = sprintf('demands_14nodes_matlab_%d.mat', j);
    load(demandName);
    demandPair = m(1:nDemands, :);
    demandPairMatrix = demandPairMatrix(1:nDemands, :);
    clear m;
    tic;
    result(j) = fcn_bm(freqMax, demandPair, demandPairMatrix, 15, 1);
    runtime = toc;
    fprintf('\nTraffic %d is finished using %.2f seconds.\n', j, runtime);
    if ~exist('results', 'dir')
        mkdir('results')
    end
    resultName = sprintf('results/benchmark_N%d_%d.mat', nDemands, j);
    save(resultName)
end

exit;
