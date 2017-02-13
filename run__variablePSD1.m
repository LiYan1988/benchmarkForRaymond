clc;
clear;
close all;
yalmip('clear')

%% Generate traffic matrix
freqMax = 40;

%% Run optimization - 1. preprocessing
Fdistance = 1;

i = 1;
for j = 1 : 1
    load('trafficMultiload_30.mat')
    demandPairs = demandPairs(1:10, :);
    demandPairMatrix = demandPairMatrix(1:10, :);
    result(i) = fcn_ofdm_benchmark_variablePSD_gurobi(freqMax, Fdistance, demandPairs, demandPairMatrix);
    display_status = sprintf('\nTraffic %d in batch %d is finished.\n', j);
    disp(display_status);
    i = i + 1;
    if ~exist('results', 'dir')
        mkdir('results')
    end
    fn = sprintf('results/variablePSD_N30_%d.mat', i);
    save(fn)
end

