% baseline 25 demands in DT network
clc;
clear;
close all;

yalmip('clear')

%% Generate traffic matrix
freqMax = 800;

%% Run optimization - 1. preprocessing
Fdistance = 1;

i = 1;
for j = 1 : 1
    load('30demandsInput.mat');
    demandPairs = demandPairs(1:25,:);
    demandPairMatrix = demandPairMatrix(1:25, :);
    tic;
    result(i) = fcn_bm(freqMax, demandPairs, demandPairMatrix, 15);
    runtime = toc;
    display_status = sprintf('\nTraffic %d in batch %d is finished.\n', j);
    disp(display_status);
    i = i + 1;
    if ~exist('results', 'dir')
        mkdir('results')
    end
    fn = sprintf('results/variablePSD_N25_%d.mat', i);
    save(fn)
end
