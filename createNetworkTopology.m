clc;
clear;
close all;

% create network adjacent matrix and network cost matrix
load networkTopology.mat
load benchmarkNetworkSimple.mat
clearvars networkCoordinates

networkCostMatrix = sparse(networkTopology(:, 1), networkTopology(:, 2), networkTopology(:, 3), 14, 14);
networkCostMatrix = full(networkCostMatrix);
networkCostMatrix(networkCostMatrix==0) = Inf;

networkAdjacentMatrix = full(sparse(networkTopology(:, 1), networkTopology(:, 2), 1, 14, 14));

save('GermenNetworkTopology.mat', 'networkCostMatrix', 'networkAdjacentMatrix', 'networkTopology')