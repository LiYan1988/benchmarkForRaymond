clc;
clear;
close all;

% create 15 traffic demands
load 15demands.mat
demandPairs = demands;
demandPairMatrix = zeros(size(demands, 1), 14);
for i=1:size(demandPairs, 1)
    demandPairMatrix(i, demandPairs(i, 1)) = 1;
    demandPairMatrix(i, demandPairs(i, 2)) = -1;
end

save('15demandsInput.mat', 'demandPairs', 'demandPairMatrix')