clc;
clear;
close all;

% create 15 traffic demands
load 15demands.mat

n = 14;
p = combnk(1:n, 2);
N = 15;

newDemands = demands;
while size(newDemands, 1)<30
    d = randi(size(p, 1));
    m = p(d, :);
    flag = 1;
    for i=1:size(newDemands, 1)
        if m(1)==newDemands(i, 1) && m(2)==newDemands(i, 2)
            flag = 0;
            break
        end
    end
    if flag==1
        newDemands(end+1, :) = zeros(1, 3);
        newDemands(end, 1:2) = m;
        newDemands(end, 3) = randi([30, 100]);
    end
end

%%
demandPairs = newDemands;
demandPairMatrix = zeros(size(demands, 1), 14);
for i=1:size(demandPairs, 1)
    demandPairMatrix(i, demandPairs(i, 1)) = 1;
    demandPairMatrix(i, demandPairs(i, 2)) = -1;
end

save('30demandsInput.mat', 'demandPairs', 'demandPairMatrix')