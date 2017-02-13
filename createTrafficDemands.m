function [ trafficDemand, vin ] = createTrafficDemands( nodeNum, ...
        demandTotalPair, trafficDemandIntensity )
    %Creates random traffic demands
    %   Detailed explanation goes here
    
    % randomly choose demandNum pairs from the all-to-all pairs
    sdPair = nchoosek(1:nodeNum, 2);
    % first sample without replacement, then with replacement if
    % demandTotalPair > nchoosek(nodeNum, 2) 
    vidx = randperm(size(sdPair, 1), min(demandTotalPair, size(sdPair, 1))); 
    if demandTotalPair > size(sdPair, 1)
        vidx = [vidx, randi(size(sdPair, 1), 1, demandTotalPair-size(sdPair, 1))];
    end
    % count values of each node pair
    [pairCounts, pairIdx] = hist(vidx, unique(vidx));
    % generate random slot demands
    datarate = unifrnd(trafficDemandIntensity(1), trafficDemandIntensity(2), ...
        demandTotalPair, 1);
    trafficDemand = [sdPair(vidx, :), datarate, pairCounts(vidx)'];
    trafficDemand = sortrows(trafficDemand);
    
    % Create traffic demand matrix v_{i, n}
    rowidx = 1:demandTotalPair;
    rowidx = [rowidx; rowidx];
    rowidx = rowidx(:);
    colidx = trafficDemand(:, 1:2);
    colidx = colidx.';
    colidx = colidx(:);
    vin = full(sparse(rowidx, colidx, 1, demandTotalPair, nodeNum));
end

