function [ trafficVector ] = CreateTrafficVector( networkCostMatrix, ...
        cDistance, totalTraffic, subcarrierBandwidth )
    %Create a traffic matrix according to David's JLT paper
    
    %   Input:
    %   networkCostMatrix: distance between every node pair
    %   cDistance: characteristic distance, e.g., inf, 1500km, 3000km
    
    %   Output:
    %   trafficMatrix: traffic matrix   
    
    nodeNum = size(networkCostMatrix, 1);
    distanceShortest = zeros(nodeNum);
    
    for i = 1 : nodeNum
        for j = i+1 : nodeNum
            [~, distanceShortest(i, j)] = ...
                kShortestPath(networkCostMatrix, i, j, 1);
        end
    end
    
    trafficMatrix = (distanceShortest+distanceShortest');
    trafficMatrix = exp(-trafficMatrix/cDistance);
    trafficMatrix = trafficMatrix-diag(diag(trafficMatrix));
    trafficMatrix = totalTraffic*trafficMatrix/sum(trafficMatrix(:));
    trafficMatrix = ceil(trafficMatrix/subcarrierBandwidth/2)*2;
    
    trafficVector = zeros(nchoosek(nodeNum, 2), 1);
    cnt = 1;
    for i = 1 : nodeNum
        for j = i+1 : nodeNum
            trafficVector(cnt) = trafficMatrix(i, j);
            cnt = cnt+1;
        end
    end
    
end

