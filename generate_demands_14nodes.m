% generate 100 sets of demands, each set has 100 demands
% the network has 14 nodes
% node index starts from 1

% reset random seed such that we get the same results every time
rng(0)

Nsimu = 100;
demands = cell(Nsimu, 1);
n = 14; % number of nodes
p = combnk(1:n, 2);
N = 100; % number of demands

for i=1:Nsimu
    d = randi(size(p, 1), [N, 1]);
    demands{i} = p(d, :);
    r = randi([30, 100], [N, 1]); % bandwidths are in [30, 100] GHz
    demands{i}(:, end+1) = r;
    fileName = sprintf('demands/demands_14nodes_matlab_%d.csv', i);
    csvwrite(fileName, demands{i});
    
    demandPairMatrix = zeros(size(demands{i}, 1), n);
    for j=1:size(demands{i}, 1)
        demandPairMatrix(j, demands{i}(j, 1)) = 1;
        demandPairMatrix(j, demands{i}(j, 2)) = -1;
    end
    fileName = sprintf('demands/demands_14nodes_matlab_%d.mat', i);
    m = demands{i};
    save(fileName, 'm', 'demandPairMatrix');
end
