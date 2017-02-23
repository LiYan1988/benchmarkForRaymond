% generate 100 sets of demands, each set has 100 demands
% the network has 24 nodes
% node index starts from 1

% reset random seed such that we get the same results every time
rng(0)

Nsimu = 100;
demands = cell(Nsimu, 1);
n = 24; % number of nodes
p = combnk(1:n, 2);
N = 100; % number of demands

for i=1:Nsimu
    d = randi(size(p, 1), [N, 1]);
    demands{i} = p(d, :);
    r = randi([30, 100], [N, 1]); % bandwidths are in [30, 100] GHz
    demands{i}(:, end+1) = r;
    fileName = sprintf('demands/demands_24nodes_matlab_%d.csv', i);
    csvwrite(fileName, demands{i});
end



% demandPairs = m;
% demandPairMatrix = zeros(size(demands, 1), n);
% for i=1:size(demandPairs, 1)
%     demandPairMatrix(i, demandPairs(i, 1)) = 1;
%     demandPairMatrix(i, demandPairs(i, 2)) = -1;
% end