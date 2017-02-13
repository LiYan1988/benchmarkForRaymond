clear;
p_Lambda = struct();

rng('default')
maxTraffic = 180;
p_Lambda = randi([18, maxTraffic], 15, 300);

fn = sprintf('p_Lambda.mat');
save(fn, 'p_Lambda');