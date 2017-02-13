n = 14;
p = combnk(1:n, 2);
N = 15;
d = randi(size(p, 1), [1, N]);
m = p(d, :);
r = randi([30, 100], [N, 1]);
m(:, end+1) = r;
save('demands15.mat', 'm')