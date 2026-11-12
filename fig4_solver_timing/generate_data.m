function [A, y, X] = generate_data(cfg, n, p, r)
%GENERATE_DATA Data set r for sample size n and dimension p (r = 0: lambda-selection pilot).
%   Same model as Fig. 1: symmetric Gaussian matrices, entries z-scored across subjects,
%   zero diagonals, three blocks of strength gamma, y ~ Bernoulli(logistic(<A, B> + beta)).
%   The blocks span 32 nodes; for p < 32, B is cropped to its first p rows and columns.

stream = RandStream('twister', 'Seed', cfg.seed + 1e5*p + 100*n + r);

B1 = cfg.gamma * ones(5, 5);
B2 = -cfg.gamma * ones(5, 5);
B3 = cfg.gamma * ones(4, 4);
sNods = size(B1, 1) + size(B2, 1) + size(B3, 1);
Btrue = blkdiag(zeros(5, 5), B1, zeros(6, 6), B2, zeros(7, 7), B3, ...
    zeros(max(p - sNods - 18, 0)));
Btrue = Btrue(1:p, 1:p);

AA = randn(stream, [p, p, n]);
AA = AA + permute(AA, [2 1 3]);
AA2 = zscore(reshape(AA, [p^2, n])')';
AA2(logical(reshape(eye(p), [p^2, 1])), :) = 0;
A = reshape(AA2, [p, p, n]);

X = ones(n, 1);
eta = AA2' * Btrue(:) + cfg.betaTrue;
y = double(rand(stream, n, 1) < 1 ./ (1 + exp(-eta)));
end
