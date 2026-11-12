function [t, df] = edge_tstats(E, y)
%EDGE_TSTATS Pooled-variance two-sample t statistic for every column of E.
%   E is n x m (subjects x edges), y is n x 1 with values 0/1. Positive t means a larger
%   mean in group y = 1 (FH+). Vectorized so that it can be called inside permutation loops.

y = logical(y(:));
n1 = sum(y);
n0 = numel(y) - n1;
E1 = E(y, :);
E0 = E(~y, :);
m1 = sum(E1, 1) / n1;
m0 = sum(E0, 1) / n0;
ss = sum((E1 - m1).^2, 1) + sum((E0 - m0).^2, 1);
df = n1 + n0 - 2;
s2 = ss / df;
t = (m1 - m0) ./ sqrt(s2 * (1/n1 + 1/n0));
t(~isfinite(t)) = 0;   % constant edges (zero variance) carry no evidence
end
