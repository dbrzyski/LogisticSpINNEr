function out = cpm_fold(Etrain, ytrain, Etest, pThr)
%CPM_FOLD Connectome-based Predictive Modeling for a binary outcome on one train/test split.
%   Edge selection on the training subjects only: two-sample t-test of each edge between the
%   classes (equivalent to testing the point-biserial correlation), keeping edges with two-sided
%   p < pThr. Positive-network strength = sum of the selected edges with t > 0, negative-network
%   strength = sum of those with t < 0 (Shen et al. 2017). A logistic regression of the class on
%   the strengths is fitted to the training subjects and gives the test scores (linear predictor):
%     score.combined  both strengths;  score.pos  positive network only;  score.neg  negative only.

[t, df] = edge_tstats(Etrain, ytrain);
p = 2 * tcdf(-abs(t), df);
selPos = p < pThr & t > 0;
selNeg = p < pThr & t < 0;

strength = @(E) [sum(E(:, selPos), 2), sum(E(:, selNeg), 2)];
Strain = strength(Etrain);
Stest  = strength(Etest);

out.nPos = nnz(selPos);
out.nNeg = nnz(selNeg);
out.selPos = find(selPos);
out.selNeg = find(selNeg);
out.score.combined = logisticScore(Strain, ytrain, Stest, [1 2]);
out.score.pos      = logisticScore(Strain, ytrain, Stest, 1);
out.score.neg      = logisticScore(Strain, ytrain, Stest, 2);
end

function s = logisticScore(Strain, ytrain, Stest, cols)
% Columns without selected edges are constant (all zero) and are dropped. With no usable
% column the score is constant, which gives AUC = 0.5.
keep = cols(std(Strain(:, cols), 0, 1) > 0);
if isempty(keep)
    s = zeros(size(Stest, 1), 1);
    return
end
mu = mean(Strain(:, keep), 1);
sd = std(Strain(:, keep), 0, 1);
Z  = (Strain(:, keep) - mu) ./ sd;
Zt = (Stest(:, keep) - mu) ./ sd;
warnState = warning('off', 'all');
restore = onCleanup(@() warning(warnState));
b = glmfit(Z, ytrain, 'binomial', 'link', 'logit');
s = [ones(size(Zt, 1), 1), Zt] * b;
end
