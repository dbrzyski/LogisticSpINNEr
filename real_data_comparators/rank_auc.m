function a = rank_auc(score, y)
%RANK_AUC Area under the ROC curve from ranks (Mann-Whitney), ties counted as 1/2.
y = logical(y(:));
score = score(:);
n1 = sum(y);
n0 = numel(y) - n1;
if n1 == 0 || n0 == 0
    a = NaN;
    return
end
r = tiedrank(score);
a = (sum(r(y)) - n1 * (n1 + 1) / 2) / (n1 * n0);
end
