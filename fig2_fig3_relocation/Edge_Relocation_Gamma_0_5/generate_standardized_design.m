function A = generate_standardized_design(p, n, dataSeed)
%GENERATE_STANDARDIZED_DESIGN Reproduce the matrix-covariate generation code.

previousState = rng;
restoreState = onCleanup(@() rng(previousState)); %#ok<NASGU>
rng(dataSeed, 'twister');

AA = randn([p, p, n]);
AA = AA + permute(AA, [2, 1, 3]);
AA2 = zscore(reshape(AA, [p^2, n])')';
diagonalIndices = logical(reshape(eye(p), [p^2, 1]));
AA2(diagonalIndices, :) = 0;
A = reshape(AA2, [p, p, n]);
end
