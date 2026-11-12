function A = generate_relocation_standardized_design(p, n, dataSeed)
%GENERATE_RELOCATION_STANDARDIZED_DESIGN Reproduce the standardized matrix design.

previousState = rng;
restoreState = onCleanup(@() rng(previousState)); %#ok<NASGU>
rng(dataSeed, 'twister');

rawA = randn([p, p, n]);
rawA = rawA + permute(rawA, [2, 1, 3]);
standardizedA = zscore(reshape(rawA, [p^2, n])')';
diagonalIndices = logical(reshape(eye(p), [p^2, 1]));
standardizedA(diagonalIndices, :) = 0;
A = reshape(standardizedA, [p, p, n]);
end
