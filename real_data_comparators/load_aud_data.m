function [A, y, info] = load_aud_data(cfg)
%LOAD_AUD_DATA Load the AUD connectivity data exactly as in the bootstrap analysis (Figs. 5-6).
%   [A, y, info] = load_aud_data(cfg) reads, from cfg.dataDir,
%     GSRandNonGSRarrays.mat  array AA_NOGSR_low: 214 x 214 x 161 correlation matrices
%                             (no global signal regression, low-sucrose condition);
%     gen_hist_alcohol.csv    column x: 'FHP' (family history positive) or 'FHN' (negative).
%   Subjects whose matrix contains missing values are dropped (subject 119), giving n = 160.
%   The matrices are symmetrized and their diagonals set to zero. y is 1 for FH+ and 0 for FH-.
%   cfg.roiSubset, if not empty, restricts the matrices to those regions (smoke test).

data = load(fullfile(cfg.dataDir, 'GSRandNonGSRarrays.mat'), 'AA_NOGSR_low');
genHist = readtable(fullfile(cfg.dataDir, 'gen_hist_alcohol.csv'));

A0 = data.AA_NOGSR_low;
y0 = double(ismember(genHist.x, 'FHP'));

keep = squeeze(~any(any(isnan(A0), 1), 2));
A = A0(:, :, keep);
y = y0(keep);

if ~isempty(cfg.roiSubset)
    A = A(cfg.roiSubset, cfg.roiSubset, :);
end

[p, ~, n] = size(A);
A = (A + permute(A, [2 1 3])) / 2;
A(repmat(logical(eye(p)), [1 1 n])) = 0;

info.n = n;
info.p = p;
info.nFHP = sum(y == 1);
info.nFHN = sum(y == 0);
info.excludedSubjects = find(~keep);
end
