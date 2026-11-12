function run_rda_comparators_final(analyses, smokeTest)
%RUN_RDA_COMPARATORS_FINAL Edgewise FDR tests and CPM on the AUD connectivity data.
%   run_rda_comparators_final()              both comparators, full experiment
%   run_rda_comparators_final({'fdr'})       edgewise tests only
%   run_rda_comparators_final({'cpm'})       CPM nested CV and its permutation test only
%   run_rda_comparators_final([], true)      quick pipeline check (results_smoke/)
%
%   Data are loaded by load_aud_data exactly as in the bootstrap analysis (non-GSR, low sucrose,
%   subject 119 dropped, n = 160, 214 regions, symmetric with zero diagonals); edge values are
%   Fisher z-transformed. Every result file is saved when complete and skipped on re-run.

if nargin < 2 || isempty(smokeTest)
    smokeTest = false;
end
allAnalyses = {'fdr', 'cpm'};
if nargin < 1 || isempty(analyses)
    analyses = allAnalyses;
end
assert(all(ismember(analyses, allAnalyses)), 'Unknown analysis name.');

cfg = config_rda_comparators_final(smokeTest);
here = fileparts(mfilename('fullpath'));
for f = {'load_aud_data', 'edge_tstats', 'rank_auc', 'cpm_fold'}
    assert(strcmp(which(f{1}), fullfile(here, [f{1} '.m'])), ...
        '%s does not resolve to this folder; run from %s.', f{1}, here);
end
for f = {'GSRandNonGSRarrays.mat', 'gen_hist_alcohol.csv'}
    assert(isfile(fullfile(cfg.dataDir, f{1})), 'Data file %s not found in %s.', f{1}, cfg.dataDir);
end
if ~isfolder(cfg.resultsDir)
    mkdir(cfg.resultsDir);
end

[A, y, info] = load_aud_data(cfg);
p = info.p;
n = info.n;
fprintf('Data: n = %d (FH+ %d, FH- %d), p = %d, excluded subjects: %s\n', ...
    n, info.nFHP, info.nFHN, p, mat2str(info.excludedSubjects'));

edgeIdx = find(triu(true(p), 1));
[rows, cols] = ind2sub([p p], edgeIdx);
E = reshape(A, p^2, n)';
E = atanh(min(max(E(:, edgeIdx), -0.999999), 0.999999));   % n x m Fisher z values

if ismember('fdr', analyses)
    runFDR(cfg, E, y, rows, cols, edgeIdx, info);
end
if ismember('cpm', analyses)
    runCPM(cfg, E, y, info);
end
end

%% Edgewise tests with Benjamini-Hochberg control -----------------------------------------------
function runFDR(cfg, E, y, rows, cols, edgeIdx, info)
outFile = fullfile(cfg.resultsDir, 'fdr_results.mat');
if isfile(outFile)
    fprintf('Skipping %s (already done)\n', outFile);
    return
end
fprintf('[%s] Edgewise t-tests with BH control ...\n', char(datetime('now')));
tic;
[t, df] = edge_tstats(E, y);
pval = 2 * tcdf(-abs(t), df);

result = struct();
result.t        = t;
result.df       = df;
result.p        = pval;
result.q        = bh_adjust(pval);
result.edgeIdx  = edgeIdx;
result.rows     = rows;
result.cols     = cols;
result.subjectMeanZ = mean(E, 2);   % subject-level mean connectivity (global-shift diagnostic)
result.y        = y;
result.info     = info;
result.cfg      = cfg;
result.seconds  = toc;
result.matlab   = version;
result.finished = char(datetime('now'));
save(outFile, 'result');
fprintf('  saved %s\n', outFile);
end

function q = bh_adjust(p)
% Benjamini-Hochberg adjusted p-values (q-values).
[ps, order] = sort(p(:));
m = numel(ps);
qs = ps .* m ./ (1:m)';
qs = min(1, flipud(cummin(flipud(qs))));
q = zeros(size(p));
q(order) = qs;
end

%% CPM ---------------------------------------------------------------------------------------
function runCPM(cfg, E, y, info)
n = numel(y);
% Outer partitions: same construction and seed as ../RDA_prediction, so repeat 1 is its split.
rng(cfg.foldSeed, 'twister');
foldId = zeros(n, cfg.nRepeats);
for r = 1:cfg.nRepeats
    partition = cvpartition(y, 'KFold', cfg.nOuterFolds);
    for k = 1:cfg.nOuterFolds
        foldId(test(partition, k), r) = k;
    end
end
if ~cfg.smokeTest && isfile(cfg.predictionFolds)
    saved = load(cfg.predictionFolds, 'foldId', 'y');
    assert(isequal(saved.y, y) && isequal(saved.foldId(:, 1), foldId(:, 1)), ...
        'Repeat 1 does not reproduce the ../RDA_prediction outer folds.');
end
foldFile = fullfile(cfg.resultsDir, 'cpm_outer_folds.mat');
if ~isfile(foldFile)
    save(foldFile, 'foldId', 'y', 'info');
end

models = {'combined', 'pos', 'neg'};
for r = 1:cfg.nRepeats
    outFile = fullfile(cfg.resultsDir, sprintf('cpm_rep%02d.mat', r));
    if isfile(outFile)
        fprintf('Skipping %s (already done)\n', outFile);
        continue
    end
    tic;
    result = struct('repeat', r, 'pThresholds', cfg.cpmPThresholds, 'foldId', foldId(:, r), 'y', y);
    for j = 1:numel(cfg.cpmPThresholds)
        pThr = cfg.cpmPThresholds(j);
        oof = struct('combined', zeros(n, 1), 'pos', zeros(n, 1), 'neg', zeros(n, 1));
        foldAUC = zeros(cfg.nOuterFolds, numel(models));
        nSel = zeros(cfg.nOuterFolds, 2);
        sel = cell(cfg.nOuterFolds, 2);
        for k = 1:cfg.nOuterFolds
            isTest = foldId(:, r) == k;
            out = cpm_fold(E(~isTest, :), y(~isTest), E(isTest, :), pThr);
            for m = 1:numel(models)
                oof.(models{m})(isTest) = out.score.(models{m});
                foldAUC(k, m) = rank_auc(out.score.(models{m}), y(isTest));
            end
            nSel(k, :) = [out.nPos, out.nNeg];
            sel(k, :) = {out.selPos, out.selNeg};
        end
        res = struct();
        res.pThr = pThr;
        res.models = models;
        res.oofScore = oof;
        res.foldAUC = foldAUC;
        res.pooledAUC = cellfun(@(mm) rank_auc(oof.(mm), y), models);
        res.nSelected = nSel;           % per fold: [positive network, negative network]
        res.selected = sel;             % per fold: edge indices (into the upper-triangle list)
        result.byThreshold(j) = res;
    end
    result.seconds = toc;
    result.matlab = version;
    result.finished = char(datetime('now'));
    save(outFile, 'result');
    fprintf('  saved %s (%.1f s)\n', outFile, result.seconds);
end

% Permutation test of the whole CV pipeline (repeat-1 folds, main threshold, combined model).
outFile = fullfile(cfg.resultsDir, 'cpm_perm.mat');
if isfile(outFile)
    fprintf('Skipping %s (already done)\n', outFile);
    return
end
fprintf('[%s] CPM permutation test, %d permutations ...\n', char(datetime('now')), cfg.nPermCPM);
tic;
rng(cfg.cpmPermSeed, 'twister');
nullPooled = zeros(cfg.nPermCPM, 1);
nullMeanFold = zeros(cfg.nPermCPM, 1);
for b = 1:cfg.nPermCPM
    yb = y(randperm(n));
    oofb = zeros(n, 1);
    fa = zeros(cfg.nOuterFolds, 1);
    for k = 1:cfg.nOuterFolds
        isTest = foldId(:, 1) == k;
        out = cpm_fold(E(~isTest, :), yb(~isTest), E(isTest, :), cfg.cpmMain);
        oofb(isTest) = out.score.combined;
        fa(k) = rank_auc(out.score.combined, yb(isTest));
    end
    nullPooled(b) = rank_auc(oofb, yb);
    nullMeanFold(b) = mean(fa);
end
result = struct('pThr', cfg.cpmMain, 'model', 'combined', ...
    'nullPooledAUC', nullPooled, 'nullMeanFoldAUC', nullMeanFold, ...
    'seconds', toc, 'matlab', version, 'finished', char(datetime('now')));
save(outFile, 'result');
fprintf('  saved %s (%.0f s)\n', outFile, result.seconds);
end
