function summarize_rda_comparators_final(smokeTest)
%SUMMARIZE_RDA_COMPARATORS_FINAL Tables for the edgewise FDR tests and CPM.
%   summarize_rda_comparators_final()       full results in results/
%   summarize_rda_comparators_final(true)   smoke results in results_smoke/
%
%   Writes to the results folder:
%     summary_fdr.csv          edgewise tests and the global-shift diagnostic
%     summary_cpm.csv          CPM by threshold and model
%     summary_prediction.csv   CPM next to the four penalized methods on the same outer folds
%     summary_overlap.csv      Logistic SpINNEr bootstrap edges against both comparators
%     summary_spinner_edges.csv  per-edge detail for those edges
%     summary.txt              everything printed

if nargin < 1
    smokeTest = false;
end
cfg = config_rda_comparators_final(smokeTest);
logFile = fullfile(cfg.resultsDir, 'summary.txt');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
cleanup = onCleanup(@() diary('off'));
fprintf('RDA comparators, final version: summary (%s)\n\n', char(datetime('now')));

%% Edgewise tests ----------------------------------------------------------------------------
F = loadResult(cfg, 'fdr_results.mat');
m = numel(F.p);
y = F.y;
g = F.subjectMeanZ;
[tg, dfg] = edge_tstats(g, y);
pg = 2 * tcdf(-abs(tg), dfg);
fdr = table(m, min(F.p), max(abs(F.t)), nnz(F.p < 0.001), 0.001 * m, ...
    nnz(F.q < cfg.fdrLevels(1)), nnz(F.q < cfg.fdrLevels(2)), min(F.q), ...
    mean(F.t), mean(F.t > 0), nnz(F.t > 2.5), nnz(F.t < -2.5), ...
    mean(g(y == 1)), mean(g(y == 0)), tg, pg, ...
    'VariableNames', {'nEdges', 'MinP', 'MaxAbsT', 'nP_lt_0p001', 'Expected_P_lt_0p001', ...
    'nFDR05', 'nFDR10', 'MinQ', 'MeanT', 'FracTpos', 'nT_gt_2p5', 'nT_lt_m2p5', ...
    'MeanZ_FHpos', 'MeanZ_FHneg', 'tMeanZ', 'pMeanZ'});
fprintf('== Edgewise t-tests (FH+ vs FH-, Fisher z), Benjamini-Hochberg ==\n');
disp(fdr);
writetable(fdr, fullfile(cfg.resultsDir, 'summary_fdr.csv'));

%% CPM ---------------------------------------------------------------------------------------
reps = cell(cfg.nRepeats, 1);
for r = 1:cfg.nRepeats
    reps{r} = loadResult(cfg, sprintf('cpm_rep%02d.mat', r));
end
perm = loadResult(cfg, 'cpm_perm.mat');
cpm = table();
for j = 1:numel(cfg.cpmPThresholds)
    models = reps{1}.byThreshold(j).models;
    for k = 1:numel(models)
        fa1 = reps{1}.byThreshold(j).foldAUC(:, k);
        pooledAll = cellfun(@(R) R.byThreshold(j).pooledAUC(k), reps);
        meanFoldAll = cellfun(@(R) mean(R.byThreshold(j).foldAUC(:, k)), reps);
        nSel = cell2mat(cellfun(@(R) R.byThreshold(j).nSelected, reps, 'UniformOutput', false));
        pPermPooled = NaN;
        pPermMeanFold = NaN;
        if cfg.cpmPThresholds(j) == perm.pThr && strcmp(models{k}, perm.model)
            nb = numel(perm.nullPooledAUC);
            pPermPooled = (1 + sum(perm.nullPooledAUC >= reps{1}.byThreshold(j).pooledAUC(k))) / (1 + nb);
            pPermMeanFold = (1 + sum(perm.nullMeanFoldAUC >= mean(fa1))) / (1 + nb);
        end
        row = table(cfg.cpmPThresholds(j), string(models{k}), ...
            mean(fa1), std(fa1) / sqrt(numel(fa1)), reps{1}.byThreshold(j).pooledAUC(k), ...
            pPermMeanFold, pPermPooled, ...
            mean(meanFoldAll), mean(pooledAll), min(pooledAll), max(pooledAll), ...
            median(nSel(:, 1)), median(nSel(:, 2)), ...
            'VariableNames', {'pThreshold', 'Model', 'Rep1MeanFoldAUC', 'Rep1SEFoldAUC', ...
            'Rep1PooledAUC', 'PermP_MeanFold', 'PermP_Pooled', 'AllRepsMeanFoldAUC', ...
            'AllRepsMeanPooledAUC', 'AllRepsMinPooledAUC', 'AllRepsMaxPooledAUC', ...
            'MedianPosEdges', 'MedianNegEdges'});
        cpm = [cpm; row]; %#ok<AGROW>
    end
end
fprintf('== CPM, nested CV (%d repeats x %d folds; permutation test on repeat 1, %d permutations) ==\n', ...
    cfg.nRepeats, cfg.nOuterFolds, numel(perm.nullPooledAUC));
disp(cpm);
fprintf('Permutation null (p < %g, combined): pooled AUC mean %.3f, 95th percentile %.3f; mean fold AUC 95th percentile %.3f\n\n', ...
    perm.pThr, mean(perm.nullPooledAUC), prctile(perm.nullPooledAUC, 95), prctile(perm.nullMeanFoldAUC, 95));
writetable(cpm, fullfile(cfg.resultsDir, 'summary_cpm.csv'));

if cfg.smokeTest
    return
end

%% Prediction on the shared outer folds (repeat 1), only if ../RDA_prediction is present -----------
if isfile(cfg.predictionSummary)
    prev = readtable(cfg.predictionSummary, 'TextType', 'string');
    main = cpm(cpm.pThreshold == cfg.cpmMain & cpm.Model == "combined", :);
    pred = [table(prev.Method, prev.MeanFoldAUC, prev.SEFoldAUC, prev.PooledAUC, ...
            'VariableNames', {'Method', 'MeanFoldAUC', 'SEFoldAUC', 'PooledAUC'}); ...
            table("CPM", main.Rep1MeanFoldAUC, main.Rep1SEFoldAUC, main.Rep1PooledAUC, ...
            'VariableNames', {'Method', 'MeanFoldAUC', 'SEFoldAUC', 'PooledAUC'})];
    fprintf('== Prediction on the ../RDA_prediction outer folds (repeat 1) ==\n');
    disp(pred);
    writetable(pred, fullfile(cfg.resultsDir, 'summary_prediction.csv'));
else
    fprintf('(%s not found: prediction comparison skipped)\n\n', cfg.predictionSummary);
end

%% Logistic SpINNEr bootstrap edges against the comparators ---------------------------------------
assert(isfile(cfg.bootstrapFile), ['Missing %s: run ../fig5_fig6_real_data/generate_bootstrap_results.m ' ...
    'and summarize_results.m first.'], cfg.bootstrapFile);
boot = load(cfg.bootstrapFile, 'estimAll');
freq = mean(boot.estimAll ~= 0, 3);
freqEdge = freq(F.edgeIdx)';                      % 1 x m, same order as F.t
spin = freqEdge >= cfg.bootstrapThreshold;
nSpin = nnz(spin);

rankAbsT = zeros(1, m);
[~, order] = sort(abs(F.t), 'descend');
rankAbsT(order) = 1:m;
strong = F.p < 0.001;
nStrong = nnz(strong);
nBoth = nnz(spin & strong);
pHyper = hygecdf(nBoth - 1, m, nStrong, nSpin, 'upper');   % P(X >= nBoth) for a random set

% CPM selection frequency of every edge over all training folds of the main setting
cpmCount = zeros(1, m);
nFits = 0;
for r = 1:cfg.nRepeats
    j = find(reps{r}.pThresholds == cfg.cpmMain, 1);
    for kk = 1:size(reps{r}.byThreshold(j).selected, 1)
        idx = [reps{r}.byThreshold(j).selected{kk, 1}(:); reps{r}.byThreshold(j).selected{kk, 2}(:)];
        cpmCount(idx) = cpmCount(idx) + 1;
        nFits = nFits + 1;
    end
end
cpmFreq = cpmCount / nFits;
rho = corr(freqEdge', abs(F.t'), 'Type', 'Spearman');

ov = table(nSpin, nnz(spin & F.rows' <= 200 & F.cols' <= 200), median(rankAbsT(spin)), ...
    min(rankAbsT(spin)), nStrong, nBoth, nSpin * nStrong / m, pHyper, nnz(spin & F.q < 0.05), ...
    nnz(spin & sign(F.t) > 0), median(cpmFreq(spin)), mean(cpmFreq(~spin)), rho, ...
    'VariableNames', {'nSpinnerEdges', 'nCortical', 'MedianRankAbsT', 'BestRankAbsT', ...
    'nEdges_p_lt_0p001', 'nSpinner_p_lt_0p001', 'ExpectedIfRandom', 'HypergeomP', ...
    'nSpinnerFDR05', 'nSpinnerTpos', 'MedianCPMSelFreq', 'MeanCPMSelFreqOthers', 'SpearmanFreqVsAbsT'});
fprintf('== Logistic SpINNEr bootstrap edges (>= %.0f%%) against the comparators ==\n', ...
    100 * cfg.bootstrapThreshold);
disp(ov);
writetable(ov, fullfile(cfg.resultsDir, 'summary_overlap.csv'));

% Per-edge detail with network labels (regions 201-214 are subcortical)
parc = readtable(cfg.parcellationFile, 'TextType', 'string');
label = @(i) networkLabel(i, parc);
r1 = F.rows(spin); c1 = F.cols(spin);
coef = median(boot.estimAll, 3);
edges = table(r1, c1, arrayfun(label, r1), arrayfun(label, c1), freqEdge(spin)', ...
    coef(F.edgeIdx(spin)), F.t(spin)', F.p(spin)', F.q(spin)', rankAbsT(spin)', cpmFreq(spin)', ...
    'VariableNames', {'RegionA', 'RegionB', 'NetworkA', 'NetworkB', 'BootFreq', 'MedianCoef', ...
    't', 'p', 'qBH', 'RankAbsT', 'CPMSelFreq'});
edges = sortrows(edges, 'RankAbsT');
fprintf('Per-edge detail (sorted by |t| rank among %d edges):\n', m);
disp(edges);
writetable(edges, fullfile(cfg.resultsDir, 'summary_spinner_edges.csv'));
end

function s = networkLabel(i, parc)
if i <= height(parc)
    s = parc.descrip1(parc.node == i);
else
    s = "subcortical";
end
end

function R = loadResult(cfg, name)
path = fullfile(cfg.resultsDir, name);
assert(isfile(path), 'Missing result file %s; run run_rda_comparators_final first.', path);
S = load(path, 'result');
R = S.result;
end
