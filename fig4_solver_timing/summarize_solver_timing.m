function summary = summarize_solver_timing(smokeTest)
%SUMMARIZE_SOLVER_TIMING Fit times of the three IRLS solvers per (p, n).
%   Median fit times, geometric-mean speed-ups, time per IRLS step, setup times and agreement
%   of the estimates. Writes solver_timing_summary.csv to cfg.resultsDir.

if nargin < 1
    smokeTest = false;
end
cfg = config_solver_timing(smokeTest);

files = dir(fullfile(cfg.resultsDir, 'time_p*_n*_rep*.mat'));
res = arrayfun(@(f) load(fullfile(f.folder, f.name), 'result').result, files);
assert(~isempty(res), 'No results in %s.', cfg.resultsDir);

pAll = [res.p]'; nAll = [res.n]';
sec = [res.seconds]; adm = [res.admmIters]; irl = [res.irlsSteps];
rdB = [res.relDiffB]; dBt = [res.diffBeta];
tS = [sec.svd]'; tW = [sec.woodbury]'; tD = [sec.direct]';
iS = [irl.svd]'; iW = [irl.woodbury]'; iD = [irl.direct]';
gm = @(x) exp(mean(log(x)));

cells = unique([pAll, nAll], 'rows');
rows = cell(size(cells, 1), 18);
for c = 1:size(cells, 1)
    k = pAll == cells(c, 1) & nAll == cells(c, 2);
    rc = res(k);
    rows(c, :) = {cells(c, 1), cells(c, 2), numel(rc), rc(1).lambdaN, rc(1).lambdaL, ...
        median(tS(k)), median(tW(k)), median(tD(k)), ...
        gm(tD(k) ./ tW(k)), gm(tD(k) ./ tS(k)), gm(tS(k) ./ tW(k)), ...
        median(tS(k) ./ iS(k)), median(tW(k) ./ iW(k)), median(tD(k) ./ iD(k)), ...
        median([rc.svdSeconds]), median([rc.kSeconds]), median(iW(k)), median([rc.nnzB])};
end
summary = cell2table(rows, 'VariableNames', {'p', 'n', 'nReps', 'lambdaN', 'lambdaL', ...
    'medSecSVD', 'medSecWoodbury', 'medSecDirect', ...
    'gmDirectOverWoodbury', 'gmDirectOverSVD', 'gmSVDOverWoodbury', ...
    'secPerIrlsSVD', 'secPerIrlsWoodbury', 'secPerIrlsDirect', ...
    'medSvdSeconds', 'medKSeconds', 'medIrlsSteps', 'medNnzB'});
disp(summary);
fprintf('ADMM iterations identical across methods: %d of %d data sets\n', ...
    sum([adm.svd] == [adm.woodbury] & [adm.direct] == [adm.woodbury]), numel(res));
fprintf('IRLS steps identical across methods: %d of %d data sets\n', ...
    sum(iS == iW & iD == iW), numel(res));
fprintf('Largest relative difference in B vs woodbury: svd %.1e, direct %.1e\n', ...
    max([rdB.svd]), max([rdB.direct]));
fprintf('Largest difference in beta vs woodbury: svd %.1e, direct %.1e\n', ...
    max([dBt.svd]), max([dBt.direct]));
fprintf('Overall geometric-mean speed-up: direct/woodbury %.1f, direct/svd %.1f, svd/woodbury %.2f\n', ...
    gm(tD ./ tW), gm(tD ./ tS), gm(tS ./ tW));
writetable(summary, fullfile(cfg.resultsDir, 'solver_timing_summary.csv'));
end
