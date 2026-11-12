function run_sim_prediction(methodsToRun, smokeTest)
%RUN_SIM_PREDICTION Simulation study of out-of-sample prediction for all four methods.
%   run_sim_prediction()                         all methods, full experiment
%   run_sim_prediction({'LogisticSpINNEr'})      selected methods only
%   run_sim_prediction([], true)                 quick pipeline check
%
%   For each replicate, training data (n = 300) and an independent test set (nTest = 5000)
%   are drawn from the same model (generate_replicate.m). Each method selects lambda by
%   spinner's own 5-fold CV on the training data; the fit is scored on the test set
%   (prediction_measures.m). Each (replicate, method) result is saved immediately;
%   finished fits are skipped when the function is run again.

if nargin < 2 || isempty(smokeTest)
    smokeTest = false;
end
cfg = config_sim_prediction(smokeTest);
if nargin < 1 || isempty(methodsToRun)
    methodsToRun = cfg.methodNames;
end
assert(all(ismember(methodsToRun, cfg.methodNames)), 'Unknown method name.');

addpath(cfg.softwareDir, '-begin');
expected = fullfile(cfg.softwareDir, 'spinner.m');
assert(strcmp(which('spinner'), expected), ...
    'spinner resolves to %s instead of the pinned copy %s.', which('spinner'), expected);

if ~isfolder(cfg.resultsDir)
    mkdir(cfg.resultsDir);
end
cfgFile = fullfile(cfg.resultsDir, 'config.mat');
if isfile(cfgFile)
    saved = load(cfgFile, 'cfg').cfg;
    assert(isequal(rmfield(saved, 'settings'), rmfield(cfg, 'settings')) && ...
        isequal(saved.settings, cfg.settings), ...
        'Stored configuration differs from the current one; use a new results folder.');
else
    save(cfgFile, 'cfg');
end

for r = 1:cfg.nReps
    methods = find(ismember(cfg.methodNames, methodsToRun));
    outFiles = arrayfun(@(m) fullfile(cfg.resultsDir, ...
        sprintf('fit_%s_rep%03d.mat', cfg.methodNames{m}, r)), methods, 'UniformOutput', false);
    if all(isfile(outFiles))
        continue
    end
    data = generate_replicate(cfg, r);

    for i = 1:numel(methods)
        m = methods(i);
        method = cfg.methodNames{m};
        if isfile(outFiles{i})
            continue
        end
        fprintf('[%s] replicate %d/%d, %s ...\n', char(datetime('now')), r, cfg.nReps, method);

        tic;
        est = fit_method(data.y, data.X, data.A, cfg.settings{m});
        elapsed = toc;

        result = prediction_measures(est, data, method, cfg);
        result.method      = method;
        result.replicate   = r;
        result.B           = sparse(est.B);
        result.beta        = est.beta;
        result.bestLambdaN = est.bestLambdaN;
        result.bestLambdaL = est.bestLambdaL;
        result.trainFH     = mean(data.y);
        result.seconds     = elapsed;
        result.settings    = cfg.settings{m};
        result.matlab      = version;
        result.finished    = char(datetime('now'));
        save(outFiles{i}, 'result');
        fprintf('   done in %.1f min: AUC %.3f (oracle %.3f), excess KL %.4f\n', ...
            elapsed / 60, result.AUC, result.oracleAUC, result.excessKL);
    end
end
end
