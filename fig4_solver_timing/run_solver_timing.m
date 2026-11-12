function run_solver_timing(smokeTest)
%RUN_SOLVER_TIMING Computation time of Logistic SpINNEr with the three IRLS solvers.
%   run_solver_timing()       full experiment
%   run_solver_timing(true)   quick pipeline check
%
%   For each of cfg.nReps data sets per (p, n), one Logistic SpINNEr fit at the stored lambdas
%   is timed with each SP1_SolveMethod (the order rotates between replicates). The time of the
%   setup alone (SVD; K = A*A') is recorded as well. Each (p, n, replicate) result is saved
%   immediately; finished ones are skipped when the function is run again.

if nargin < 1 || isempty(smokeTest)
    smokeTest = false;
end
cfg = config_solver_timing(smokeTest);

addpath(cfg.softwareDir, '-begin');
expected = fullfile(cfg.softwareDir, 'LogisticSpinner.m');
assert(strcmp(which('LogisticSpinner'), expected), ...
    'LogisticSpinner resolves to %s instead of %s.', which('LogisticSpinner'), expected);
assert(all(isfield(ParseArgumentsFit({}), {'SP1_SolveMethod', 'singleThreadCalc', 'eigSymmetric'})), ...
    'This SpINNEr version lacks the SP1_SolveMethod / singleThreadCalc / eigSymmetric options.');
assert(all(cfg.pGrid .* (cfg.pGrid - 1) / 2 > cellfun(@max, cfg.nGrids)), ...
    'The SVD code requires n < p(p-1)/2.');

if ~isfolder(cfg.resultsDir)
    mkdir(cfg.resultsDir);
end
cfg.nThreads = maxNumCompThreads;
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
save(fullfile(cfg.resultsDir, sprintf('config_%s.mat', stamp)), 'cfg');
copyfile(fullfile(cfg.softwareDir, 'ProxyLogistic.m'), ...
    fullfile(cfg.resultsDir, sprintf('ProxyLogistic_used_%s.m', stamp)));   % not yet in a commit

% Warm-up (first calls include JIT compilation)
[A, y, X] = generate_data(cfg, 50, cfg.pGrid(1), 999);
for m = 1:numel(cfg.methods)
    LogisticSpinner(y, X, A, 1, 1, ones(cfg.pGrid(1)) - eye(cfg.pGrid(1)), ...
        ParseArgumentsFit({cfg.fit, 'SP1_SolveMethod', cfg.methods{m}, 'maxIters', 2}));
end

for i = 1:numel(cfg.pGrid)
    p = cfg.pGrid(i);
    W = ones(p) - eye(p);
    for n = cfg.nGrids{i}
        lambdas = load(fullfile(cfg.lambdaDir, sprintf('lambdas_p%02d_n%03d.mat', p, n)), 'lambdas').lambdas;
        for r = 1:cfg.nReps
            outFile = fullfile(cfg.resultsDir, sprintf('time_p%02d_n%03d_rep%02d.mat', p, n, r));
            if isfile(outFile)
                continue
            end
            [A, y, X] = generate_data(cfg, n, p, r);

            % Setup alone, as computed inside LogisticSpinner
            idxs = logical(reshape(triu(ones(p, p), 1), [p^2, 1]));
            Avec = reshape(A, [p^2, n])';
            Avec = 2*Avec(:, idxs);
            nt = maxNumCompThreads(1);
            tic;
            [U, S, V] = svd(Avec, 'econ'); %#ok<ASGLU>
            svdSeconds = toc;
            tic;
            K = Avec*Avec'; %#ok<NASGU>
            kSeconds = toc;
            maxNumCompThreads(nt);

            % Timed fits; order rotates to balance drift
            order = circshift(1:numel(cfg.methods), r - 1);
            fits = struct();
            for m = order
                method = cfg.methods{m};
                params = ParseArgumentsFit({cfg.fit, 'SP1_SolveMethod', method});
                fprintf('[%s] p = %d, n = %d, replicate %d, %s ...\n', ...
                    char(datetime('now')), p, n, r, method);
                tic;
                out = LogisticSpinner(y, X, A, lambdas.N, lambdas.L, W, params);
                elapsed = toc;
                fprintf('   %.2f s, %d ADMM iterations\n', elapsed, out.count);
                fits.(method).seconds   = elapsed;
                fits.(method).admmIters = out.count;
                fits.(method).irlsSteps = sum(out.LogProxyCounts);
                fits.(method).optVal    = out.optVal;
                fits.(method).B         = out.B;
                fits.(method).beta      = out.beta;
            end

            result = struct();
            result.p           = p;
            result.n           = n;
            result.replicate   = r;
            result.lambdaN     = lambdas.N;
            result.lambdaL     = lambdas.L;
            result.svdSeconds  = svdSeconds;
            result.kSeconds    = kSeconds;
            ref = fits.woodbury;
            for m = 1:numel(cfg.methods)
                method = cfg.methods{m};
                result.seconds.(method)   = fits.(method).seconds;
                result.admmIters.(method) = fits.(method).admmIters;
                result.irlsSteps.(method) = fits.(method).irlsSteps;
                result.optVal.(method)    = fits.(method).optVal;
                result.relDiffB.(method)  = norm(fits.(method).B - ref.B, 'fro') / max(norm(ref.B, 'fro'), eps);
                result.diffBeta.(method)  = norm(fits.(method).beta - ref.beta);
            end
            result.nnzB        = nnz(triu(ref.B, 1));
            result.order       = cfg.methods(order);
            result.nThreads    = maxNumCompThreads;
            result.finished    = char(datetime('now'));
            save(outFile, 'result');
            fprintf('   svd %.3f s, woodbury %.3f s, direct %.2f s\n', ...
                result.seconds.svd, result.seconds.woodbury, result.seconds.direct);
        end
    end
end
end
