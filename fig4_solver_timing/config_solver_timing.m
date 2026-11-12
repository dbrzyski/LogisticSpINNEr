function cfg = config_solver_timing(smokeTest)
%CONFIG_SOLVER_TIMING Settings for the timing study of the IRLS solvers in Logistic SpINNEr (Fig. 4).
%   cfg = config_solver_timing()      full experiment
%   cfg = config_solver_timing(true)  quick pipeline check (2 sample sizes, p = 30, 1 replicate)
%
%   Each fit is a single Logistic SpINNEr fit (both penalties) at fixed lambdas, run with
%   SP1_SolveMethod 'svd', 'woodbury' and 'direct' on the same data, on one computational
%   thread (singleThreadCalc = 2) with eigSymmetric = true. Scenarios (grid, data seeds,
%   lambdas, fit settings) are those of an earlier timing study (SVD against direct).

if nargin < 1
    smokeTest = false;
end

here = fileparts(mfilename('fullpath'));
cfg.softwareDir     = fullfile(fileparts(here), 'software');   % needs SP1_SolveMethod, singleThreadCalc, eigSymmetric
cfg.softwareVersion = 'b4fc36b';
cfg.lambdaDir       = fullfile(here, 'lambdas');   % lambdas_p*_n*.mat selected in the earlier study
cfg.smokeTest       = smokeTest;

% Design
cfg.methods  = {'svd', 'woodbury', 'direct'};
cfg.pGrid    = [30, 45, 60, 75, 90];   % p(p-1)/2 = 435, 990, 1770, 2775, 4005 free entries
cfg.nGrids   = repmat({[50, 100, 150]}, 1, numel(cfg.pGrid));   % sample sizes per p; each n < p(p-1)/2
cfg.nReps    = 5;               % timed data sets per (n, p)
cfg.gamma    = 0.5;             % block signal strength, as in Fig. 1(a)-(e)
cfg.betaTrue = 1;
cfg.seed     = 20261024;        % data set (n, p, r) uses seed + 1e5*p + 100*n + r; r = 0 was the pilot

% Fit settings (Fig. 1 values; maxIters = 100 ADMM iterations), one thread per fit
cfg.fit = struct( ...
    'Family', 'Binomial', ...
    'maxIters', 100, ...
    'singleThreadCalc', 2, ...
    'eigSymmetric', true, ...
    'displayStatus', false);

if smokeTest
    cfg.pGrid  = 30;
    cfg.nGrids = {[50, 100]};
    cfg.nReps = 1;
    cfg.fit.maxIters = 10;
    cfg.resultsDir = fullfile(here, 'results_smoke');
else
    cfg.resultsDir = fullfile(here, 'results');
end
end
