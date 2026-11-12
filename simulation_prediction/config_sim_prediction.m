function cfg = config_sim_prediction(smokeTest)
%CONFIG_SIM_PREDICTION Settings for the simulation study of out-of-sample prediction.
%   cfg = config_sim_prediction()      full experiment
%   cfg = config_sim_prediction(true)  quick pipeline check (2 replicates, 3 x 3 grids)
%
%   Data-generating model and method settings follow Fig. 1(a)-(e)
%   (fig1_simulations/a-e_heatmaps/getEstimates.m): n = 300, gamma = 0.5, p = 60, beta = 1.

if nargin < 1
    smokeTest = false;
end

here = fileparts(mfilename('fullpath'));
cfg.codeRoot    = fileparts(here);                                   % repository root
cfg.softwareDir = fullfile(cfg.codeRoot, 'software_84509f9');        % SpINNEr @ 84509f9
cfg.smokeTest   = smokeTest;

% Data-generating model
cfg.n        = 300;        % training sample size
cfg.nTest    = 5000;       % independent test sample size (same model)
cfg.p        = 60;         % matrix dimension
cfg.gamma    = 0.5;        % signal strength of the blocks in B_true
cfg.betaTrue = 1;          % intercept
cfg.nReps    = 50;         % replicates
cfg.seed     = 20261023;   % replicate r uses seeds seed + r (training) and seed + 100000 + r (test)

B1 = cfg.gamma * ones(5, 5);
B2 = -cfg.gamma * ones(5, 5);
B3 = cfg.gamma * ones(4, 4);
sNods = size(B1, 1) + size(B2, 1) + size(B3, 1);
cfg.Btrue = blkdiag(zeros(5, 5), B1, zeros(6, 6), B2, zeros(7, 7), B3, ...
    zeros(cfg.p - sNods - 18, cfg.p - sNods - 18));

cfg.methodNames = {'LogisticSpINNEr', 'SpINNEr', 'LogisticLASSO', 'LogisticNUCLEAR'};

% Method settings copied from getEstimates.m (Fig. 1)
logistic = struct( ...
    'Family', 'Binomial', ...
    'Method', 'CV', ...
    'UseParallel', true, ...
    'gridLengthN', 10, ...
    'gridLengthL', 10, ...
    'LambdaL', [], ...
    'LambdaN', [], ...
    'maxIters', 100, ...
    'displayStatus', false);

settings = cell(numel(cfg.methodNames), 1);
settings{1} = logistic;                                % Logistic SpINNEr
settings{2} = logistic;  settings{2}.Family  = 'Gaussian';   % SpINNEr (least squares on 0/1)
settings{3} = logistic;  settings{3}.LambdaN = 0;      % Logistic LASSO
settings{4} = logistic;  settings{4}.LambdaL = 0;      % Logistic NUCLEAR

if smokeTest
    cfg.nReps = 2;
    cfg.nTest = 1000;
    for m = 1:numel(settings)
        settings{m}.gridLengthN = 3;
        settings{m}.gridLengthL = 3;
        settings{m}.maxIters    = 20;
        settings{m}.UseParallel = false;
    end
    cfg.resultsDir = fullfile(here, 'results_smoke');
else
    cfg.resultsDir = fullfile(here, 'results');
end
cfg.settings = settings;
end
