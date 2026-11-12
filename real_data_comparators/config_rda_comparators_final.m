function cfg = config_rda_comparators_final(smokeTest)
%CONFIG_RDA_COMPARATORS_FINAL Settings for the comparators reported with Logistic SpINNEr.
%   cfg = config_rda_comparators_final()      full experiment
%   cfg = config_rda_comparators_final(true)  quick pipeline check (30 ROIs, few permutations)
%
%   Two comparators on the AUD connectivity data (FH+ vs FH-):
%     FDR  edgewise two-sample t-tests with Benjamini-Hochberg control;
%     CPM  Connectome-based Predictive Modeling (Shen et al. 2017), adapted to a binary outcome,
%          evaluated by nested CV whose first repeat is the outer split of ../RDA_prediction.
%   Edge values are Fisher z-transformed correlations.
%   Final version of the exploratory ../RDA_comparators (which also has NBS and raw correlations).

if nargin < 1
    smokeTest = false;
end

% Paths are relative to this folder. The data are not distributed: they are read from
% '../../../AUD data' (the project folder) if it exists, otherwise from '../fig5_fig6_real_data/data'
% (the layout of the public code repository).
here = fileparts(mfilename('fullpath'));
cfg.codeRoot = fileparts(here);
candidates = {fullfile(fileparts(fileparts(cfg.codeRoot)), 'AUD data'), ...
              fullfile(cfg.codeRoot, 'fig5_fig6_real_data', 'data')};
found = find(cellfun(@isfolder, candidates), 1);
if isempty(found)
    found = numel(candidates);   % run_rda_comparators_final reports the missing files
end
cfg.dataDir = candidates{found};
cfg.parcellationFile = fullfile(cfg.dataDir, 'Schaefer_parcellation.csv');   % Yeo labels of regions 1-200

% Logistic SpINNEr bootstrap estimates (214 x 214 x 300), written by ../fig5_fig6_real_data
cfg.bootstrapFile = fullfile(cfg.codeRoot, 'fig5_fig6_real_data', 'results', 'EstimatesBootstrap.mat');

% Optional: the nested-CV prediction experiment of the four penalized methods. If present, repeat 1 of
% the CPM folds is checked against its outer folds and its AUCs are tabulated next to CPM.
cfg.predictionFolds   = fullfile(cfg.codeRoot, 'RDA_prediction', 'results', 'outer_folds.mat');
cfg.predictionSummary = fullfile(cfg.codeRoot, 'RDA_prediction', 'results', 'rda_prediction_summary.csv');
cfg.smokeTest = smokeTest;
cfg.roiSubset = [];                    % [] = all 214 regions, as in the main analysis

% Edgewise tests
cfg.fdrLevels = [0.05 0.10];

% CPM
cfg.cpmPThresholds = [0.01 0.05 0.001];   % first entry is the main setting (Shen et al. 2017)
cfg.cpmMain        = 0.01;
cfg.nOuterFolds    = 5;
cfg.nRepeats       = 10;               % repeat 1 reproduces the ../RDA_prediction partition
cfg.foldSeed       = 20261022;         % same seed as ../RDA_prediction
cfg.nPermCPM       = 5000;             % label permutations of the whole CV pipeline (repeat 1)
cfg.cpmPermSeed    = 20261113;

% Logistic SpINNEr bootstrap selection used for the comparison
cfg.bootstrapThreshold = 0.90;

if smokeTest
    cfg.roiSubset  = 1:30;
    cfg.nRepeats   = 2;
    cfg.nPermCPM   = 20;
    cfg.resultsDir = fullfile(here, 'results_smoke');
else
    cfg.resultsDir = fullfile(here, 'results');
end
end
