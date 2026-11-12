function summary = summarize_sim_prediction(smokeTest)
%SUMMARIZE_SIM_PREDICTION Test-set prediction measures over replicates for every method.
%   Mean and SE over replicates of each measure from prediction_measures.m, and the
%   paired difference of each method from Logistic SpINNEr (same replicates).
%   Writes results/sim_prediction_summary.csv, results/sim_prediction_paired.csv and
%   results/sim_prediction_boxplots.pdf.

if nargin < 1
    smokeTest = false;
end
cfg = config_sim_prediction(smokeTest);
measures = {'AUC', 'oracleAUC', 'excessKL', 'logLoss', 'bayesLogLoss', 'brier', ...
    'bayesBrier', 'probMSE', 'corrEta', 'calSlope', 'relFrobError', 'seconds'};
nMethods = numel(cfg.methodNames);

% values(r, m, k): replicate r, method m, measure k (NaN where not finished)
values = NaN(cfg.nReps, nMethods, numel(measures));
for m = 1:nMethods
    for r = 1:cfg.nReps
        file = fullfile(cfg.resultsDir, sprintf('fit_%s_rep%03d.mat', cfg.methodNames{m}, r));
        if isfile(file)
            res = load(file, 'result').result;
            for k = 1:numel(measures)
                values(r, m, k) = res.(measures{k});
            end
        end
    end
end
done = all(~isnan(values(:, :, strcmp(measures, 'AUC'))), 2);
fprintf('Replicates finished for all methods: %d of %d\n', sum(done), cfg.nReps);
values = values(done, :, :);
nDone = sum(done);

% Mean (SE) per method
rows = cell(nMethods, 1 + 2 * numel(measures));
names = {'Method'};
for k = 1:numel(measures)
    names = [names, measures(k), {['SE_' measures{k}]}]; %#ok<AGROW>
end
for m = 1:nMethods
    rows{m, 1} = cfg.methodNames{m};
    for k = 1:numel(measures)
        v = values(:, m, k);
        rows{m, 2*k}     = mean(v, 'omitnan');
        rows{m, 2*k + 1} = std(v, 'omitnan') / sqrt(sum(~isnan(v)));
    end
end
summary = cell2table(rows, 'VariableNames', names);

% Paired differences from Logistic SpINNEr
paired = {'AUC', 'excessKL', 'corrEta'};
prow = {};
for m = 2:nMethods
    for k = 1:numel(paired)
        kk = strcmp(measures, paired{k});
        d = values(:, m, kk) - values(:, 1, kk);
        if strcmp(paired{k}, 'excessKL')                   % lower is better for KL
            lsBetter = mean(d > 0);
        else
            lsBetter = mean(d < 0);
        end
        prow(end+1, :) = {cfg.methodNames{m}, paired{k}, mean(d), std(d) / sqrt(nDone), ...
            lsBetter, signrank_p(d)}; %#ok<AGROW>
    end
end
pairedTable = cell2table(prow, 'VariableNames', ...
    {'Method', 'Measure', 'MeanDiffVsLogSpINNEr', 'SE', 'ShareLogSpINNErBetter', 'SignrankP'});

disp(summary(:, {'Method', 'AUC', 'SE_AUC', 'oracleAUC', 'excessKL', 'SE_excessKL', ...
    'corrEta', 'calSlope', 'brier', 'bayesBrier', 'relFrobError'}));
disp(pairedTable);
writetable(summary, fullfile(cfg.resultsDir, 'sim_prediction_summary.csv'));
writetable(pairedTable, fullfile(cfg.resultsDir, 'sim_prediction_paired.csv'));

% Box plots
labels = {'Log. SpINNEr', 'SpINNEr', 'Log. LASSO', 'Log. NUCLEAR'};
panels = {'AUC', 'Test AUC'; 'excessKL', 'Excess KL (log scale)'; ...
          'corrEta', 'corr(\eta-hat, \eta)'; 'calSlope', 'Calibration slope'};
figure('Units', 'inches', 'Position', [1 1 9 7], 'Color', 'w');
for k = 1:size(panels, 1)
    subplot(2, 2, k);
    kk = strcmp(measures, panels{k, 1});
    boxplot(values(:, :, kk), 'Labels', labels);
    hold on;
    if strcmp(panels{k, 1}, 'AUC')
        yline(mean(values(:, 1, strcmp(measures, 'oracleAUC'))), 'k--', 'oracle');
        ylim([min(0.5, min(ylim)), 1]);
    elseif strcmp(panels{k, 1}, 'excessKL')
        set(gca, 'YScale', 'log');
    elseif strcmp(panels{k, 1}, 'calSlope')
        yline(1, 'k--');
    end
    title(panels{k, 2}); grid on;
end
sgtitle(sprintf('n = %d, \\gamma = %g, p = %d, %d replicates, test n = %d', ...
    cfg.n, cfg.gamma, cfg.p, nDone, cfg.nTest), 'FontSize', 11);
exportgraphics(gcf, fullfile(cfg.resultsDir, 'sim_prediction_boxplots.pdf'), 'ContentType', 'vector');
end

function p = signrank_p(d)
d = d(~isnan(d));
if numel(d) < 2 || all(d == 0)
    p = NaN;
else
    p = signrank(d);
end
end
