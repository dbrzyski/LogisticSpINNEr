%% Clear the workspace
clear

%% Paths
here = fileparts(mfilename('fullpath'));
SampleSizeRes_Path = fullfile(here, '..', 'f_sample_size', 'results.mat');
SignalStrengthRes_Path = fullfile(here, '..', 'g_signal_strength', 'results_sigStr_sweep_relative_iter50.mat');
outDir = fullfile(here, 'output');
if ~isfolder(outDir), mkdir(outDir); end
SampleSizeSavepath = fullfile(outDir, 'sampleSize_performance.pdf');
SignalStrengthSavepath = fullfile(outDir, 'signalStr_performance.pdf');

%% Signal used in simulation
B1 = ones(5,5); B2 = -ones(5,5); B3 = ones(4,4); sNods = size(B1,1) + size(B2,1) + size(B3,1);
B_true = blkdiag(zeros(5,5), B1, zeros(6,6), B2, zeros(7,7), B3, zeros(60 - sNods - 18, 60- sNods - 18));
B_Fnorm = norm(B_true,'F');

%% Read the data
data1 = load(SampleSizeRes_Path);
data2 = load(SignalStrengthRes_Path);

%% Generate plots
ConsideredMethods = {'Logistic SpINNEr', 'SpINNEr', 'Logistic LASSO', 'Logistic NUCLEAR'};
MethodOrderInPlot = {'SpINNEr', 'Logistic NUCLEAR', 'Logistic LASSO', 'Logistic SpINNEr'};

% find permutation for the results reflecting the methods order
[~,orderPerm] = ismember(MethodOrderInPlot, ConsideredMethods);

% signal strength
%----------------------------------------------------------------------------
iter = 50;
nFixed = 300;
sigStrGrid = [0.1, 0.25, 0.5, 1, 1.5, 2]; 
nMethods = length(ConsideredMethods);

ACC = data2.ACC;
ACC = ACC(:,:,orderPerm); % adjust the order of methods in plots
mean_ACC = squeeze(mean(ACC, 2));
sem_ACC  = squeeze(std(ACC, 0, 2)) ./ sqrt(iter); % Standard Error

figure('Units', 'inches', 'Position', [1, 1, 8.5, 5], 'Color', 'w');
hold on;
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
colors = colors(orderPerm,:);
markers = {'o', 's', 'd', '^'};

for ee = 1:nMethods
    errorbar(sigStrGrid, mean_ACC(:,ee), 2*sem_ACC(:,ee), ...
        'Color', colors(ee,:), 'LineWidth', 1.9, ...
        'Marker', markers{ee}, 'MarkerFaceColor', colors(ee,:), ...
        'MarkerSize', 9, 'CapSize', 10, 'DisplayName', MethodOrderInPlot{ee});
end

% xlabel('Signal Strength (\gamma)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Relative Frobenius Norm Error', 'FontSize', 12, 'FontWeight', 'bold');
% title(sprintf('Effect of Signal Strength (n=%d)', nFixed), 'FontSize', 14);
xlabel('Signal Strength (\gamma)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Frobenius Error', 'FontSize', 12, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 14, 'TickDir', 'out', 'Box', 'on');
ax = gca;
xt = ax.XTick;
ax.FontSize = 17;
ax.XTick = xt;
ax.XTickLabel = string(xt);
ax.LineWidth = 1.7; % thicker tick/box lines
% legend({'Method A', '\vspace{16pt}Method B', 'Method C', 'Method D'}, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 17);
ylim([0.8 1.11])
set(gca, 'TickDir', 'in');
xl = get(gca,'XLabel');
xl.Position(2) = xl.Position(2) + 0.015; 
yl = get(gca,'YLabel');
yl.Position(1) = yl.Position(1) + 0.015; 
ax = gca;                             % your target axes
ax.XRuler.TickLabelGapMultiplier = -0.1;   % < 1 brings X tick labels closer (default = 1)
ax.YRuler.TickLabelGapMultiplier = 0;   % < 1 brings Y tick labels closer

xlabel('Signal Strength (\gamma)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Frobenius Error', 'FontSize', 12, 'FontWeight', 'bold');

%Save as Vector PDF
exportgraphics(gcf, SignalStrengthSavepath, 'ContentType', 'vector');
close all

% sample size
%----------------------------------------------------------------------------
iter = 50;
nGrid = [100, 120, 140, 160, 180, 200, 220, 240];
nMethods = length(ConsideredMethods);

ACC = data1.ACC;
ACC = ACC(:,:,orderPerm); % adjust the order of methods in plots
ACC = ACC./B_Fnorm; % to make the error relative as in other scenario
mean_ACC = squeeze(mean(ACC, 2));
sem_ACC  = squeeze(std(ACC, 0, 2)) ./ sqrt(iter); % Standard Error

figure('Units', 'inches', 'Position', [1, 1, 8.5, 5], 'Color', 'w');
hold on;
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
colors = colors(orderPerm,:);
markers = {'o', 's', 'd', '^'};

for ee = 1:nMethods
    errorbar(nGrid, mean_ACC(:,ee), 2*sem_ACC(:,ee), ...
        'Color', colors(ee,:), 'LineWidth', 1.9, ...
        'Marker', markers{ee}, 'MarkerFaceColor', colors(ee,:), ...
        'MarkerSize', 9, 'CapSize', 10, 'DisplayName', MethodOrderInPlot{ee});
end

xlabel('Sample size (n)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Frobenius Error', 'FontSize', 12, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 14, 'TickDir', 'out', 'Box', 'on');
ax = gca;
xt = ax.XTick;
ax.FontSize = 17;
ax.XTick = xt;
ax.XTickLabel = string(xt);
ax.LineWidth = 1.7; % thicker tick/box lines

legend('Location', 'southwest', 'FontSize', 17);
ylim([0.9 1.01])
set(gca, 'TickDir', 'in');
xl = get(gca,'XLabel');
xl.Position(2) = xl.Position(2) + 0.003; 
yl = get(gca,'YLabel');
yl.Position(1) = yl.Position(1) + 0.015; 
ax = gca;                             % your target axes
ax.XRuler.TickLabelGapMultiplier = -0.1;   % < 1 brings X tick labels closer (default = 1)
ax.YRuler.TickLabelGapMultiplier = 0;   % < 1 brings Y tick labels closer

xlabel('Sample size (n)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Frobenius Error', 'FontSize', 12, 'FontWeight', 'bold');

%Save as Vector PDF
exportgraphics(gcf, SampleSizeSavepath, 'ContentType', 'vector');
close all
