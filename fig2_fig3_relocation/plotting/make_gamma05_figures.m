%% Figs. 2 and 3 of the manuscript at gamma = 0.5.
% Run this entire script. Requires MATLAB R2020a or later; no toolboxes.
% Data:
%   Edge_Relocation_Gamma_0_5   blocks (5,5,4): Fig. 2 first point and Fig. 3(a)
%   Block_Relocation_Gamma_0_5  blocks (5,5,5)-(16,16,16): rest of Fig. 2, Fig. 3(b)-(f)
% Both runs use n = 300, p = 60, base seed 20250705 and the same node layout.
% MeanACC is already relative Frobenius error and TwoSEM is already twice the
% standard error: neither is rescaled here.
% Output: vector PDFs in the figures folder beside this script. The manuscript
% uses fig2-intact-block-gamma05.pdf and the one-row fig3-row-gamma05_*.pdf files.

scriptDir = fileparts(mfilename('fullpath'));
dataRoot = fileparts(scriptDir);
figureDir = fullfile(scriptDir, 'figures');
if ~isfolder(figureDir)
    mkdir(figureDir);
end

edge = readSummary(fullfile(dataRoot, 'Edge_Relocation_Gamma_0_5', ...
    'results', 'swap_pilot_ACC_summary.csv'), false);
blocks = readSummary(fullfile(dataRoot, 'Block_Relocation_Gamma_0_5', ...
    'results', 'block_relocation_gamma0p5_summary_rep50_methods.csv'), true);
assert(all(blocks.CoefficientMagnitude == 0.5), 'Expected gamma = 0.5 block results.');

% (5,5,4) intact comes from the edge relocation run (k = 0)
edgeIntact = edge(edge.SwapProportion == 0, {'Method', 'MeanACC', 'TwoSEM'});
edgeIntact.BlockLabel = repmat("5-5-4", height(edgeIntact), 1);
edgeIntact.ActualSwapProportion = zeros(height(edgeIntact), 1);
blockIntact = blocks(blocks.ActualSwapProportion == 0, ...
    {'Method', 'MeanACC', 'TwoSEM', 'BlockLabel', 'ActualSwapProportion'});
blockOrder = ["5-5-4", "5-5-5", "6-6-6", "7-7-7", "8-8-8", "12-12-12", "16-16-16"];

% Exact order, RGB values, and markers from generatePlots.m after orderPerm.
style.Keys = ["SpINNEr", "LogisticNUCLEAR", "LogisticLASSO", "LogisticSpINNEr"];
style.Labels = {'SpINNEr', 'Logistic NUCLEAR', 'Logistic LASSO', 'Logistic SpINNEr'};
style.Colors = [0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560; ...
    0.9290 0.6940 0.1250; 0 0.4470 0.7410];
style.Markers = {'o', 's', 'd', '^'};

%% Fig. 2: intact blocks.
intact = [edgeIntact; blockIntact];
intact.PlotX = zeros(height(intact), 1);
for b = 1:numel(blockOrder)
    intact.PlotX(intact.BlockLabel == blockOrder(b)) = b;
end
assert(all(intact.PlotX > 0) && height(intact) == 4*numel(blockOrder), ...
    'Unexpected intact block configurations.');
fig = newFigure(12, 5); % Extra width keeps seven block labels readable.
ax = axes('Parent', fig);
drawSeries(ax, intact, 'PlotX', style);
formatAxes(ax);
% Printed at about 29% (one column): thicker curves and larger text than the
% shared style. The tick font is limited by the seven block labels fitting side by side.
set(findobj(ax, 'Type', 'ErrorBar'), 'LineWidth', 3.2, 'MarkerSize', 14, 'CapSize', 14);
set(ax, 'FontSize', 20, 'LineWidth', 2.4);
ax.XRuler.TickLabelGapMultiplier = 0.2;
xlim(ax, [0.65 numel(blockOrder) + 0.35]);
xticks(ax, 1:numel(blockOrder));
xticklabels(ax, "(" + replace(blockOrder, "-", ",") + ")");
ax.XTickLabelRotation = 0;   % MATLAB would tilt the labels at this font size
ylim(ax, [0.80 1.00]);
yticks(ax, 0.80:0.05:1.00);
xlabel(ax, 'Block sizes', 'FontSize', 23, 'FontWeight', 'bold');
ylabel(ax, 'Relative Frobenius Error', 'FontSize', 23, 'FontWeight', 'bold');
legend(ax, 'Location', 'southeast', 'FontSize', 22);
saveFigure(fig, figureDir, 'fig2-intact-block-gamma05');

%% Fig. 3(a): edge relocation, blocks (5,5,4).
fig = newFigure(8.5, 5);
relocationPanel(axes('Parent', fig), edge, 'SwapProportion', ...
    [0.75 1.00], 0.75:0.05:1.00, style);
saveFigure(fig, figureDir, 'fig3-relocation-gamma05_a');

%% Fig. 3(b)-(f): block relocation ((5,5,5) is not shown in the manuscript).
labels = ["6-6-6", "7-7-7", "8-8-8", "12-12-12", "16-16-16"];
panels = 'bcdef';
yLimits = [0.86 1.00; 0.90 1.00; 0.92 1.00; 0.95 1.00; 0.96 1.00];
ySteps = [0.02 0.02 0.02 0.01 0.01];
for b = 1:numel(labels)
    fig = newFigure(8.5, 5);
    current = blocks(blocks.BlockLabel == labels(b), :);
    relocationPanel(axes('Parent', fig), current, 'ActualSwapProportion', ...
        yLimits(b, :), yLimits(b, 1):ySteps(b):yLimits(b, 2), style);
    saveFigure(fig, figureDir, sprintf('fig3-relocation-gamma05_%s', panels(b)));
end

%% Fig. 3, one-row version: four square panels without legends plus a shared legend strip.
rowLabels = ["6-6-6", "7-7-7", "12-12-12"];
rowLimits = [0.80 1.00; 0.88 1.00; 0.90 1.00; 0.95 1.00];
rowSteps = [0.05 0.02 0.02 0.01];
rowData = {edge, blocks(blocks.BlockLabel == rowLabels(1), :), ...
    blocks(blocks.BlockLabel == rowLabels(2), :), blocks(blocks.BlockLabel == rowLabels(3), :)};
rowX = {'SwapProportion', 'ActualSwapProportion', 'ActualSwapProportion', 'ActualSwapProportion'};
for b = 1:4
    fig = newFigure(4.2, 4.2);
    ax = axes('Parent', fig);
    relocationPanel(ax, rowData{b}, rowX{b}, rowLimits(b, :), ...
        rowLimits(b, 1):rowSteps(b):rowLimits(b, 2), style);
    legend(ax, 'off');
    % Smaller tick numbers (17 -> 14) so the panels look less packed; axis titles stay at 12.
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    ax.XLabel.FontSize = 12;
    ax.YLabel.FontSize = 12;
    saveFigure(fig, figureDir, sprintf('fig3-row-gamma05_%s', panels(b) - 1));
end
fig = newFigure(16.8, 0.8);
ax = axes('Parent', fig, 'Visible', 'off');
hold(ax, 'on');
for m = numel(style.Keys):-1:1   % Logistic SpINNEr first
    label = style.Labels{m};
    if m > 1
        label = [label, repmat(char(160), 1, 6)];   % extra space between items (plain spaces are trimmed)
    end
    errorbar(ax, NaN, NaN, NaN, 'Color', style.Colors(m, :), 'LineWidth', 1.9, ...
        'Marker', style.Markers{m}, 'MarkerFaceColor', style.Colors(m, :), ...
        'MarkerSize', 9, 'CapSize', 10, 'DisplayName', label);
end
legend(ax, 'Orientation', 'horizontal', 'FontSize', 17, 'Location', 'north', 'LineWidth', 1.2);
saveFigure(fig, figureDir, 'fig3-row-gamma05_legend');

%% Local helpers.
function data = readSummary(path, isBlock)
    assert(isfile(path), 'Missing source data: %s', path);
    opts = detectImportOptions(path, 'TextType', 'string');
    opts = setvartype(opts, 'Method', 'string');
    if isBlock
        opts = setvartype(opts, 'BlockLabel', 'string');
    end
    data = readtable(path, opts);
    required = {'Method', 'MeanACC', 'TwoSEM'};
    if isBlock
        required = [required, {'BlockLabel', 'ActualSwapProportion'}];
    else
        required = [required, {'SwapProportion'}];
    end
    assert(all(ismember(required, data.Properties.VariableNames)), ...
        'Missing required columns in %s.', path);
    assert(all(isfinite(data.MeanACC)) && all(isfinite(data.TwoSEM)) ...
        && all(data.TwoSEM >= 0), 'Invalid means or error bars in %s.', path);
end

function fig = newFigure(width, heightInches)
    fig = figure('Units', 'inches', 'Position', [1 1 width heightInches], ...
        'Color', 'w', 'Visible', 'off');
end

function relocationPanel(ax, data, xName, yLimits, yTicks, style)
    drawSeries(ax, data, xName, style);
    formatAxes(ax);
    xlim(ax, [-0.01 0.245]);   % the largest relocated proportion is 0.238, for (7,7,7)
    xticks(ax, [0 0.1 0.2]);
    xticklabels(ax, {'0', '0.1', '0.2'});
    ylim(ax, yLimits);
    yticks(ax, yTicks);
    xlabel(ax, 'Relocated edge proportion', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'Relative Frobenius Error', 'FontSize', 12, 'FontWeight', 'bold');
    legend(ax, 'Location', 'southeast', 'FontSize', 17);
end

function drawSeries(ax, data, xName, style)
    assert(all(ismember(data.Method, style.Keys)), 'Unrecognized method in input.');
    hold(ax, 'on');
    referenceX = [];
    for m = 1:numel(style.Keys)
        current = sortrows(data(data.Method == style.Keys(m), :), xName);
        assert(~isempty(current), 'Missing method: %s.', style.Keys(m));
        x = current.(xName);
        assert(all(isfinite(x)) && numel(unique(x)) == numel(x), ...
            'Invalid or duplicate x values for %s.', style.Keys(m));
        if m == 1
            referenceX = x;
        else
            assert(isequal(x, referenceX), 'Methods have different x grids.');
        end
        errorbar(ax, x, current.MeanACC, current.TwoSEM, ...
            'Color', style.Colors(m, :), 'LineWidth', 1.9, ...
            'Marker', style.Markers{m}, 'MarkerFaceColor', style.Colors(m, :), ...
            'MarkerSize', 9, 'CapSize', 10, 'DisplayName', style.Labels{m});
    end
end

function formatAxes(ax)
    grid(ax, 'on');
    set(ax, 'FontName', 'Helvetica', 'FontSize', 17, ...
        'LineWidth', 1.7, 'TickDir', 'in', 'Box', 'on');
    % Match the close tick-label spacing used by the reference MATLAB script.
    if isprop(ax.XRuler, 'TickLabelGapMultiplier')
        ax.XRuler.TickLabelGapMultiplier = -0.1;
    end
    if isprop(ax.YRuler, 'TickLabelGapMultiplier')
        ax.YRuler.TickLabelGapMultiplier = 0;
    end
end

function saveFigure(fig, directory, stem)
    cleanup = onCleanup(@() close(fig));
    drawnow;
    path = fullfile(directory, [stem '.pdf']);
    exportgraphics(fig, path, 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('Created: %s\n', path);
end
