function plot_solver_timing_paper(smokeTest)
%PLOT_SOLVER_TIMING_PAPER One-column manuscript figure of the single-thread solver timing study.
%   Two separate 1.72 x 2.1 in panels, placed side by side as subfigures in the manuscript:
%   (a) median fit time against p, one colour and marker per n, one line style per solver;
%   (b) speed-up of Woodbury over direct against p (geometric mean over data sets).
%   The SVD solver is timed by run_solver_timing but not plotted (not in the manuscript).
%   Reads results/time_*.mat; writes
%   results/solver_timing_paper_a.pdf and _b.pdf. The manuscript's figure files are not touched.

if nargin < 1
    smokeTest = false;
end
cfg = config_solver_timing(smokeTest);
files = dir(fullfile(cfg.resultsDir, 'time_p*_n*_rep*.mat'));
res = arrayfun(@(f) load(fullfile(f.folder, f.name), 'result').result, files);
pAll = [res.p]'; nAll = [res.n]';
sec = [res.seconds];
tW = [sec.woodbury]'; tD = [sec.direct]';
pVals = unique(pAll)'; nVals = unique(nAll)';

% Per (n, p): median times and geometric-mean speed-ups
gm = @(x) exp(mean(log(x)));
medW = zeros(numel(nVals), numel(pVals)); medD = medW; gmDW = medW;
for i = 1:numel(nVals)
    for j = 1:numel(pVals)
        k = nAll == nVals(i) & pAll == pVals(j);
        medW(i, j) = median(tW(k));
        medD(i, j) = median(tD(k));
        gmDW(i, j) = gm(tD(k) ./ tW(k));
    end
end

tickFont  = 7;
labelFont = 6.5;
lineWidth = 0.8;
axisWidth = 0.7;
markSize  = 3.6;
colors  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
markers = {'o', 's', 'd'};
panelSize = [1.72, 2.1];

% (a) fit times; legend entries per (n, solver)
[fig1, layout1] = new_panel();
ax1 = nexttile(layout1); hold(ax1, 'on');
h = gobjects(2*numel(nVals), 1);
for i = 1:numel(nVals)
    h(2*i-1) = plot(ax1, pVals, medW(i, :), '-', 'Color', colors(i, :), 'LineWidth', lineWidth, ...
        'Marker', markers{i}, 'MarkerSize', markSize, 'MarkerFaceColor', colors(i, :), ...
        'DisplayName', sprintf('n = %d, Alg. 1', nVals(i)));
    h(2*i) = plot(ax1, pVals, medD(i, :), '--', 'Color', colors(i, :), 'LineWidth', lineWidth, ...
        'Marker', markers{i}, 'MarkerSize', markSize, 'MarkerFaceColor', 'w', ...
        'DisplayName', sprintf('n = %d, direct', nVals(i)));
end
yLim1 = [1e-2, 1e3];
style_axes(ax1, pVals, 10.^(-2:3), yLim1);
xlabel(ax1, 'p', 'FontSize', labelFont, 'FontWeight', 'bold');
ylabel(ax1, 'Median fit time [s]', 'FontSize', labelFont, 'FontWeight', 'bold');
lgd1 = legend(ax1, h, 'FontSize', 5.5, 'LineWidth', axisWidth);
lgd1.ItemTokenSize = [9, 4];
% right-hand side, in the gap between the Algorithm 1 and the direct curves
lgd1.Units = 'normalized';
axPos = ax1.Position;
yFrac = @(v) (log10(v) - log10(yLim1(1))) / diff(log10(yLim1));   % data value -> fraction of axis height
lgd1.Position(1) = axPos(1) + axPos(3) - lgd1.Position(3) - 0.035;
lgd1.Position(2) = axPos(2) + axPos(4) * (yFrac(max(medW(:))) + 0.04);
save_panel(fig1, 'solver_timing_paper_a.pdf');

% (b) speed-up of Woodbury over direct
[fig2, layout2] = new_panel();
ax2 = nexttile(layout2); hold(ax2, 'on');
for i = 1:numel(nVals)
    plot(ax2, pVals, gmDW(i, :), '-', 'Color', colors(i, :), 'LineWidth', lineWidth, ...
        'Marker', markers{i}, 'MarkerSize', markSize, 'MarkerFaceColor', colors(i, :), ...
        'DisplayName', sprintf('n = %d', nVals(i)));
end
style_axes(ax2, pVals, 10.^(1:4), [10, 1e4]);
xlabel(ax2, 'p', 'FontSize', labelFont, 'FontWeight', 'bold');
ylabel(ax2, 'Speed-up', 'FontSize', labelFont, 'FontWeight', 'bold');
lgd2 = legend(ax2, 'Location', 'northwest', 'FontSize', 6, 'LineWidth', axisWidth);
lgd2.ItemTokenSize = [10, 4];
save_panel(fig2, 'solver_timing_paper_b.pdf');

    function [fig, layout] = new_panel()
        fig = figure('Units', 'inches', 'Position', [1 1 panelSize], 'Color', 'w');
        layout = tiledlayout(fig, 1, 1, 'Padding', 'tight');   % minimal outer whitespace
    end

    function save_panel(fig, fileName)
        % Exact panelSize page (exportgraphics would crop to the content)
        set(fig, 'PaperUnits', 'inches', 'PaperSize', panelSize, 'PaperPosition', [0, 0, panelSize]);
        print(fig, fullfile(cfg.resultsDir, fileName), '-dpdf', '-vector');
    end

    function style_axes(ax, xTicks, yTicks, yLim)
        xLim = [xTicks(1) - 5, xTicks(end) + 5];
        if isscalar(xTicks)   % smoke test
            xLim = xTicks + [-5, 5];
        end
        set(ax, 'YScale', 'log', 'XTick', xTicks, 'XLim', xLim, ...
            'YTick', yTicks, 'YLim', yLim, 'FontSize', tickFont, 'LineWidth', axisWidth, ...
            'TickDir', 'in', 'Box', 'on', 'YMinorTick', 'off');
        grid(ax, 'on');
        ax.YMinorGrid = 'off';
        ax.XRuler.TickLabelGapMultiplier = 0.2;
        ax.YRuler.TickLabelGapMultiplier = 0.2;
    end
end
