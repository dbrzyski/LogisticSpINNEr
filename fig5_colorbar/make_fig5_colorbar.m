function make_fig5_colorbar()
%MAKE_FIG5_COLORBAR Readable colorbar for the glass-brain figure (Fig. 5).
%   The BrainNet Viewer rendering in glass_brain_final.pdf colors each edge by its median bootstrap
%   coefficient on the jet colormap, from the smallest (-0.1210) to the largest (0.1433) coefficient
%   of the 41 cortical edges; its own colorbar is unreadable at print size. This script recomputes
%   that range from the bootstrap estimates and draws the same colormap as a vector PDF with larger
%   labels. The manuscript crops the original colorbar and places this one next to the brains.
%
%   Output: fig5-colorbar.pdf in this folder.

here = fileparts(mfilename('fullpath'));
bootFile = fullfile(here, '..', 'fig5_fig6_real_data', 'results', 'EstimatesBootstrap.mat');
S = load(bootFile, 'estimAll');
E = S.estimAll;                                  % 214 x 214 x 300

% Edges of Figs. 5-6: nonzero in >= 90% of the fits, both regions cortical (1-200)
freq = mean(E ~= 0, 3);
coef = median(E, 3);
sel = triu(freq >= 0.90, 1);
sel(201:end, :) = false;
sel(:, 201:end) = false;
vals = coef(sel);
lims = [min(vals) max(vals)];
fprintf('%d edges, median coefficients from %.4f to %.4f\n', nnz(sel), lims);
assert(nnz(sel) == 41, 'Expected the 41 cortical edges of Figs. 5-6.');
assert(all(abs(round(lims, 4) - [-0.1210 0.1433]) < 1e-12), ...
    'Range differs from the colorbar printed in glass_brain_final.pdf.');

fig = figure('Units', 'inches', 'Position', [1 1 0.9 1.6], 'Color', 'w', 'Visible', 'off');
ax = axes('Parent', fig, 'Units', 'inches', 'Position', [0.08 0.08 0.16 1.44]);
imagesc(ax, [0 1], lims, linspace(lims(1), lims(2), 256)');
set(ax, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'on', ...
    'FontName', 'Helvetica', 'FontSize', 9, 'LineWidth', 0.5, 'TickLength', [0.04 0.04]);
colormap(ax, jet(256));
yticks(ax, [-0.1 -0.05 0 0.05 0.1]);
yticklabels(ax, {'-0.10', '-0.05', '0', '0.05', '0.10'});
ylabel(ax, 'Median coefficient', 'FontSize', 9);

out = fullfile(here, 'fig5-colorbar.pdf');
exportgraphics(fig, out, 'ContentType', 'vector', 'BackgroundColor', 'white');
close(fig);
fprintf('Created: %s\n', out);
end
