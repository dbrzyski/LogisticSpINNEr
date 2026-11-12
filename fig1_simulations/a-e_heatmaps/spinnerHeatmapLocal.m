function spinnerHeatmapLocal(b, tit, bottom, top, colMap)
if nargin == 1
    tit ='';
end

if nargin <= 2
    b       = double(b);
    bottom  = min(min(b));
    top     = max(max(b));
end

if nargin <= 4
    generateColorMap = true;
else
    generateColorMap = false;
end

%% heatmap
margins = [0.05, 0.05];
figure('name','Estimates');
axes('Position', [margins(2), margins(1), (1-2*margins(2)), (1-2*margins(1))]);
imagesc(b);
if generateColorMap
    colormap(spinnerColormap(bottom ,top));
else
    colormap(colMap);
end
title(tit);
axis equal;
axis tight;
clim manual
clim([bottom  top]);
set(gca, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []);
cb = colorbar('Location','northoutside'); 
cb.FontSize = 14; 
cminmax = clim;                 % [vmin vmax]
cb.Ticks = linspace(cminmax(1), cminmax(2), 4);   % min, mid, max 

% Round to 2 digits after comma
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);

% --- Make positions stable and current before editing ---
ax = gca;
drawnow;                              % let MATLAB settle layout
ax.Units = 'normalized';
cb.Units = 'normalized';

% --- Save axes position to lock layout ---
axPos = ax.Position;                  % [x y w h] of the axes

% --- Shorten the colorbar horizontally and place it just above the axes ---
frac = 0.55;                          % keep 60% of the axes width (tune as needed)
gap  = 0.010;                         % vertical gap above the axes (normalized)

cbPos      = cb.Position;             % start from current [x y w h]
cbPos(3)   = frac * axPos(3);         % width = fraction of axes width
cbPos(1)   = axPos(1) + (axPos(3)-cbPos(3))/2; % center under axes
cbPos(4)   = cbPos(4);                % keep the current height
cbPos(2)   = axPos(2) + axPos(4) + gap;        % place above axes, not inside
cb.Position = cbPos;

% --- Restore axes position in case MATLAB nudged it during the change ---
ax.Position = axPos;


end
