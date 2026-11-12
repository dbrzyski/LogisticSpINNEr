function [BTrue, swapInfo] = make_nested_swapped_signals(BBaseline, kGrid, swapSeed)
%MAKE_NESTED_SWAPPED_SIGNALS Move nested sets of upper-triangular signal edges.

if ~ismatrix(BBaseline) || size(BBaseline, 1) ~= size(BBaseline, 2)
    error('SwapPilot:InvalidSignal', 'BBaseline must be square.');
end
if norm(BBaseline - BBaseline.', 'fro') > 1e-12
    error('SwapPilot:InvalidSignal', 'BBaseline must be symmetric.');
end
if isempty(kGrid) || any(kGrid < 0) || any(kGrid ~= floor(kGrid)) || ...
        any(diff(kGrid) < 0)
    error('SwapPilot:InvalidKGrid', ...
        'kGrid must contain nonnegative integers in nondecreasing order.');
end

p = size(BBaseline, 1);
upperMask = triu(true(p), 1);
sourceCandidates = find(upperMask & (BBaseline ~= 0));
targetCandidates = find(upperMask & (BBaseline == 0));

if max(kGrid) > min(numel(sourceCandidates), numel(targetCandidates))
    error('SwapPilot:InvalidKGrid', ...
        'The requested number of swaps exceeds the available edge count.');
end

previousState = rng;
restoreState = onCleanup(@() rng(previousState)); %#ok<NASGU>
rng(swapSeed, 'twister');
sourceOrder = sourceCandidates(randperm(numel(sourceCandidates)));
targetOrder = targetCandidates(randperm(numel(targetCandidates)));
sourceOrder = sourceOrder(:);
targetOrder = targetOrder(:);

[sourceRows, sourceColumns] = ind2sub([p, p], sourceOrder);
[targetRows, targetColumns] = ind2sub([p, p], targetOrder);
sourceValues = BBaseline(sourceOrder);

nConditions = numel(kGrid);
BTrue = zeros(p, p, nConditions);
removedEdges = cell(nConditions, 1);
addedEdges = cell(nConditions, 1);

for conditionIndex = 1:nConditions
    k = kGrid(conditionIndex);
    currentSignal = BBaseline;

    for edgeIndex = 1:k
        sourceRow = sourceRows(edgeIndex);
        sourceColumn = sourceColumns(edgeIndex);
        targetRow = targetRows(edgeIndex);
        targetColumn = targetColumns(edgeIndex);
        edgeValue = sourceValues(edgeIndex);

        currentSignal(sourceRow, sourceColumn) = 0;
        currentSignal(sourceColumn, sourceRow) = 0;
        currentSignal(targetRow, targetColumn) = edgeValue;
        currentSignal(targetColumn, targetRow) = edgeValue;
    end

    BTrue(:, :, conditionIndex) = currentSignal;
    removedEdges{conditionIndex} = [ ...
        sourceRows(1:k), sourceColumns(1:k), sourceValues(1:k)];
    addedEdges{conditionIndex} = [ ...
        targetRows(1:k), targetColumns(1:k), sourceValues(1:k)];
end

swapInfo.swapSeed = swapSeed;
swapInfo.kGrid = kGrid;
swapInfo.swapProportions = kGrid ./ numel(sourceCandidates);
swapInfo.sourceCandidateCount = numel(sourceCandidates);
swapInfo.targetCandidateCount = numel(targetCandidates);
swapInfo.sourceEdgeOrder = [sourceRows, sourceColumns, sourceValues];
swapInfo.targetEdgeOrder = [targetRows, targetColumns];
swapInfo.removedEdges = removedEdges;
swapInfo.addedEdges = addedEdges;
end
