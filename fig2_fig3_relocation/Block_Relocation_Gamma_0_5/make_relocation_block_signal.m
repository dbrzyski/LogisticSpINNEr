function [B, info] = make_relocation_block_signal(p, blockSizes, coefficientMagnitude)
%MAKE_RELOCATION_BLOCK_SIGNAL Reuse the original node pools at a fixed magnitude.
%   The original fixed node order makes larger blocks contain smaller ones.
%   Only coefficient magnitude changes; no normalization by block size is used.

if ~isscalar(coefficientMagnitude) || ~isfinite(coefficientMagnitude) || ...
        coefficientMagnitude <= 0
    error('BlockRelocation:InvalidMagnitude', ...
        'coefficientMagnitude must be a positive finite scalar.');
end

if p ~= 60
    error('BlockRelocation:InvalidDimension', ...
        'This experiment is fixed at p=60.');
end
if ~isnumeric(blockSizes) || ~isequal(size(blockSizes), [1, 3]) || ...
        any(blockSizes < [5, 5, 4]) || any(blockSizes > 20) || ...
        any(blockSizes ~= floor(blockSizes))
    error('BlockRelocation:InvalidBlockSizes', ...
        'blockSizes must be a 1-by-3 integer vector between [5,5,4] and 20.');
end

nodePools = { ...
    [6:10, 1:5, 11:15, 33:37], ...
    [17:21, 16, 22:28, 38:44], ...
    [29:32, 45:60]};
blockSigns = [1, -1, 1];

allNodes = [nodePools{:}];
if numel(unique(allNodes)) ~= p || ~isequal(sort(allNodes), 1:p)
    error('BlockRelocation:InvalidNodePools', ...
        'The fixed node pools must form a partition of 1:p.');
end

B = zeros(p, p);
blockNodes = cell(3, 1);
for blockIndex = 1:3
    currentNodes = nodePools{blockIndex}(1:blockSizes(blockIndex));
    blockNodes{blockIndex} = currentNodes(:);
    B(currentNodes, currentNodes) = coefficientMagnitude * blockSigns(blockIndex);
end

nonzeroValues = B(B ~= 0);
if any(abs(nonzeroValues) ~= coefficientMagnitude)
    error('BlockRelocation:CoefficientMismatch', ...
        'Every nonzero coefficient must have the configured magnitude.');
end

upperMask = triu(true(p), 1);
nonzeroUpperEdgeCount = nnz(upperMask & (B ~= 0));
info.blockSizes = blockSizes;
info.blockSigns = blockSigns;
info.nodePools = nodePools;
info.blockNodes = blockNodes;
info.permutedNodeOrder = allNodes;
info.coefficientMagnitude = coefficientMagnitude;
info.nonzeroUpperEdgeCount = nonzeroUpperEdgeCount;
info.totalUpperEdgeCount = p * (p - 1) / 2;
info.upperEdgeDensity = nonzeroUpperEdgeCount / info.totalUpperEdgeCount;
info.fullMatrixNonzeroCount = nnz(B);
info.fullMatrixDensity = nnz(B) / numel(B);
info.upperSignalEnergy = sum(B(upperMask).^2);
info.rank = rank(B);
info.frobeniusNorm = norm(B, 'fro');
info.nuclearNorm = sum(svd(B));
info.l1Norm = sum(abs(B(:)));
end
