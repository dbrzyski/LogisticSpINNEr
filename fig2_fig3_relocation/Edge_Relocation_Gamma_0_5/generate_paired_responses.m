function [responses, responseUniforms, probabilities] = generate_paired_responses( ...
    A, X, BTrue, intercept, responseSeed)
%GENERATE_PAIRED_RESPONSES Use one uniform draw per subject across all k values.

[p1, p2, n] = size(A);
if p1 ~= p2 || size(BTrue, 1) ~= p1 || size(BTrue, 2) ~= p1
    error('SwapPilot:DimensionMismatch', ...
        'A and BTrue must contain conformable square matrices.');
end
if size(X, 1) ~= n
    error('SwapPilot:DimensionMismatch', ...
        'X and A must contain the same number of observations.');
end

previousState = rng;
restoreState = onCleanup(@() rng(previousState)); %#ok<NASGU>
rng(responseSeed, 'twister');
responseUniforms = rand(n, 1);

nConditions = size(BTrue, 3);
responses = zeros(n, nConditions);
probabilities = zeros(n, nConditions);

for conditionIndex = 1:nConditions
    currentSignal = BTrue(:, :, conditionIndex);
    for observationIndex = 1:n
        eta = trace(A(:, :, observationIndex)' * currentSignal) + ...
            X(observationIndex, :) * intercept;
        probability = 1 / (1 + exp(-eta));
        probabilities(observationIndex, conditionIndex) = probability;
        responses(observationIndex, conditionIndex) = ...
            responseUniforms(observationIndex) < probability;
    end
end
end
