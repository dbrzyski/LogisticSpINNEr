function [responses, responseUniforms, probabilities] = generate_relocation_paired_responses( ...
        A, X, BTrue, intercept, responseSeed)
%GENERATE_RELOCATION_PAIRED_RESPONSES Reuse one Bernoulli uniform across conditions.

[p1, p2, n] = size(A);
if p1 ~= p2 || size(BTrue, 1) ~= p1 || size(BTrue, 2) ~= p1
    error('BlockRelocation:DimensionMismatch', ...
        'A and BTrue must contain conformable square matrices.');
end
if size(X, 1) ~= n
    error('BlockRelocation:DimensionMismatch', ...
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
        if eta >= 0
            probability = 1 / (1 + exp(-eta));
        else
            expEta = exp(eta);
            probability = expEta / (1 + expEta);
        end
        probabilities(observationIndex, conditionIndex) = probability;
        responses(observationIndex, conditionIndex) = ...
            responseUniforms(observationIndex) < probability;
    end
end
end
