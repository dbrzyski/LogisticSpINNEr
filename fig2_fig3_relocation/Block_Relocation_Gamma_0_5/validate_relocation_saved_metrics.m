function validate_relocation_saved_metrics(result, config)
%VALIDATE_RELOCATION_SAVED_METRICS Validate completed fits, including on resume.

expectedSize = [config.nConditions, config.nMethods];
if ~isequal(size(result.completed), expectedSize) || ...
        ~isequal(size(result.AUC), expectedSize) || ...
        ~isequal(size(result.Euclid), expectedSize) || ...
        ~isequal(size(result.Probs), [config.n, expectedSize]) || ...
        ~isequal(size(result.betaHat), expectedSize)
    error('BlockRelocation:InvalidSavedMetricShape', ...
        'Checkpoint completion flags or prediction metrics have incorrect dimensions.');
end
for methodIndex = 1:config.nMethods
    for conditionIndex = 1:config.nConditions
        probs = result.Probs(:, conditionIndex, methodIndex);
        auc = result.AUC(conditionIndex, methodIndex);
        euclid = result.Euclid(conditionIndex, methodIndex);
        if ~ismember(methodIndex, config.predictionMethodIndices)
            if ~isnan(auc) || ~isnan(euclid) || any(~isnan(probs))
                error('BlockRelocation:UnexpectedGaussianMetrics', ...
                    'Gaussian SpINNEr must not have probability metrics.');
            end
        elseif result.completed(conditionIndex, methodIndex)
            estimate = struct('AUC', auc, 'Euclid', euclid, 'Probs', probs);
            read_relocation_prediction_metrics( ...
                estimate, result.responses(:, conditionIndex));
        end
        if result.completed(conditionIndex, methodIndex) && ...
                (~isfinite(result.ACC(conditionIndex, methodIndex)) || ...
                ~isfinite(result.betaHat(conditionIndex, methodIndex)))
            error('BlockRelocation:InvalidCompletedFit', ...
                'Completed fits must contain finite RFE and intercept estimates.');
        end
    end
end
end
