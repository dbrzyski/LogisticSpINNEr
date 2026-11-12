function estimate = fit_spinner_method(y, X, A, settings)
%FIT_SPINNER_METHOD Fit one method through the current SpINNEr class API.

model = spinner(y, X, A);

try
    model.fit(settings);
catch fitError
    % SpINNEr-main stores the completed estimate immediately before its
    % History update. The current History property is a cell array, whereas
    % the completed estimate is a structure. Some MATLAB releases reject
    % that mixed-type concatenation. Continue only when the solver has
    % already returned a complete estimate; all earlier errors are rethrown.
    if ~is_complete_estimate(model.Estimate, size(A, 1))
        rethrow(fitError);
    end
    warning('SwapPilot:SpinnerHistoryUpdate', ...
        ['SpINNEr completed the fit, but its History update failed. ', ...
        'The completed Estimate is being used without changing SpINNEr-main.']);
end

estimate = model.Estimate;
if ~is_complete_estimate(estimate, size(A, 1))
    error('SwapPilot:IncompleteEstimate', ...
        'SpINNEr did not return a valid coefficient matrix.');
end
end

function isComplete = is_complete_estimate(estimate, p)
isComplete = isstruct(estimate) && isfield(estimate, 'B') && ...
    isequal(size(estimate.B), [p, p]) && all(isfinite(estimate.B(:)));
end
