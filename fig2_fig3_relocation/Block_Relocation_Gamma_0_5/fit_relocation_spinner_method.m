function estimate = fit_relocation_spinner_method(y, X, A, settings)
%FIT_RELOCATION_SPINNER_METHOD Fit one method through the current spinner class.

model = spinner(y, X, A);
try
    model.fit(settings);
catch fitError
    % A completed estimate may exist before the known History type mismatch.
    % Continue only in that narrow case; every solver/CV error is rethrown.
    if ~is_complete_estimate(model.Estimate, size(A, 1))
        rethrow(fitError);
    end
    warning('BlockRelocation:SpinnerHistoryUpdate', ...
        ['SpINNEr completed the fit, but its History update failed. ', ...
        'The completed Estimate is used without changing SpINNEr-main.']);
end

estimate = model.Estimate;
if ~is_complete_estimate(estimate, size(A, 1))
    error('BlockRelocation:IncompleteEstimate', ...
        'SpINNEr did not return a valid coefficient matrix.');
end
end

function isComplete = is_complete_estimate(estimate, p)
isComplete = isstruct(estimate) && isfield(estimate, 'B') && ...
    isequal(size(estimate.B), [p, p]) && all(isfinite(estimate.B(:)));
end
