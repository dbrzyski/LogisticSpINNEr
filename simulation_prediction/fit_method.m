function estimate = fit_method(y, X, A, settings)
%FIT_METHOD Fit one method through the spinner class and return its Estimate.
%   spinner.fit stores the estimate and then appends
%   it to History (a cell array); some MATLAB releases reject that concatenation. A
%   completed estimate is used in that case; any earlier error is rethrown.

model = spinner(y, X, A);
try
    model.fit(settings);
catch fitError
    est = model.Estimate;
    if ~(isstruct(est) && isfield(est, 'B') && all(isfinite(est.B(:))))
        rethrow(fitError);
    end
    warning('SimPrediction:HistoryUpdate', ...
        'Fit completed but the History update failed; using the completed Estimate.');
end
estimate = model.Estimate;
end
