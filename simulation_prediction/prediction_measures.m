function meas = prediction_measures(est, data, method, cfg)
%PREDICTION_MEASURES Test-set prediction measures of one fitted method.
%   The fitted linear predictor etaHat = <A, B> + beta is turned into a probability
%   pHat = logistic(etaHat) for the logistic methods. SpINNEr's least-squares score is on
%   the probability scale, so it is recalibrated by Platt scaling: a logistic regression
%   of the responses on the scores of an independent calibration sample, then applied to
%   the test scores. (Training scores cannot be used: they often separate the training
%   responses perfectly, and the Platt coefficients diverge.)
%   A score with SD below 1e-8 is treated as constant (empty model): AUC 0.5, corrEta 0,
%   calSlope NaN.
%
%   Measures (true probabilities p = logistic(eta) are known in simulation):
%     AUC            test AUC of the score (rank-based, no calibration needed)
%     oracleAUC      test AUC of the true eta (best achievable)
%     excessKL       mean KL(p || pHat) = expected test log-loss minus Bayes log-loss
%     logLoss        test log-loss on the observed y;  bayesLogLoss for the true p
%     brier          test Brier score on y;  bayesBrier for the true p
%     probMSE        mean (pHat - p)^2
%     corrEta        Pearson correlation of the score with the true eta
%     calSlope       slope of the logistic regression of yTest on logit(pHat) (1 = calibrated)

p = cfg.p;
tol = 1e-8;
scoreCal  = reshape(data.Acal, p^2, [])' * est.B(:) + est.beta(1);
scoreTest = reshape(data.Atest, p^2, [])' * est.B(:) + est.beta(1);
isConstant = std(scoreTest) < tol;

meas.plattCoef = [NaN; NaN];
if strcmp(method, 'SpINNEr')
    if ~isConstant
        meas.plattCoef = glmfit(scoreCal, data.yCal, 'binomial');
        etaHat = meas.plattCoef(1) + meas.plattCoef(2) * scoreTest;
    else
        etaHat = repmat(log(mean(data.yCal) / (1 - mean(data.yCal))), cfg.nTest, 1);
    end
else
    etaHat = scoreTest;
end

clip = 1e-12;
pHat  = min(max(1 ./ (1 + exp(-etaHat)), clip), 1 - clip);
pTrue = 1 ./ (1 + exp(-data.etaTest));
y = data.yTest;

if isConstant
    meas.AUC = 0.5;
else
    [~, ~, ~, meas.AUC] = perfcurve(y, scoreTest, 1);
end
[~, ~, ~, meas.oracleAUC] = perfcurve(y, data.etaTest, 1);
meas.excessKL  = mean(pTrue .* log(pTrue ./ pHat) + (1 - pTrue) .* log((1 - pTrue) ./ (1 - pHat)));
meas.logLoss      = -mean(y .* log(pHat) + (1 - y) .* log(1 - pHat));
meas.bayesLogLoss = -mean(y .* log(pTrue) + (1 - y) .* log(1 - pTrue));
meas.brier      = mean((pHat - y).^2);
meas.bayesBrier = mean((pTrue - y).^2);
meas.probMSE    = mean((pHat - pTrue).^2);
if isConstant
    meas.corrEta  = 0;          % no association; keeps empty fits in the averages
    meas.calSlope = NaN;
else
    meas.corrEta  = corr(scoreTest, data.etaTest);
    coef          = glmfit(etaHat, y, 'binomial');
    meas.calSlope = coef(2);
end
meas.relFrobError = norm(est.B - cfg.Btrue, 'fro') / norm(cfg.Btrue, 'fro');
end
