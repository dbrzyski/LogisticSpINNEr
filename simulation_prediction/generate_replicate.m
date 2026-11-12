function data = generate_replicate(cfg, r)
%GENERATE_REPLICATE Training and test data for replicate r (deterministic in cfg.seed and r).
%   Training data are generated as in Fig. 1: symmetric Gaussian matrices, each entry
%   z-scored across the training subjects, zero diagonals, y ~ Bernoulli(logistic(eta)).
%   Test matrices come from a separate random stream and are standardized with the
%   training means and SDs; their responses follow the same model. A calibration sample
%   (size n, third stream) is generated the same way; it is used only to Platt-scale
%   SpINNEr's least-squares score (see prediction_measures.m).

trainStream = RandStream('twister', 'Seed', cfg.seed + r);
testStream  = RandStream('twister', 'Seed', cfg.seed + 100000 + r);
calStream   = RandStream('twister', 'Seed', cfg.seed + 200000 + r);

[Araw, yUnif]         = draw(trainStream, cfg.p, cfg.n);
[AtestRaw, yTestUnif] = draw(testStream, cfg.p, cfg.nTest);
[AcalRaw, yCalUnif]   = draw(calStream, cfg.p, cfg.n);

% Standardization with training statistics (zscore of the training set, as in Fig. 1)
p = cfg.p;
Xtrain = reshape(Araw, p^2, cfg.n);
centre = mean(Xtrain, 2);
scale  = std(Xtrain, 0, 2);
isDiag = logical(reshape(eye(p), p^2, 1));
standardize = @(Araw) reshape((reshape(Araw, p^2, []) - centre) ./ scale .* ~isDiag, ...
    p, p, []);

data.A     = standardize(Araw);
data.Atest = standardize(AtestRaw);
data.Acal  = standardize(AcalRaw);
data.X     = ones(cfg.n, 1);
linPred = @(A) reshape(A, p^2, [])' * cfg.Btrue(:) + cfg.betaTrue;
data.eta     = linPred(data.A);
data.etaTest = linPred(data.Atest);
etaCal       = linPred(data.Acal);
data.y      = double(yUnif     < 1 ./ (1 + exp(-data.eta)));
data.yTest  = double(yTestUnif < 1 ./ (1 + exp(-data.etaTest)));
data.yCal   = double(yCalUnif  < 1 ./ (1 + exp(-etaCal)));
end

function [A, u] = draw(stream, p, n)
A = randn(stream, [p, p, n]);
A = A + permute(A, [2 1 3]);
u = rand(stream, n, 1);
end
