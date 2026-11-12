%% Clear the workspace
clear; close all; clc;

%% Paths
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'software_84509f9'))   % SpINNEr toolbox at commit 84509f9, used for Fig. 1
SavePath = fullfile(here, 'output');   % B_ESTIMS.mat and the heatmap PDFs
if ~isfolder(SavePath), mkdir(SavePath); end

%% Simulations Settings
% Fixed sample size as requested for signal structure experiments
nFixed = 300; 
% Grid of Signal Strength (sigStr) to experiment with alternative structures
sigStrGrid = 0.5; 
iter = 1;     % Increased for statistical stability
p = 60;        % Matrix dimension
sigStrBeta = 1;

%% Considered Methods
ConsideredMethods = {'LogisticSpINNEr', 'SpINNEr', 'LogisticLASSO', 'LogisticNUCLEAR'};
nMethods = length(ConsideredMethods);
SETTINGS = cell(length(ConsideredMethods), 1);

% ------ LogisticSpINNEr settings ------------
SETTINGS{1} = struct(...
    'Family', 'Binomial',...
    'Method', 'CV' ,...
    'UseParallel', true, ...
    'gridLengthN', 10, ...
    'gridLengthL', 10, ...
    'LambdaL', [], ...
    'LambdaN', [], ...
    'maxIters', 100, ...
    'displayStatus', false ...
    );

% ------ SpINNEr settings ------------
SETTINGS{2} = struct(...
    'Family', 'Gaussian',...
    'Method', 'CV' ,...
    'UseParallel', true, ...
    'gridLengthN', 10, ...
    'gridLengthL', 10, ...
    'LambdaL', [], ...
    'LambdaN', [], ...
    'maxIters', 100, ...
    'displayStatus', false ...
    );

% ------ LogisticLASSO settings ------------
SETTINGS{3} = struct(...
    'Family', 'Binomial',...
    'Method', 'CV' ,...
    'UseParallel', true, ...
    'gridLengthN', 10, ...
    'gridLengthL', 10, ...
    'LambdaL', [], ...
    'LambdaN', 0, ...
    'maxIters', 100, ...
    'displayStatus', false ...
    );

% ------ LogisticNUCLEAR settings ------------
SETTINGS{4} = struct(...
    'Family', 'Binomial',...
    'Method', 'CV' ,...
    'UseParallel', true, ...
    'gridLengthN', 10, ...
    'gridLengthL', 10, ...
    'LambdaL', 0, ...
    'LambdaN', [], ...
    'maxIters', 100, ...
    'displayStatus', false ...
    );

%% Main loops
NgridLength = length(sigStrGrid);
ACC = zeros(NgridLength, iter, nMethods);
LAM_N = zeros(NgridLength, iter, length(ConsideredMethods));
LAM_L = zeros(NgridLength, iter, length(ConsideredMethods));
OBJVAL_LOGLIK = zeros(NgridLength, iter, length(ConsideredMethods));
OBJVAL_NUCLEAR = zeros(NgridLength, iter, length(ConsideredMethods));
OBJVAL_LASSO = zeros(NgridLength, iter, length(ConsideredMethods));
TIME = zeros(NgridLength, iter, nMethods);
B_ESTIMS = cell(NgridLength, iter, nMethods);

% Initialize Progress Bar
h = waitbar(0, 'Starting simulations...', 'Name', 'Experiment Progress');
totalSteps = NgridLength * iter;

for ii = 1:NgridLength
    currentSigStr = sigStrGrid(ii);
    
    % --- Update True Signal Structure ---
    % Experimenting with block strength as per SpINNEr paper
    B1      = currentSigStr * ones(5,5);
    B2      = -currentSigStr * ones(5,5);
    B3      = currentSigStr * ones(4,4);
    sNods   = size(B1,1) + size(B2,1) + size(B3,1);
    B_true  = blkdiag(zeros(5,5), B1, zeros(6,6), B2, zeros(7,7), B3, zeros(p - sNods - 18, p - sNods - 18));
    Beta_true = sigStrBeta;

    for jj = 1:iter
        % Update waitbar
        step = ((ii-1) * iter) + jj;
        waitbar(step/totalSteps, h, sprintf('Signal Strength: %.2f | Iteration: %d/%d', currentSigStr, jj, iter));
        
        % Generate Data (Fixed n)
        AA = randn([p, p, nFixed]);
        AA = AA + permute(AA, [2 1 3]); 
        AA2 = zscore(reshape(AA, [p^2, nFixed])')';
        idxsDiag = logical(reshape(eye(p,p),[p^2,1]));
        AA2(idxsDiag,:) = 0;
        A = reshape(AA2, [p, p, nFixed]);
        X = ones(nFixed, 1);
        
        y = zeros(nFixed, 1);
        for i = 1:nFixed
            eta = trace(A(:,:,i)' * B_true) + X(i,:) * Beta_true;
            prob = 1 / (1 + exp(-eta));
            y(i) = rand < prob;
        end

        % Fit Methods
        for ee = 1:nMethods
            csettings = SETTINGS{ee};
            out = spinner(y, X, A);
            tic();
            out.fit(csettings);
            ctime = toc();
            % Collect the results
            out_est = out.Estimate;
            ACC(ii,jj,ee) = norm(B_true-out_est.B, 'fro');
            if ~isempty(out_est.bestLambdaN)
                LAM_N(ii,jj,ee) = out_est.bestLambdaN;
            else
                LAM_N(ii,jj,ee) = NaN;
            end
            if ~isempty(out_est.bestLambdaL)
                LAM_L(ii,jj,ee) = out_est.bestLambdaL;
            else
                LAM_L(ii,jj,ee) = NaN;
            end
            OBJVAL_LOGLIK(ii,jj,ee) = out_est.optVal_loglik;
            OBJVAL_NUCLEAR(ii,jj,ee) = out_est.optVal_lamNucl;
            OBJVAL_LASSO(ii,jj,ee) = out_est.optVal_lamLasso;
            TIME(ii,jj,ee) = ctime;
            B_ESTIMS{ii,jj,ee} = out_est.B;
        end
    end
end
close(h);

%% Save results
save(fullfile(SavePath, 'B_ESTIMS.mat'), 'B_ESTIMS');

%% Heatmaps
% original size
plotNames = {'logSpinner.pdf', 'Spinner.pdf', 'logLasso.pdf', 'logNuclear.pdf'};
for ii = 1:4
    spinnerHeatmapLocal(B_ESTIMS{1,1,ii})
    exportgraphics(gcf, fullfile(SavePath, plotNames{ii}), 'ContentType', 'vector');
    close all
end
spinnerHeatmapLocal(B_true)
exportgraphics(gcf, fullfile(SavePath, 'Btrue.pdf'), 'ContentType', 'vector');
close all

% reduced size
plotNames = {'logSpinnerRed.pdf', 'SpinnerRed.pdf', 'logLassoRed.pdf', 'logNuclearRed.pdf'};
for ii = 1:4
    cEstim = B_ESTIMS{1,1,ii};
    spinnerHeatmapLocal(cEstim(1:35, 1:35))
    exportgraphics(gcf, fullfile(SavePath, plotNames{ii}), 'ContentType', 'vector');
    close all
end
spinnerHeatmapLocal(B_true(1:35, 1:35))
exportgraphics(gcf, fullfile(SavePath, 'BtrueRed.pdf'), 'ContentType', 'vector');
close all






