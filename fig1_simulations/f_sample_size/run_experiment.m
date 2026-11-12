%% Clear the workspace
clear

%% Paths
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'software_84509f9'))   % SpINNEr toolbox at commit 84509f9, used for Fig. 1
SavePath = here;   % results.mat is written next to this script

%% Simulations Settings
nGrid = [100, 120, 140, 160, 180, 200, 220, 240];
% nGrid = 150;      % considered samples numbers
iter = 50;     % number of repetitions for each n
p = 60;        % matrix dimension
sigStr = 1;
sigStrBeta = 1;

%% True signal 
rng(2)
B1      = sigStr*ones(5,5);
B2      = -sigStr*ones(5,5);
B3      = sigStr*ones(4,4);
sNods   = size(B1,1) + size(B2,1) + size(B3,1);
B_true       = blkdiag(zeros(5,5), B1, zeros(6,6), B2, zeros(7,7), B3, zeros(p - sNods - 18, p- sNods - 18));
Beta_true = sigStrBeta;

%% Considered Methods
ConsideredMethods = {'LogisticSpINNEr', 'SpINNEr', 'LogisticLASSO', 'LogisticNUCLEAR'};
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
% Placeholders for the results
NgridLength = length(nGrid);
ACC = zeros(NgridLength, iter, length(ConsideredMethods));
LAM_N = zeros(NgridLength, iter, length(ConsideredMethods));
LAM_L = zeros(NgridLength, iter, length(ConsideredMethods));
OBJVAL_LOGLIK = zeros(NgridLength, iter, length(ConsideredMethods));
OBJVAL_NUCLEAR = zeros(NgridLength, iter, length(ConsideredMethods));
OBJVAL_LASSO = zeros(NgridLength, iter, length(ConsideredMethods));
TIME = zeros(NgridLength, iter, length(ConsideredMethods));

for ii = 1:NgridLength
    n = nGrid(ii);
    for jj = 1:iter
        % regressor matrices 
        AA               = randn([p, p, n]);
        AA               = AA + permute(AA, [2 1 3]); % matrices are symmetric
        AA2              = reshape(AA, [p^2, n]);
        AA2              = zscore(AA2')';             % standardize data
        idxsDiag         = logical(reshape(eye(p,p),[p^2,1]));
        AA2(idxsDiag,:)  = 0;                         % zeros on the diagonals
        A               = reshape(AA2, [p, p, n]);   % 3-way tensor of symmetric matrices with zeros on diagonals and standardized across 3rd dimension
        
        % Generate X as intercept column
        X = ones(n, 1);
        
        % Compute linear predictors and generate binary responses
        y = zeros(n, 1);
        for i = 1:n
            Ai = A(:, :, i);
            xi = X(i, :)';
            eta = trace(Ai' * B_true) + xi' * Beta_true;
            prob = 1 / (1 + exp(-eta));
            y(i) = rand < prob;  % Bernoulli draw
        end

        % Get estimates
        for ee = 1:length(ConsideredMethods)
            csettings = SETTINGS{ee};
            out = spinner(y, X, A);
            tic()
            out.fit(csettings(:));
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
        end
    end
    disp(n)
end

%% Save results
save(fullfile(SavePath, 'results.mat'), 'ACC', 'LAM_N', 'LAM_L', 'OBJVAL_LOGLIK', 'OBJVAL_NUCLEAR', 'OBJVAL_LASSO', 'TIME');

%% Plots
% Compute mean ACC over iterations (for generality, though iter=1)
mean_ACC = mean(ACC, 2); % average over jj dimension
mean_ACC = squeeze(mean_ACC); % now NgridLength x length(ConsideredMethods)

figure;
hold on;
colors = lines(4); % Use MATLAB's default color scheme for nice, distinct colors
linestyles = {'-', '--', ':', '-.'}; % Different line styles
markers = {'o', 's', 'd', '^'}; % Different markers
for ee = 1:length(ConsideredMethods)
    plot(nGrid, mean_ACC(:,ee), 'Color', colors(ee,:), 'LineWidth', 2, 'LineStyle', linestyles{ee}, 'Marker', markers{ee}, 'MarkerFaceColor', colors(ee,:));
end
legend(ConsideredMethods, 'Location', 'southwest');
xlabel('Sample size n', 'FontSize', 14);
ylabel('Frobenius norm error', 'FontSize', 14);
% title('Performance comparison of methods', 'FontSize', 16);
grid on;
hold off;

% Save the plot as high-quality PNG
exportgraphics(gcf, fullfile(SavePath, 'performance_plot_notitle.pdf'), 'Resolution', 300);