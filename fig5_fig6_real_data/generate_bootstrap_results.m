% Bootstrap fits of Logistic SpINNEr on the real data (Figs. 5-6).
%
% Data (not included in this repository): data/GSRandNonGSRarrays.mat with
% 214 x 214 x 161 resting-state connectivity arrays (200 Schaefer cortical
% regions + 14 subcortical), of which this script uses AA_NOGSR_low (no global
% signal regression, low-sucrose session), and data/gen_hist_alcohol.csv with
% family history of alcoholism (FHP / FHN) in column x. One subject's matrix
% has missing values and is dropped, giving n = 160 (78 FH+, 82 FH-).
%
% The 300 bootstrap fits were run in five parts (35, 35, 30, 90 and 110
% bootstrap samples). This is the part-5 version: it draws and discards 190
% samples, then fits numBT = 110. No random seed was set.
%--------------------------------------------------------------------------
% Clear the workspace
clear

% Paths
here = fileparts(mfilename('fullpath'));
DataPath = fullfile(here, 'data');                                % put the data files here
addpath(fullfile(here, '..', 'software_84509f9'))                 % SpINNEr toolbox at commit 84509f9
savePath = fullfile(here, 'results');
if ~isfolder(savePath), mkdir(savePath); end

% Load the data
data = load(fullfile(DataPath, 'GSRandNonGSRarrays.mat'));
gen_alcohol = readtable(fullfile(DataPath, 'gen_hist_alcohol.csv')); %    'FHN' and 'FHP' entries only

% check distribution
FHNidxs = ismember(gen_alcohol.x,'FHN'); % 83
FHPidxs = ismember(gen_alcohol.x,'FHP'); % 78

% number of Bootstrap samples
numBT = 110;

%-----------------------------------------------
% SpINNEr settings
Settings = struct(...
    'displayStatus', false, ...
    'Family', 'Binomial',...
    'Method', 'CV' ,...
    'UseParallel', true, ...
    'gridLengthN', 15, ...
    'gridLengthL', 15, ...
    'searchDensity', 'New', ...
    'expParameterL', 3, ...
    'expParameterN', 3, ...
    'gridMaxN', 50, ...
    'maxIters', 100, ...
    'addZerosInGrids', false, ...
    'VectorizedParallelization', false ...
    );

%% Run logistic SpINNEr
%--------------------------------------------------------------------------
%             NoN - Global Signal Regression  &   Low Sucrose
%--------------------------------------------------------------------------

% Prepare the data
n0 = size(data.AA_NOGSR_low,3);
A0 = data.AA_NOGSR_low;
y0 = FHPidxs; % family history of alcoholism: 1 --> Positive, 0 --> Negative 

% Identify predictor matrices with missing data
containsMissing = zeros(n0,1);
for ii = 1:n0
    containsMissing(ii) = nnz(isnan(A0(:,:,ii)));
end

% final sample size
n = sum(~containsMissing);

% eliminate predictors and observations with missing data
A = A0(:,:,~containsMissing);
y = y0(~containsMissing);
X = ones(n,1);

% make matrices symmetric and with zeros on diagonals
for ii = 1:n
    A(:,:,ii) = (A(:,:,ii) + A(:,:,ii)')/2;
    A(:,:,ii) = A(:,:,ii) - diag(diag(A(:,:,ii)));
end

p = size(A,1);

%% Generate bootstrap results
estimatesBT = zeros(p, p, numBT);
allResults = cell(numBT,1);
FLAG_INTERUPTED = 0;

for bb = 1:190
    sample_BT = randsample(n, n, true);
end

try
    for bb = 1:numBT
        % random permutation
        sample_BT = randsample(n, n, true);
        A_BT = A(:, :, sample_BT);
        y_BT = y(sample_BT);
        
        % standardize the data
        for ii = 1:p
            for jj = 1:p
                cvect = squeeze(A_BT(ii,jj,:));
                if any(cvect)
                    cvect = cvect-mean(cvect);
                    cvect = cvect./sqrt(sum(cvect.^2));
                    A_BT(ii,jj,:) = cvect;
                end
            end
        end
            
        sp = spinner(y_BT, X, A_BT);
        sp.fit(Settings(:));
    
        estimatesBT(:,:,bb) = sp.Estimate.B;
        allResults{bb} = sp.Estimate;
        sprintf('*** FINALIZED BOOTSTRAP ITERATIONS: %d  ***', bb)
    end
catch
    FLAG_INTERUPTED = 1;
    save(fullfile(savePath, 'bootstrap_results_part5_stopped.mat'), 'estimatesBT')
    save(fullfile(savePath, 'allResults_part5_stopped.mat'), 'allResults')
end

if FLAG_INTERUPTED == 0
    save(fullfile(savePath, 'bootstrap_results_part5.mat'), 'estimatesBT')
    save(fullfile(savePath, 'allResults_part5.mat'), 'allResults')
end


% 
% 
% Settings = struct(...
%     'Family', 'Binomial',...
%     'Method', 'CV' ,...
%     'UseParallel', true, ...
%     'gridLengthN', 15, ...
%     'gridLengthL', 15, ...
%     'searchDensity', 'New', ...
%     'expParameterL', 3, ...
%     'expParameterN', 3, ...
%     'gridMaxN', 50, ...
%     'maxIters', 100, ...
%     'addZerosInGrids', false, ...
%     'VectorizedParallelization', false, ...
%     'LambdaL', 10, ...
%     'LambdaN', 6 ...
%     );
% 
% sp = spinner(y_BT, X, A_BT);
% sp.fit(Settings(:))




% spinnerHeatmap(sp.Estimate.B)
% spinnerHeatmap(sp.Estimate.logliksCV)
% spinnerHeatmap(sp.Estimate.logliksCV == min(min(sp.Estimate.logliksCV)))

% save the results
% save(fullfile(savePath,'LogSpinner_NGSR_LOW.mat'), "sp")