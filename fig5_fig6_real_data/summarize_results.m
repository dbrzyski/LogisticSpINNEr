clear

%% Paths
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'software_84509f9'))   % for spinnerHeatmap
resPath = fullfile(here, 'results');                % output of the five parts of generate_bootstrap_results.m

%% Load the data
estim1 = load(fullfile(resPath,'bootstrap_results_part1.mat'));
estim2 = load(fullfile(resPath,'bootstrap_results_part2.mat'));
estim3 = load(fullfile(resPath,'bootstrap_results_part3.mat'));
estim4 = load(fullfile(resPath,'bootstrap_results_part4.mat'));
estim5 = load(fullfile(resPath,'bootstrap_results_part5.mat'));

results1 = load(fullfile(resPath,'allResults_part1.mat'));
results2 = load(fullfile(resPath,'allResults_part2.mat'));
results3 = load(fullfile(resPath,'allResults_part3.mat'));
results4 = load(fullfile(resPath,'allResults_part4.mat'));
results5 = load(fullfile(resPath,'allResults_part5.mat'));

%% Merge the results
estimAll = cat(3, estim1.estimatesBT, estim2.estimatesBT, estim3.estimatesBT, estim4.estimatesBT, estim5.estimatesBT);
resultsAll = cat(1, results1.allResults, results2.allResults, results3.allResults, results4.allResults, results5.allResults);

save(fullfile(resPath, 'EstimatesBootstrap.mat'), 'estimAll')
save(fullfile(resPath, 'ResultsBootstrap.mat'), 'resultsAll')

estimAll1 = estimAll(:,:,1:150);
estimAll2 = estimAll(:,:,151:300);
save(fullfile(resPath, 'EstimatesBootstrap_P1.mat'), 'estimAll1')
save(fullfile(resPath, 'EstimatesBootstrap_P2.mat'), 'estimAll2')


%% Analysis
n = size(estimAll,3);
SIGNIFF_LEVEL = floor(0.95*n);
aa = sum(estimAll~=0,3);
aa = aa-diag(diag(aa));
[bb,cc]=find(aa>=SIGNIFF_LEVEL);

iii = unique([bb;cc]);


estimAllMed = median(estimAll,3);





%% Clustering
ClusterMatrix = estimAllMed(iii, iii);

PosNodes = sum(ClusterMatrix,1) > 0;
NegNodes = sum(ClusterMatrix,1) <= 0;

% ----------- Positive Nodes
% LINKAGE_TYPE = 'complete';
LINKAGE_TYPE = 'average';
% LINKAGE_TYPE = 'single';
THRESH_PROP = 0;

% Thresholding
cB0 = ClusterMatrix(PosNodes, PosNodes);
cB = abs(cB0);
maxEntry = max(max(cB));
FiltThresh = THRESH_PROP*maxEntry;
cB(cB < FiltThresh) = 0;

% Similarity matrix
similarity = cB;
similarity = similarity - min(min(similarity));
similarity(similarity>1) = 1;
similarity = similarity/max(max(similarity));

% Distance matrix
distance = 1 - similarity;
distance = distance - diag(diag(distance));

% Hierarchical clustering: leaf order
Y = squareform(distance);
Z = linkage(Y, LINKAGE_TYPE);
leafOrder = optimalleaforder(Z, Y);

spinnerHeatmap(cB0(leafOrder,leafOrder))


% ----------- Negative Nodes
% LINKAGE_TYPE = 'complete';
LINKAGE_TYPE = 'average';
% LINKAGE_TYPE = 'single';
THRESH_PROP = 0;

% Thresholding
cB0 = ClusterMatrix(NegNodes, NegNodes);
cB = abs(cB0);
maxEntry = max(max(cB));
FiltThresh = THRESH_PROP*maxEntry;
cB(cB < FiltThresh) = 0;

% Similarity matrix
similarity = cB;
similarity = similarity - min(min(similarity));
similarity(similarity>1) = 1;
similarity = similarity/max(max(similarity));

% Distance matrix
distance = 1 - similarity;
distance = distance - diag(diag(distance));

% Hierarchical clustering: leaf order
Y = squareform(distance);
Z = linkage(Y, LINKAGE_TYPE);
leafOrder = optimalleaforder(Z, Y);

spinnerHeatmap(cB0(leafOrder,leafOrder))

