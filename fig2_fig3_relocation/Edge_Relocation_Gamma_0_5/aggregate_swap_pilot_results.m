function summaryFile = aggregate_swap_pilot_results(outputDirectory)
%AGGREGATE_SWAP_PILOT_RESULTS Combine all completed ACC result arrays.

config = swap_pilot_config();
if nargin < 1 || isempty(outputDirectory)
    outputDirectory = fullfile(config.simulationDir, 'recomputed_summaries');
end

nConditions = numel(config.kGrid);
nMethods = numel(config.methodNames);
ACC = nan(nConditions, config.nReplicates, nMethods);

for replicateId = 1:config.nReplicates
    replicateFile = fullfile(config.resultsDir, ...
        sprintf('swap_pilot_replicate_%02d.mat', replicateId));
    if ~isfile(replicateFile)
        error('SwapPilot:MissingReplicate', ...
            'Missing replicate result: %s', replicateFile);
    end

    loaded = load(replicateFile, 'result');
    result = loaded.result;
    if ~strcmp(result.configSignature, config.signature) || ...
            ~all(result.completed(:))
        error('SwapPilot:IncompleteReplicate', ...
            'Replicate %d is incomplete or has a different configuration.', ...
            replicateId);
    end
    ACC(:, replicateId, :) = reshape( ...
        result.ACC, [nConditions, 1, nMethods]);
end

meanACC = squeeze(mean(ACC, 2));
semACC = squeeze(std(ACC, 0, 2)) ./ sqrt(config.nReplicates);
twiceSemACC = 2 .* semACC;
kGrid = config.kGrid;
swapProportions = config.swapProportions;
methodNames = config.methodNames;
nReplicates = config.nReplicates;

% Delivery output routing; supplied input checkpoints stay unchanged.
config.resultsDir = outputDirectory;
if ~isfolder(config.resultsDir)
    mkdir(config.resultsDir);
end
summaryFile = fullfile(config.resultsDir, 'swap_pilot_ACC_summary.mat');
save(summaryFile, 'ACC', 'meanACC', 'semACC', 'twiceSemACC', ...
    'kGrid', 'swapProportions', 'methodNames', 'nReplicates');

write_long_acc_table(config, ACC);
write_summary_acc_table(config, meanACC, semACC, twiceSemACC);

fprintf('ACC summary written to %s.\n', summaryFile);
end

function write_long_acc_table(config, ACC)
nRows = numel(config.kGrid) * config.nReplicates * numel(config.methodNames);
replicateColumn = zeros(nRows, 1);
kColumn = zeros(nRows, 1);
proportionColumn = zeros(nRows, 1);
methodColumn = cell(nRows, 1);
accColumn = zeros(nRows, 1);

rowIndex = 0;
for conditionIndex = 1:numel(config.kGrid)
    for replicateId = 1:config.nReplicates
        for methodIndex = 1:numel(config.methodNames)
            rowIndex = rowIndex + 1;
            replicateColumn(rowIndex) = replicateId;
            kColumn(rowIndex) = config.kGrid(conditionIndex);
            proportionColumn(rowIndex) = config.swapProportions(conditionIndex);
            methodColumn{rowIndex} = config.methodNames{methodIndex};
            accColumn(rowIndex) = ACC(conditionIndex, replicateId, methodIndex);
        end
    end
end

accTable = table(replicateColumn, kColumn, proportionColumn, ...
    methodColumn, accColumn, 'VariableNames', ...
    {'Replicate', 'K', 'SwapProportion', 'Method', 'ACC'});
writetable(accTable, fullfile(config.resultsDir, 'swap_pilot_ACC_long.csv'));
end

function write_summary_acc_table(config, meanACC, semACC, twiceSemACC)
nRows = numel(config.kGrid) * numel(config.methodNames);
kColumn = zeros(nRows, 1);
proportionColumn = zeros(nRows, 1);
methodColumn = cell(nRows, 1);
meanColumn = zeros(nRows, 1);
semColumn = zeros(nRows, 1);
twiceSemColumn = zeros(nRows, 1);

rowIndex = 0;
for conditionIndex = 1:numel(config.kGrid)
    for methodIndex = 1:numel(config.methodNames)
        rowIndex = rowIndex + 1;
        kColumn(rowIndex) = config.kGrid(conditionIndex);
        proportionColumn(rowIndex) = config.swapProportions(conditionIndex);
        methodColumn{rowIndex} = config.methodNames{methodIndex};
        meanColumn(rowIndex) = meanACC(conditionIndex, methodIndex);
        semColumn(rowIndex) = semACC(conditionIndex, methodIndex);
        twiceSemColumn(rowIndex) = twiceSemACC(conditionIndex, methodIndex);
    end
end

summaryTable = table(kColumn, proportionColumn, methodColumn, ...
    meanColumn, semColumn, twiceSemColumn, 'VariableNames', ...
    {'K', 'SwapProportion', 'Method', 'MeanACC', 'SEM', 'TwoSEM'});
writetable(summaryTable, fullfile(config.resultsDir, 'swap_pilot_ACC_summary.csv'));
end

