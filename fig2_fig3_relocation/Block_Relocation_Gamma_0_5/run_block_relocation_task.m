function resultFile = run_block_relocation_task(taskId, numWorkers)
%RUN_BLOCK_RELOCATION_TASK Run one block-setting/replicate at gamma=0.5.
%   Block setting varies fastest: tasks 1:6 are replicate 1 and tasks
%   295:300 are replicate 50.

config = block_relocation_config();
[blockIndex, replicateId] = decode_task_id(taskId, config);
if nargin < 2 || isempty(numWorkers)
    numWorkers = config.defaultWorkers;
end
if ~isscalar(numWorkers) || numWorkers < 1 || numWorkers ~= floor(numWorkers)
    error('BlockRelocation:InvalidWorkerCount', ...
        'numWorkers must be a positive integer.');
end

configure_relocation_spinner_path(config.spinnerDir);
if ~isfolder(config.resultsDir)
    mkdir(config.resultsDir);
end

blockSizes = config.blockSizeGrid(blockIndex, :);
blockLabel = config.blockLabels{blockIndex};
resultFile = fullfile(config.resultsDir, sprintf( ...
    'relocation_gamma0p5_blocks_%s_replicate_%02d.mat', blockLabel, replicateId));

seeds = make_relocation_replicate_seeds(config.baseSeed, replicateId);
[BBaseline, signalInfo] = make_relocation_block_signal( ...
    config.p, blockSizes, config.coefficientMagnitude);
kGrid = round(config.swapProportionTargets .* ...
    signalInfo.nonzeroUpperEdgeCount);
if numel(unique(kGrid)) ~= numel(kGrid)
    error('BlockRelocation:DuplicateSwapCounts', ...
        'Rounded swap counts are not unique for blocks %s.', blockLabel);
end
[BTrue, swapInfo] = make_nested_relocation_signals( ...
    BBaseline, kGrid, seeds.swapSeed);
actualSwapProportions = swapInfo.swapProportions;

A = generate_relocation_standardized_design(config.p, config.n, seeds.dataSeed);
X = ones(config.n, 1);
[responses, responseUniforms, probabilities] = ...
    generate_relocation_paired_responses( ...
    A, X, BTrue, config.intercept, seeds.responseSeed);
validate_relocation_generated_inputs(config, blockSizes, BBaseline, signalInfo, ...
    BTrue, swapInfo, kGrid, A, responses, probabilities);

if isfile(resultFile)
    loaded = load(resultFile, 'result');
    if ~isfield(loaded, 'result')
        error('BlockRelocation:InvalidCheckpoint', ...
            'Checkpoint does not contain result: %s', resultFile);
    end
    result = loaded.result;
    validate_checkpoint(result, config, blockIndex, replicateId, ...
        BTrue, responses, swapInfo);
    fprintf('Resuming blocks %s, replicate %d from %s.\n', ...
        blockLabel, replicateId, resultFile);
else
    result = initialize_result(config, blockIndex, replicateId, seeds, ...
        signalInfo, kGrid, actualSwapProportions, BTrue, responses, ...
        responseUniforms, probabilities, swapInfo);
    result.numWorkersRequested = numWorkers;
    save_checkpoint(resultFile, result);
end

if all(result.completed(:))
    fprintf('Blocks %s, replicate %d are already complete.\n', ...
        blockLabel, replicateId);
    return;
end

ensure_parallel_pool(numWorkers, taskId);
upperMask = triu(true(config.p), 1);
for conditionIndex = 1:config.nConditions
    for methodIndex = 1:config.nMethods
        if result.completed(conditionIndex, methodIndex)
            continue;
        end

        methodName = config.methodNames{methodIndex};
        fprintf(['Blocks %s, replicate %d/%d, k=%d (%.6f), ', ...
            'method %d/%d (%s) started at %s.\n'], ...
            blockLabel, replicateId, config.nReplicates, ...
            kGrid(conditionIndex), actualSwapProportions(conditionIndex), ...
            methodIndex, config.nMethods, methodName, timestamp_utc());

        rng(seeds.cvSeed, 'twister');
        try
            estimate = fit_relocation_spinner_method( ...
                responses(:, conditionIndex), X, A, ...
                config.settings{methodIndex});
            currentBTrue = BTrue(:, :, conditionIndex);
            currentBHat = estimate.B;
            currentAbsBHat = abs(currentBHat);

            if ~isfield(estimate, 'beta') || ~isscalar(estimate.beta) || ...
                    ~isfinite(estimate.beta)
                error('BlockRelocation:InvalidIntercept', ...
                    'A finite scalar fitted intercept is required.');
            end
            result.betaHat(conditionIndex, methodIndex) = estimate.beta;
            if ismember(methodIndex, config.predictionMethodIndices)
                [auc, euclid, probs] = read_relocation_prediction_metrics( ...
                    estimate, responses(:, conditionIndex));
                result.AUC(conditionIndex, methodIndex) = auc;
                result.Euclid(conditionIndex, methodIndex) = euclid;
                result.Probs(:, conditionIndex, methodIndex) = probs;
            end

            result.BHat(:, :, conditionIndex, methodIndex) = currentBHat;
            result.ACC(conditionIndex, methodIndex) = ...
                norm(currentBTrue - currentBHat, 'fro') / ...
                norm(currentBTrue, 'fro');
            result.bestLambdaN(conditionIndex, methodIndex) = ...
                scalar_or_nan(estimate, 'bestLambdaN');
            result.bestLambdaL(conditionIndex, methodIndex) = ...
                scalar_or_nan(estimate, 'bestLambdaL');
            result.estimateRank(conditionIndex, methodIndex) = rank(currentBHat);
            result.estimateFrobeniusNorm(conditionIndex, methodIndex) = ...
                norm(currentBHat, 'fro');
            result.maxAbsCoefficient(conditionIndex, methodIndex) = ...
                max(currentAbsBHat(:));
            result.allZeroEstimate(conditionIndex, methodIndex) = ...
                result.maxAbsCoefficient(conditionIndex, methodIndex) <= ...
                config.zeroTolerance;
            result.zeroUpperEdgeCount(conditionIndex, methodIndex) = ...
                sum(currentAbsBHat(upperMask) <= config.zeroTolerance);
            result.completed(conditionIndex, methodIndex) = true;
            result.lastFailure = struct([]);
            result.updatedAtUtc = timestamp_utc();
            save_checkpoint(resultFile, result);
        catch methodError
            result.lastFailure = struct( ...
                'conditionIndex', conditionIndex, ...
                'k', kGrid(conditionIndex), ...
                'swapProportion', actualSwapProportions(conditionIndex), ...
                'methodIndex', methodIndex, ...
                'methodName', methodName, ...
                'identifier', methodError.identifier, ...
                'message', methodError.message, ...
                'report', getReport(methodError, 'extended', ...
                'hyperlinks', 'off'), ...
                'failedAtUtc', timestamp_utc());
            result.updatedAtUtc = timestamp_utc();
            save_checkpoint(resultFile, result);
            rethrow(methodError);
        end

        fprintf(['Blocks %s, replicate %d, k=%d, method %s ', ...
            'completed; ACC=%.10g.\n'], blockLabel, replicateId, ...
            kGrid(conditionIndex), methodName, ...
            result.ACC(conditionIndex, methodIndex));
    end
end

validate_relocation_saved_metrics(result, config);
result.finishedAtUtc = timestamp_utc();
result.updatedAtUtc = result.finishedAtUtc;
save_checkpoint(resultFile, result);
fprintf('Blocks %s, replicate %d completed at %s.\n', ...
    blockLabel, replicateId, result.finishedAtUtc);
end

function [blockIndex, replicateId] = decode_task_id(taskId, config)
if ~isscalar(taskId) || taskId < 1 || taskId > config.totalTasks || ...
        taskId ~= floor(taskId)
    error('BlockRelocation:InvalidTaskId', ...
        'taskId must be an integer from 1 through %d.', config.totalTasks);
end
blockIndex = mod(taskId - 1, config.nBlockSettings) + 1;
replicateId = floor((taskId - 1) / config.nBlockSettings) + 1;
end

function result = initialize_result(config, blockIndex, replicateId, seeds, ...
        signalInfo, kGrid, actualSwapProportions, BTrue, responses, ...
        responseUniforms, probabilities, swapInfo)
nConditions = config.nConditions;
nMethods = config.nMethods;

result.schemaVersion = config.schemaVersion;
result.configSignature = config.signature;
result.blockIndex = blockIndex;
result.blockSizes = config.blockSizeGrid(blockIndex, :);
result.blockLabel = config.blockLabels{blockIndex};
result.replicateId = replicateId;
result.n = config.n;
result.p = config.p;
result.intercept = config.intercept;
result.coefficientMagnitude = config.coefficientMagnitude;
result.signalInfo = signalInfo;
result.targetSwapProportions = config.swapProportionTargets;
result.kGrid = kGrid;
result.swapProportions = actualSwapProportions;
result.methodNames = config.methodNames;
result.settings = config.settings;
result.predictionMetricScope = config.predictionMetricScope;
result.seeds = seeds;
result.swapInfo = swapInfo;
result.BTrue = BTrue;
result.BHat = nan(config.p, config.p, nConditions, nMethods);
result.responses = responses;
result.responseUniforms = responseUniforms;
result.probabilities = probabilities;
result.responseMean = mean(responses, 1);
result.meanProbability = mean(probabilities, 1);
result.extremeProbabilityProportion = mean( ...
    probabilities < 0.01 | probabilities > 0.99, 1);
result.trueRank = zeros(1, nConditions);
result.trueFrobeniusNorm = zeros(1, nConditions);
result.trueNuclearNorm = zeros(1, nConditions);
result.trueL1Norm = zeros(1, nConditions);
for conditionIndex = 1:nConditions
    currentBTrue = BTrue(:, :, conditionIndex);
    result.trueRank(conditionIndex) = rank(currentBTrue);
    result.trueFrobeniusNorm(conditionIndex) = norm(currentBTrue, 'fro');
    result.trueNuclearNorm(conditionIndex) = sum(svd(currentBTrue));
    result.trueL1Norm(conditionIndex) = sum(abs(currentBTrue(:)));
end
result.ACC = nan(nConditions, nMethods);
result.AUC = nan(nConditions, nMethods);
result.Euclid = nan(nConditions, nMethods);
result.Probs = nan(config.n, nConditions, nMethods);
result.betaHat = nan(nConditions, nMethods);
result.bestLambdaN = nan(nConditions, nMethods);
result.bestLambdaL = nan(nConditions, nMethods);
result.estimateRank = nan(nConditions, nMethods);
result.estimateFrobeniusNorm = nan(nConditions, nMethods);
result.maxAbsCoefficient = nan(nConditions, nMethods);
result.allZeroEstimate = false(nConditions, nMethods);
result.zeroUpperEdgeCount = nan(nConditions, nMethods);
result.completed = false(nConditions, nMethods);
result.lastFailure = struct([]);
result.matlabVersion = version;
result.hostName = getenv('HOSTNAME');
result.createdAtUtc = timestamp_utc();
result.updatedAtUtc = result.createdAtUtc;
result.finishedAtUtc = '';
end

function validate_checkpoint(result, config, blockIndex, replicateId, ...
        BTrue, responses, swapInfo)
requiredFields = {'schemaVersion', 'configSignature', 'blockIndex', ...
    'blockSizes', 'blockLabel', 'replicateId', 'coefficientMagnitude', ...
    'kGrid', 'swapProportions', 'methodNames', 'swapInfo', 'BTrue', ...
    'BHat', 'responses', 'ACC', 'AUC', 'Euclid', 'Probs', 'betaHat', ...
    'predictionMetricScope', 'bestLambdaN', 'bestLambdaL', 'completed'};
for fieldIndex = 1:numel(requiredFields)
    if ~isfield(result, requiredFields{fieldIndex})
        error('BlockRelocation:InvalidCheckpoint', ...
            'Checkpoint is missing field %s.', requiredFields{fieldIndex});
    end
end

metadataMatches = ...
    result.schemaVersion == config.schemaVersion && ...
    strcmp(result.configSignature, config.signature) && ...
    result.blockIndex == blockIndex && ...
    isequal(result.blockSizes, config.blockSizeGrid(blockIndex, :)) && ...
    strcmp(result.blockLabel, config.blockLabels{blockIndex}) && ...
    result.replicateId == replicateId && ...
    result.coefficientMagnitude == config.coefficientMagnitude && ...
    strcmp(result.predictionMetricScope, config.predictionMetricScope) && ...
    isequal(result.methodNames, config.methodNames);
inputsMatch = ...
    isequaln(result.BTrue, BTrue) && ...
    isequaln(result.responses, responses) && ...
    isequaln(result.swapInfo.sourceEdgeOrder, swapInfo.sourceEdgeOrder) && ...
    isequaln(result.swapInfo.targetEdgeOrder, swapInfo.targetEdgeOrder);
if ~metadataMatches || ~inputsMatch
    error('BlockRelocation:CheckpointMismatch', ...
        'Regenerated inputs or metadata do not match the checkpoint.');
end
validate_relocation_saved_metrics(result, config);
end

function value = scalar_or_nan(estimate, fieldName)
value = NaN;
if isfield(estimate, fieldName) && ~isempty(estimate.(fieldName))
    candidate = estimate.(fieldName);
    if ~isnumeric(candidate) || ~isscalar(candidate) || ~isfinite(candidate)
        error('BlockRelocation:InvalidSelectedLambda', ...
            'Estimate field %s must be empty or a finite scalar.', fieldName);
    end
    value = candidate;
end
end

function ensure_parallel_pool(numWorkers, taskId)
if ~license('test', 'Distrib_Computing_Toolbox')
    error('BlockRelocation:MissingParallelToolbox', ...
        'Parallel Computing Toolbox is required because UseParallel is true.');
end
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= numWorkers
    error('BlockRelocation:WorkerCountMismatch', ...
        'Active pool has %d workers, but %d were requested.', ...
        pool.NumWorkers, numWorkers);
end
if ~isempty(pool)
    return;
end

maxAttempts = 3;
for attempt = 1:maxAttempts
    cluster = parcluster('local');
    jobStorageDir = make_job_storage_directory(taskId, attempt);
    if ~isfolder(jobStorageDir)
        mkdir(jobStorageDir);
    end
    cluster.JobStorageLocation = jobStorageDir;
    fprintf('Starting pool attempt %d/%d with storage %s.\n', ...
        attempt, maxAttempts, jobStorageDir);
    try
        parpool(cluster, numWorkers);
        return;
    catch poolError
        pool = gcp('nocreate');
        if ~isempty(pool)
            delete(pool);
        end
        if attempt == maxAttempts
            rethrow(poolError);
        end
        warning('BlockRelocation:ParallelPoolRetry', ...
            'Pool attempt %d/%d failed: %s', ...
            attempt, maxAttempts, poolError.message);
        pause(15 * attempt);
    end
end
end

function jobStorageDir = make_job_storage_directory(taskId, attempt)
arrayJobId = getenv('SLURM_ARRAY_JOB_ID');
arrayTaskId = getenv('SLURM_ARRAY_TASK_ID');
slurmJobId = getenv('SLURM_JOB_ID');
if isempty(arrayJobId)
    arrayJobId = 'interactive';
end
if isempty(arrayTaskId)
    arrayTaskId = sprintf('task-%03d', taskId);
end
if isempty(slurmJobId)
    slurmJobId = sprintf('matlab-%d', feature('getpid'));
end
jobStorageDir = fullfile(tempdir, sprintf( ...
    'relocation-gamma0p5-storage-%s-%s-%s-attempt-%d', ...
    arrayJobId, arrayTaskId, slurmJobId, attempt));
end

function save_checkpoint(resultFile, result)
resultsDir = fileparts(resultFile);
temporaryFile = [tempname(resultsDir), '.mat'];
cleanupTemporary = onCleanup(@() delete_if_present(temporaryFile)); %#ok<NASGU>
save(temporaryFile, 'result', '-v7.3');
[moveSucceeded, moveMessage] = movefile(temporaryFile, resultFile, 'f');
if ~moveSucceeded
    error('BlockRelocation:CheckpointWriteFailed', ...
        'Could not replace %s: %s', resultFile, moveMessage);
end
end

function delete_if_present(filePath)
if isfile(filePath)
    delete(filePath);
end
end

function value = timestamp_utc()
value = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
end
