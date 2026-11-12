function resultFile = run_swap_pilot_replicate(replicateId, numWorkers)
%RUN_SWAP_PILOT_REPLICATE Run and checkpoint one paired pilot replicate.

config = swap_pilot_config();
validate_replicate_id(replicateId, config.nReplicates);

if nargin < 2 || isempty(numWorkers)
    numWorkers = config.defaultWorkers;
end
if ~isscalar(numWorkers) || numWorkers < 1 || numWorkers ~= floor(numWorkers)
    error('SwapPilot:InvalidWorkerCount', ...
        'numWorkers must be a positive integer.');
end

configure_spinner_path(config.spinnerDir);

if ~isfolder(config.resultsDir)
    mkdir(config.resultsDir);
end
resultFile = fullfile(config.resultsDir, ...
    sprintf('swap_pilot_replicate_%02d.mat', replicateId));

seeds = make_replicate_seeds(config.baseSeed, replicateId);
BBaseline = make_baseline_signal(config.p, config.signalStrength);
[BTrue, swapInfo] = make_nested_swapped_signals( ...
    BBaseline, config.kGrid, seeds.swapSeed);
A = generate_standardized_design(config.p, config.n, seeds.dataSeed);
X = ones(config.n, 1);
[responses, responseUniforms, probabilities] = generate_paired_responses( ...
    A, X, BTrue, config.intercept, seeds.responseSeed);

validate_generated_inputs(config, BBaseline, BTrue, swapInfo, A, responses);

if isfile(resultFile)
    loaded = load(resultFile, 'result');
    if ~isfield(loaded, 'result')
        error('SwapPilot:InvalidCheckpoint', ...
            'The checkpoint does not contain a result structure: %s', resultFile);
    end
    result = loaded.result;
    validate_checkpoint(result, config, replicateId, BTrue, responses, swapInfo);
    fprintf('Resuming replicate %d from %s.\n', replicateId, resultFile);
else
    result = initialize_result(config, replicateId, seeds, BTrue, responses, ...
        responseUniforms, probabilities, swapInfo);
    result.numWorkersRequested = numWorkers;
    save_checkpoint(resultFile, result);
end

if all(result.completed(:))
    fprintf('Replicate %d is already complete: %s.\n', replicateId, resultFile);
    return;
end

ensure_parallel_pool(numWorkers, replicateId);

nConditions = numel(config.kGrid);
nMethods = numel(config.methodNames);

for conditionIndex = 1:nConditions
    for methodIndex = 1:nMethods
        if result.completed(conditionIndex, methodIndex)
            continue;
        end

        k = config.kGrid(conditionIndex);
        methodName = config.methodNames{methodIndex};
        fprintf(['Replicate %d/%d, k=%d, method %d/%d (%s) started ', ...
            'at %s.\n'], replicateId, config.nReplicates, k, methodIndex, ...
            nMethods, methodName, timestamp_utc());

        rng(seeds.cvSeed, 'twister');
        try
            estimate = fit_spinner_method( ...
                responses(:, conditionIndex), X, A, config.settings{methodIndex});
            currentBTrue = BTrue(:, :, conditionIndex);
            result.BHat(:, :, conditionIndex, methodIndex) = estimate.B;
            result.ACC(conditionIndex, methodIndex) = ...
                norm(currentBTrue - estimate.B, 'fro') / norm(currentBTrue, 'fro');
            result.completed(conditionIndex, methodIndex) = true;
            result.lastFailure = struct([]);
            result.updatedAtUtc = timestamp_utc();
            save_checkpoint(resultFile, result);
        catch methodError
            result.lastFailure = struct( ...
                'conditionIndex', conditionIndex, ...
                'k', k, ...
                'methodIndex', methodIndex, ...
                'methodName', methodName, ...
                'identifier', methodError.identifier, ...
                'message', methodError.message, ...
                'report', getReport(methodError, 'extended', 'hyperlinks', 'off'), ...
                'failedAtUtc', timestamp_utc());
            result.updatedAtUtc = timestamp_utc();
            save_checkpoint(resultFile, result);
            rethrow(methodError);
        end

        fprintf('Replicate %d, k=%d, method %s completed; ACC=%.10g.\n', ...
            replicateId, k, methodName, result.ACC(conditionIndex, methodIndex));
    end
end

result.finishedAtUtc = timestamp_utc();
result.updatedAtUtc = result.finishedAtUtc;
save_checkpoint(resultFile, result);
fprintf('Replicate %d completed at %s.\n', replicateId, result.finishedAtUtc);
end

function validate_replicate_id(replicateId, nReplicates)
if ~isscalar(replicateId) || replicateId < 1 || ...
        replicateId > nReplicates || replicateId ~= floor(replicateId)
    error('SwapPilot:InvalidReplicate', ...
        'replicateId must be an integer from 1 through %d.', nReplicates);
end
end

function ensure_parallel_pool(numWorkers, replicateId)
if ~license('test', 'Distrib_Computing_Toolbox')
    error('SwapPilot:MissingParallelToolbox', ...
        'Parallel Computing Toolbox is required because UseParallel is true.');
end

pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= numWorkers
    error('SwapPilot:WorkerCountMismatch', ...
        'The active pool has %d workers, but %d were requested.', ...
        pool.NumWorkers, numWorkers);
end
if ~isempty(pool)
    return;
end

maxAttempts = 3;
for attempt = 1:maxAttempts
    cluster = parcluster('local');
    jobStorageDir = make_job_storage_directory(replicateId, attempt);
    if ~isfolder(jobStorageDir)
        mkdir(jobStorageDir);
    end
    cluster.JobStorageLocation = jobStorageDir;

    fprintf(['Starting MATLAB parallel pool, attempt %d/%d, ', ...
        'job storage: %s.\n'], attempt, maxAttempts, jobStorageDir);
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
        warning('SwapPilot:ParallelPoolRetry', ...
            ['Parallel pool attempt %d/%d failed: %s. ', ...
            'Retrying with new job storage.'], ...
            attempt, maxAttempts, poolError.message);
        pause(15 * attempt);
    end
end
end

function jobStorageDir = make_job_storage_directory(replicateId, attempt)
arrayJobId = getenv('SLURM_ARRAY_JOB_ID');
arrayTaskId = getenv('SLURM_ARRAY_TASK_ID');
slurmJobId = getenv('SLURM_JOB_ID');

if isempty(arrayJobId)
    arrayJobId = 'interactive';
end
if isempty(arrayTaskId)
    arrayTaskId = sprintf('replicate-%02d', replicateId);
end
if isempty(slurmJobId)
    slurmJobId = sprintf('matlab-%d', feature('getpid'));
end

directoryName = sprintf( ...
    'swap-matlab-storage-%s-%s-%s-attempt-%d', ...
    arrayJobId, arrayTaskId, slurmJobId, attempt);
jobStorageDir = fullfile(tempdir, directoryName);
end

function result = initialize_result(config, replicateId, seeds, BTrue, ...
        responses, responseUniforms, probabilities, swapInfo)
nConditions = numel(config.kGrid);
nMethods = numel(config.methodNames);

result.schemaVersion = config.schemaVersion;
result.configSignature = config.signature;
result.replicateId = replicateId;
result.n = config.n;
result.p = config.p;
result.signalStrength = config.signalStrength;
result.intercept = config.intercept;
result.kGrid = config.kGrid;
result.swapProportions = config.swapProportions;
result.methodNames = config.methodNames;
result.settings = config.settings;
result.seeds = seeds;
result.swapInfo = swapInfo;
result.BTrue = BTrue;
result.BHat = nan(config.p, config.p, nConditions, nMethods);
result.responses = responses;
result.responseUniforms = responseUniforms;
result.probabilities = probabilities;
result.ACC = nan(nConditions, nMethods);
result.completed = false(nConditions, nMethods);
result.lastFailure = struct([]);
result.matlabVersion = version;
result.hostName = getenv('HOSTNAME');
result.createdAtUtc = timestamp_utc();
result.updatedAtUtc = result.createdAtUtc;
result.finishedAtUtc = '';
end

function validate_checkpoint(result, config, replicateId, BTrue, responses, swapInfo)
requiredFields = {'schemaVersion', 'configSignature', 'replicateId', ...
    'kGrid', 'swapProportions', 'methodNames', 'swapInfo', 'BTrue', ...
    'BHat', 'responses', 'ACC', 'completed'};
for fieldIndex = 1:numel(requiredFields)
    if ~isfield(result, requiredFields{fieldIndex})
        error('SwapPilot:InvalidCheckpoint', ...
            'The checkpoint is missing field %s.', requiredFields{fieldIndex});
    end
end

if ~strcmp(result.configSignature, config.signature) || ...
        result.schemaVersion ~= config.schemaVersion || ...
        result.replicateId ~= replicateId
    error('SwapPilot:CheckpointMismatch', ...
        'The checkpoint does not match this simulation configuration.');
end
if ~isequaln(result.BTrue, BTrue) || ...
        ~isequaln(result.responses, responses) || ...
        ~isequaln(result.swapInfo.sourceEdgeOrder, swapInfo.sourceEdgeOrder) || ...
        ~isequaln(result.swapInfo.targetEdgeOrder, swapInfo.targetEdgeOrder)
    error('SwapPilot:CheckpointMismatch', ...
        'Regenerated inputs do not match the saved checkpoint.');
end
end

function validate_generated_inputs(config, BBaseline, BTrue, swapInfo, A, responses)
upperMask = triu(true(config.p), 1);
baselineSupport = upperMask & (BBaseline ~= 0);
baselineZero = upperMask & (BBaseline == 0);
baselineNorm = norm(BBaseline, 'fro');

if nnz(baselineSupport) ~= config.nonzeroUpperEdgeCount
    error('SwapPilot:UnexpectedBaseline', ...
        'The baseline must have exactly %d nonzero upper-triangular edges.', ...
        config.nonzeroUpperEdgeCount);
end
if ~isequal(swapInfo.swapProportions, config.swapProportions)
    error('SwapPilot:UnexpectedSwapProportions', ...
        'The generated swap proportions do not match the configuration.');
end

for conditionIndex = 1:numel(config.kGrid)
    currentSignal = BTrue(:, :, conditionIndex);
    k = config.kGrid(conditionIndex);
    currentSupport = upperMask & (currentSignal ~= 0);

    if norm(currentSignal - currentSignal.', 'fro') > 1e-12 || ...
            nnz(currentSupport) ~= config.nonzeroUpperEdgeCount || ...
            abs(norm(currentSignal, 'fro') - baselineNorm) > 1e-12 || ...
            nnz(baselineSupport & ~currentSupport) ~= k || ...
            nnz(baselineZero & currentSupport) ~= k || ...
            ~isequal(diag(currentSignal), diag(BBaseline))
        error('SwapPilot:InvalidSwappedSignal', ...
            'The signal for k=%d failed a structural invariant.', k);
    end
end

if max(abs(A(:) - reshape(permute(A, [2, 1, 3]), [], 1))) > 1e-12
    error('SwapPilot:InvalidDesign', 'A must be symmetric.');
end
for observationIndex = 1:size(A, 3)
    if any(diag(A(:, :, observationIndex)) ~= 0)
        error('SwapPilot:InvalidDesign', 'Every A matrix must have a zero diagonal.');
    end
end
if ~all(ismember(responses(:), [0, 1]))
    error('SwapPilot:InvalidResponse', 'All responses must be binary.');
end
if any(~isfinite(A(:)))
    error('SwapPilot:InvalidDesign', 'A must contain only finite values.');
end
end

function save_checkpoint(resultFile, result)
resultsDir = fileparts(resultFile);
temporaryFile = [tempname(resultsDir), '.mat'];
cleanupTemporary = onCleanup(@() delete_if_present(temporaryFile)); %#ok<NASGU>
save(temporaryFile, 'result', '-v7.3');
[moveSucceeded, moveMessage] = movefile(temporaryFile, resultFile, 'f');
if ~moveSucceeded
    error('SwapPilot:CheckpointWriteFailed', ...
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
