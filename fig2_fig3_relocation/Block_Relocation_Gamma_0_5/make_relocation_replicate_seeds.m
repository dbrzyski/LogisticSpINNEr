function seeds = make_relocation_replicate_seeds(baseSeed, replicateId)
%MAKE_RELOCATION_REPLICATE_SEEDS Derive deterministic component-specific seeds.

if ~isscalar(replicateId) || replicateId < 1 || ...
        replicateId ~= floor(replicateId)
    error('BlockRelocation:InvalidReplicate', ...
        'replicateId must be a positive integer.');
end
replicateOffset = 1000 * replicateId;
seeds.baseSeed = baseSeed;
seeds.dataSeed = baseSeed + replicateOffset + 11;
seeds.swapSeed = baseSeed + replicateOffset + 22;
seeds.responseSeed = baseSeed + replicateOffset + 33;
seeds.cvSeed = baseSeed + replicateOffset + 44;
end
