function seeds = make_replicate_seeds(baseSeed, replicateId)
%MAKE_REPLICATE_SEEDS Derive separate deterministic seeds for one replicate.

replicateOffset = 1000 * replicateId;
seeds.baseSeed = baseSeed;
seeds.dataSeed = baseSeed + replicateOffset + 11;
seeds.swapSeed = baseSeed + replicateOffset + 22;
seeds.responseSeed = baseSeed + replicateOffset + 33;
seeds.cvSeed = baseSeed + replicateOffset + 44;
end
