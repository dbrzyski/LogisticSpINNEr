# Figs. 2–3: block size and relocation of the support

Two simulation packages and one plotting script. Both packages use n = 300, p = 60, three blocks with signs
(+, −, +) and coefficient magnitude γ = 0.5, intercept 1, base seed 20250705, 50 replicates, and the four
methods of the manuscript with 5-fold CV on 10 × 10 λ grids and `maxIters = 100`. Error is the relative
Frobenius error. The two packages share n, p, the seed, the design generator and the block positions, so they
can be combined in Fig. 2.

| Folder | Contents | Manuscript |
|---|---|---|
| `Edge_Relocation_Gamma_0_5/` | blocks (5,5,4); k = 0, …, 6 of the 26 nonzero upper-triangular edges relocated | Fig. 2 (first point), Fig. 3(a) |
| `Block_Relocation_Gamma_0_5/` | blocks (5,5,5), (6,6,6), (7,7,7), (8,8,8), (12,12,12), (16,16,16); relocated proportion 0, …, 6/26 | Fig. 2, Fig. 3(b)–(d) |
| `plotting/make_gamma05_figures.m` | Figs. 2 and 3 from the two summary CSVs in `*/results/` | writes `plotting/figures/` |

## How to run

Every task is seeded, checkpointed to `results/` and skipped if its result file exists, so runs can restart.
The original runs used a Slurm array on a cluster, one MATLAB task per call below.

```matlab
% Edge relocation: one call per replicate r = 1..50
cd fig2_fig3_relocation/Edge_Relocation_Gamma_0_5
for r = 1:50, run_swap_pilot_replicate(r, 10); end      % second argument: number of parallel workers
aggregate_swap_pilot_results();                         % writes recomputed_summaries/

% Block relocation: tasks 1..300 (block setting varies fastest; tasks 1:6 are replicate 1)
cd ../Block_Relocation_Gamma_0_5
for t = 1:300, run_block_relocation_task(t, 10); end
aggregate_block_relocation_results();                   % writes recomputed_summaries/

% Figures, from the summary CSVs in */results/
cd ../plotting
make_gamma05_figures
```

The summary CSVs of the manuscript runs are included in `*/results/`, so the plotting script runs without
refitting. To plot a rerun, copy the CSVs from `recomputed_summaries/` into `results/`.

The manuscript uses `fig2-intact-block-gamma05.pdf` and the one-row Fig. 3, `fig3-row-gamma05_a.pdf` … `_d.pdf`
(blocks (5,5,4), (6,6,6), (7,7,7) and (12,12,12)) with `fig3-row-gamma05_legend.pdf`.

## Toolbox version

The original runs used a copy of the SpINNEr toolbox downloaded in August 2026 whose commit was not recorded.
The scripts point to `../software/`. Fits of Logistic NUCLEAR, and through the λ grids possibly of the other
logistic methods, can differ slightly if the download predates commit `3f07a54` (see the top-level README).
