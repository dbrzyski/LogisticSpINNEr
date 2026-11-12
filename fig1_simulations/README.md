# Fig. 1: estimation error in simulation

Data model (all panels): p = 60; B has three blocks, 5 × 5 at +γ, 5 × 5 at −γ and 4 × 4 at +γ; intercept 1;
symmetric Gaussian matrices z-scored entrywise with zero diagonals; y ~ Bernoulli(logistic(⟨A, B⟩ + 1)).
The four methods are fitted with 5-fold CV on 10 × 10 λ grids and `maxIters = 100`. Error is the relative
Frobenius error ‖B̂ − B‖_F / ‖B‖_F. Toolbox: `../software_84509f9/`.

| Panel | Script | Setting | Output |
|---|---|---|---|
| (a)–(e) heatmaps | `a-e_heatmaps/getEstimates.m` | one replicate, n = 300, γ = 0.5; no seed | `a-e_heatmaps/output/B_ESTIMS.mat` and heatmap PDFs (the manuscript shows the 35 × 35 `*Red.pdf` crops) |
| (f) sample size | `f_sample_size/run_experiment.m` | n = 100, 120, …, 240, γ = 1, 50 replicates, `rng(2)` | `f_sample_size/results.mat` (absolute Frobenius error) |
| (g) signal strength | `g_signal_strength/run_experiment_new_structure.m` | n = 300, γ ∈ {0.1, 0.25, 0.5, 1, 1.5, 2}, 50 replicates; no seed | `results_sigStr_sweep.mat` |
| (f), (g) plots | `f-g_plots/generatePlots.m` | reads the two result files | `f-g_plots/output/*.pdf` |

The result files used for the manuscript are included, so `generatePlots.m` runs without refitting. The
signal-strength run saves `results_sigStr_sweep.mat`; the included copy was renamed to
`results_sigStr_sweep_relative_iter50.mat`, which is the name `generatePlots.m` reads. `generatePlots.m` divides
the sample-size errors by ‖B‖_F; the signal-strength errors are already relative.

Each panel used its own independent set of replicates, so results for the same condition can differ slightly
between panels. Each run fits 4 methods × 50 replicates and takes several hours.
