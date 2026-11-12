# Figs. 5–6: bootstrap fits on the real data

**The data are not included.** They are subject-level resting-state fMRI connectivity matrices and cannot be
shared publicly. To run the scripts, place the two data files in `data/` (see the header of
`generate_bootstrap_results.m` for their format).

- **Sample:** n = 160 subjects (78 with and 82 without a family history of alcoholism), after one subject with
  missing values was dropped; 214 × 214 connectivity matrices (200 Schaefer cortical regions + 14 subcortical),
  without global signal regression, low-sucrose session.
- **Model:** Logistic SpINNEr with an intercept, W with ones off the diagonal and zeros on it, λ tuned by
  5-fold CV within each bootstrap sample (15 × 15 grid, `maxIters = 100`). Toolbox: `../software_84509f9/`.
- **Bootstrap:** 300 samples, drawn with replacement and standardized within each sample.

| Script | Role |
|---|---|
| `generate_bootstrap_results.m` | bootstrap fits. The 300 fits were run in five parts (35, 35, 30, 90 and 110 samples); this is the part-5 version, which draws and discards 190 samples before fitting 110. No seed was set. Writes `results/bootstrap_results_part<k>.mat` and `allResults_part<k>.mat`. |
| `summarize_results.m` | merges the parts into `results/EstimatesBootstrap.mat` (214 × 214 × 300) and makes exploratory clustered heatmaps of frequently selected edges |

The edges in Figs. 5–6 are those nonzero in at least 90% of the bootstrap fits, with their median coefficient,
restricted to the 200 cortical regions (41 edges). The figures were drawn with BrainNet Viewer (Fig. 5) and a
circular network graph (Fig. 6); the scripts for this last selection and drawing step are not included.
