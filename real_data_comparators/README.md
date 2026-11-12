# Real data: edgewise FDR tests and CPM as comparators for Logistic SpINNEr

Two standard connectome analyses of the real data (FH+ vs FH−), used in the Real Data section of the
manuscript to place the Logistic SpINNEr edges of Figs. 5–6 among established methods:

- **Edgewise tests with FDR control:** a two-sample t-test for each of the 22,791 edges, with
  Benjamini–Hochberg q-values (Benjamini & Hochberg, 1995).
- **Connectome-based Predictive Modeling (CPM)** (Shen et al., 2017), adapted to a binary outcome.

**The data are not included** (see [`../fig5_fig6_real_data/README.md`](../fig5_fig6_real_data/README.md)).

## Design

- **Data** (`load_aud_data.m`): the same data as the bootstrap analysis of Figs. 5–6, read from
  `../fig5_fig6_real_data/data/`: `GSRandNonGSRarrays.mat` (array `AA_NOGSR_low`: 214 × 214 × 161
  correlation matrices, no global signal regression, low-sucrose condition), `gen_hist_alcohol.csv`
  (`FHP`/`FHN` per subject) and, for network labels only, `Schaefer_parcellation.csv`. Subject 119 has
  missing values and is dropped, leaving n = 160 (78 FH+, 82 FH−); matrices are symmetrized with zero
  diagonals. The features are the 22,791 upper-triangle edges, Fisher z-transformed (atanh). No covariates.
- **Edgewise tests** (`edge_tstats.m`): pooled-variance two-sample t-test per edge (agrees with `ttest2`),
  two-sided p-values, Benjamini–Hochberg q-values; an edge is significant if q < 0.05.
- **CPM** (`cpm_fold.m`), within each training part of a cross-validation split:
  1. select the edges whose t-test p < 0.01 (main setting; 0.05 and 0.001 as sensitivity settings);
  2. split them into a positive network (higher in FH+) and a negative network (lower in FH+);
  3. sum each subject's edge values over each network (network strengths);
  4. fit a logistic regression of FH status on the two standardized strengths, which scores the test
     subjects. This replaces the linear regression of the original CPM, which is for a continuous outcome.

  Stratified 5-fold cross-validation, repeated 10 times (seed 20261022); a permutation test repeats the
  whole pipeline under 5000 label permutations (seed 20261113).
- **Comparison with Logistic SpINNEr** (`summarize_rda_comparators_final.m`): the edges nonzero in at least
  90% of the 300 bootstrap fits (`../fig5_fig6_real_data/results/EstimatesBootstrap.mat`) are located among
  the edgewise statistics (|t| rank, q-values, a hypergeometric test of their overlap with the edges with
  p < 0.001, agreement of signs) and in the CPM training folds (how often CPM selects each edge).

## How to run

Place the data files in `../fig5_fig6_real_data/data/` and create
`../fig5_fig6_real_data/results/EstimatesBootstrap.mat` with the scripts in that folder. Then, in MATLAB:

```matlab
cd real_data_comparators
run_rda_comparators_final([], true);  summarize_rda_comparators_final(true);  % smoke test, ~5 s
run_rda_comparators_final();          summarize_rda_comparators_final();      % full run, ~8 min
```

Needs the Statistics and Machine Learning Toolbox only (not the SpINNEr toolbox). Each result file is saved
when complete and skipped on re-run. The smoke test (30 regions, 2 repeats, 20 permutations) writes to
`results_smoke/` and only checks the pipeline. All settings are in `config_rda_comparators_final.m`.

## Outputs (`results/`)

| File | Contents | Included |
|---|---|---|
| `summary_fdr.csv` | number of edges tested, smallest p and q, edges with p < 0.001, FDR-significant edges, direction of the t-statistics | yes |
| `summary_overlap.csv` | the Logistic SpINNEr edges against both comparators | yes |
| `summary_spinner_edges.csv` | per Logistic SpINNEr edge: regions, Yeo networks, bootstrap frequency, median coefficient, t, p, q, \|t\| rank, CPM selection frequency | yes |
| `summary_cpm.csv` | CPM by edge threshold and model | regenerated |
| `fdr_results.mat`, `cpm_rep*.mat`, `cpm_perm.mat`, `cpm_outer_folds.mat` | per-edge statistics, per-fold CPM results, permutation nulls, outer partitions (contain subject-level labels) | no |
| `summary.txt` | everything the summary prints | regenerated |

In the original project, repeat 1 of the CPM splits was also checked against, and compared with, a
nested cross-validation of the four penalized methods that is not part of the manuscript; that step is
skipped when its files are absent.

## References

- Y. Benjamini and Y. Hochberg, "Controlling the false discovery rate: a practical and powerful approach to
  multiple testing," *J. R. Stat. Soc. B*, vol. 57, no. 1, pp. 289–300, 1995.
- X. Shen et al., "Using connectome-based predictive modeling to predict individual behavior from brain
  connectivity," *Nat. Protoc.*, vol. 12, no. 3, pp. 506–518, 2017.
