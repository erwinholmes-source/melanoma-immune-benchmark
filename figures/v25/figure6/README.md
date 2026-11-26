# Main Figure 6

Place this directory at `figures/v25/figure6/` under the Additional file 3 archive root. The four panels use the within-cohort held-out AUC contrasts, the two PUCH empirical ROC curves at the original ridge penalty of 1, and the paired external AUC and Brier changes at fixed penalties 0.1, 1 and 10.

## Recreate the figure from aggregate results

Run from the archive root after installing `requirements_figures.txt`.

```sh
python figures/v25/figure6/plot_figure6.py --archive . --output ../v25_figure_replay
```

The plot uses these files within the archive.

- `results/increment/increment_all_results.csv` gives the 24 held-out AUC contrasts and their 95% intervals for panel a. Intervals use 500 outcome-stratified patient bootstrap samples, each repeating preprocessing and fitting. APM uses seven cohort-available genes in Riaz, Hugo and Liu and six in Gide.
- `figures/v25/figure6/results/PUCH_roc_vertices.csv` gives each curve's aggregated false-positive and true-positive coordinates for panel b. `external_PUCH/extension_v21/results/PUCH_frozen_model_performance.csv` gives the reference AUC and Brier results for the original penalty.
- `external_PUCH/penalty_sensitivity/ridge_penalty_sensitivity.csv` provides the paired changes and interval limits in panels c and d. The three pairs of models were fitted in 249 primary-cohort patients and evaluated in the same 55 PUCH patients. Nominal 95% intervals use 5,000 outcome-stratified external patient resamples, conditional on the fitted models and response counts.

## Rebuild the empirical ROC curves

Place the two publicly sourced PUCH CSVs in a private input directory following `external_PUCH/PUCH_README.txt`. The projection uses the original model coefficients, gene membership and training-patient scaling from `external_PUCH/extension_v21/results/frozen_models.json`. The script verifies the frozen 733-gene universe and both PUCH source files. Run from the Additional file 3 root.

```sh
python figures/v25/figure6/build_roc.py --archive . --input-dir ../v19_PUCH_raw --output ../v25_roc_replay
python figures/v25/figure6/plot_figure6.py --archive . --roc ../v25_roc_replay/PUCH_roc_vertices.csv --output ../v25_figure_replay
```

The endpoint is CR/PR versus SD/PD for 14 responders and 41 nonresponders. At each distinct probability the algorithm adds the entire tied group, then retains the endpoints and direction changes on the empirical curve. The CSV contains only model names and cumulative rates. There are 25 vertices for IFN-4 and 20 for IFN-4 plus APM6. The frozen-model probabilities for the 55 external patients and the thresholds are used only in memory. AUC from the trapezoidal area is checked independently against the fraction of correctly ordered responder–nonresponder pairs, assigning half credit to ties. The script also verifies the Brier score on the same patients. `results/PUCH_roc_verification.json` records the checks. For the AUC and Brier sensitivities, the full source replay is documented in `external_PUCH/penalty_sensitivity/README.md`.
