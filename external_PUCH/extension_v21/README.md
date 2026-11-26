# PUCH external panel scores and fixed-coefficient models

This exploratory extension evaluates the five published immune gene sets in the 55 PUCH pretreatment biopsies. It also evaluates whether a coefficient learned for antigen-presentation expression improves a fixed, single-specimen rank model in PUCH. The original two PUCH score results and the four benchmark-cohort results were already available. The analysis decisions were written before computing these new PUCH endpoints. See `method_decisions.txt`.

## Data and cohort characteristics

The pinned authors' expression and clinical files, source scale, CR/PR versus SD/PD coding and original 55-patient reconciliation are documented in `../PUCH_README.txt`. Of the 55 patients, 14 responded and 41 did not. Table 3 of the PUCH study reports 24 acral, eight mucosal, 18 cutaneous and five unknown tumor sites. The public per-patient clinical CSV used here reports treatment time, cohort and response. It has no subtype field. Results in these files describe the entire mixed-subtype series. Source: [Cui et al., 2021, Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC7862369/#Tab3).

## Equal-weight fixed gene lists

`results/PUCH_panel_results.csv` contains the four-gene score, IFN-4 and five adapted published sets. The available version uses all genes of each published set measured in PUCH. All five sets had complete coverage in PUCH. The common version uses the exact membership previously defined across the four primary cohorts. The file records every included and omitted gene. Each gene was standardized on the 55 PUCH samples using its sample SD on the deposited log scale. Scores average signed standardized values, with CD276 negative in GEP. These adaptations compare gene contents and do not reproduce the original weighted assays.

The endpoints were standalone response AUC, the paired difference from the standalone IFN-4 AUC, and Spearman correlation with IFN-4. The percentile intervals for AUC and its paired difference use the same 5,000 outcome-stratified patient resamples throughout, conditional on these scores. Correlation intervals use 5,000 unstratified paired patient resamples and recompute midranks for each draw. All intervals are nominal. Four-gene and IFN-4 point estimates, AUC intervals and their paired difference match the preceding PUCH analysis to numerical precision.

## Single-specimen rank model

All 733 genes in the pre-existing historical reference universe were measured in PUCH. The platform inventory was inspected to choose a fully covered archived universe. Target-cohort expression values, normalization statistics and outcome labels were not used to set predictor centering, scaling, weights, the gene list or the intercept. Within each patient's biopsy, all 733 genes were ranked on the source expression scale. For each score, the mean midrank over its genes was transformed to (mean rank − 1)/732 − 0.5. IFN-4 comprises IFNG, IDO1, CXCL10 and HLA-DRA. APM used the six genes available in every source cohort (B2M, HLA-A, HLA-B, HLA-C, TAP1 and TAPBP). The single-specimen rank IFN-4 and APM6 are distinct predictors from the cohort-standardized, available-gene scores reported above.

Two patient-weighted logistic ridge models used all 249 response-evaluable patients (104 responders) from Riaz, Hugo, Liu and Gide. Each had one common unpenalized intercept, with L2 penalty 1 on slopes; there were no study indicators. The means and sample SDs of the two rank scores came from the pooled training patients. The models used IFN-4 alone and IFN-4 plus APM6. Coefficients, gene membership, predictor centers, predictor scales and the rank-universe checksum are in `results/frozen_models.json`. This file is written before PUCH response labels are read. No PUCH scores or outcomes are used for fitting, probability recalibration, threshold selection or choice between several model specifications. This transfer was chosen after the earlier PUCH four-gene/IFN-4 results and the original four-cohort APM results were known. The original cohort-wise cross-validated logistic models and this new pooled rank model use different input scales.

`results/PUCH_frozen_model_performance.csv` gives AUC and Brier for both models, and the paired joint-minus-IFN4 changes. Lower Brier is better. Percentile 95% intervals use the same 5,000 outcome-stratified paired PUCH resamples as the standalone AUC analysis. Inference conditions on the fitted model and observed number of responders; training-set sampling and future prevalence differences are outside these intervals. A single-cohort test of the fixed models does not determine general clinical utility.

## Reproduction

Place this directory at `external_PUCH/extension_v21/` in Additional file 3. From the package root, prepare the original-cohort workspace following the main README (`code/fetch_sources.py`, `code/prepare_sources.py`), and obtain the PUCH files with `external_PUCH/puch_external_check.py`. Then run:

```sh
python external_PUCH/extension_v21/puch_extension.py --package-root . --workspace ../v19_work --input-dir ../v19_PUCH_raw --output ../puch_extension_results
Rscript --vanilla external_PUCH/extension_v21/verify_puch_extension.R . ../v19_work ../v19_PUCH_raw ../puch_extension_results
```

The Python analysis needs NumPy and pandas. The optional independent R point-estimate replay uses base R and jsonlite. The Python runner checks both PUCH Git blob hashes, counts, the original two PUCH score results and complete reference coverage. It exports aggregate results and model parameters. The R replay starts from gene-expression and clinical files, computes midranks and AUC independently, and refits ridges with nlminb. Aggregated validation reports are in `results/verification.json` and `results/independent_R_verification.json`. No patient identifiers, expression profiles or individual predictions are distributed with this extension.
