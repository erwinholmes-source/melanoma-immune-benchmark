# Training-cohort diagnostics and PUCH rank sensitivity

The negative APM-6 coefficient in the original joint model prompted an exploratory analysis of how its estimate depends on IFN-4 adjustment, cohort intercepts and cohort weights. All fits used the 249 response-evaluable training patients (104 responders). The study-specific fits share slopes across cohorts. The primary external test uses the original frozen model. A subsequent, exploratory ranking experiment applied the slopes from one prespecified diagnostic configuration to PUCH.

## Features and models

For each biopsy, the analysis ranked the same 733 genes, averaged tied ranks, calculated centered rank percentiles and averaged the four IFN genes or six APM genes. Both scores were standardized with the original pooled training mean and sample SD. All configurations, including omissions, used these same location and scale parameters. Panel members, full-precision scaling and reference-model coefficients are recorded in `../extension_v21/results/frozen_models.json`.

For each configuration we fitted IFN-4 alone, APM-6 alone and their joint model. The objective was the sum of weighted patient logistic losses plus half the sum of the squared slopes (λ = 1). All intercepts were unpenalized. Each patient had either unit weight or weight 249/(4 × cohort size), giving each cohort 62.25 units in the balanced analysis. Total weight was 249 in both designs. We crossed both weighting designs with a single pooled intercept or four unrestricted cohort-specific intercepts. Each slope is measured per one original pooled training SD. The four cohort-specific intercepts are nuisance parameters for the observed training cohorts.

| Intercept | Weights | APM-6 slope alone | APM-6 slope with IFN-4 |
| --- | --- | ---: | ---: |
| Common | Equal per patient | −0.024 | −0.188 |
| Common | Equal per cohort | −0.066 | −0.166 |
| Cohort specific | Equal per patient | +0.615 | +0.328 |
| Cohort specific | Equal per cohort | +0.666 | +0.528 |

The APM-6 slope changed sign when the analysis allowed cohort-specific intercepts under either weighting scheme. Equalizing the total weights of the cohorts while retaining a common intercept left the joint slope negative. For example, the Hugo mean APM-6 score was 2.40 training SDs below the pooled mean. Within Hugo, the standalone APM-6 score had AUC 0.702 (14 responders among 26 patients). The Riaz mean APM-6 score was 0.54 SDs above the pooled mean, with AUC 0.495 (10 responders among 49 patients). Cohort differences in assay, case mix, treatment or response ascertainment may contribute to shifts in score distributions or response rates. Their contributions cannot be separated in this analysis.

## Output files

- `cohort_rank_summaries.csv` gives sample and responder counts, the AUC of each prespecified rank score in its original direction, their Spearman and Pearson correlations, and their mean standardized feature values. The pooled row includes responder–nonresponder pairs from different cohorts. The individual cohort rows make comparisons only within the named cohort. These are descriptive score AUCs on the transfer feature scale.
- `training_model_coefficients.csv` lists the common slopes and applicable intercepts for all 12 fitted models. Empty fields identify parameters outside the given specification.
- `leave_one_cohort_out.csv` reports four additional common-intercept joint fits, each omitting one training cohort. Patients in the retained cohorts had equal relative weights rescaled to sum to 249, holding the penalty scale fixed. The four estimated APM-6 slopes ranged from −0.435 to −0.048. These are influence checks on training coefficients.
- `verification.json` records differences between the independent R and Python calculations, the numbers of estimates checked, and the verification of archived input hashes.

## Post hoc PUCH ranking sensitivity

The change in APM slope with a cohort-specific intercept motivated one further external comparison. After the training diagnostics and the original PUCH results were known, we recorded the choice of models and endpoints in `rank_transfer_design.txt`. We used the standardized IFN-4 slope from its single-predictor patient-weighted, cohort-intercept model (0.483). The corresponding joint model used the IFN-4 and APM-6 slopes 0.403 and 0.328. Each model used its full-precision slopes to rank all 55 PUCH patients. Adding a common constant to a linear score preserves the ranking. Absolute risk predictions would require an intercept appropriate to the external population.

| PUCH ranking or contrast | AUC or ΔAUC | 95% interval |
| --- | ---: | ---: |
| IFN-4 | 0.765 | 0.604 to 0.906 |
| IFN-4 plus APM-6 | 0.780 | 0.636 to 0.906 |
| Joint minus IFN-4 | +0.016 | −0.033 to +0.066 |

The paired increment changed direction relative to the original common-intercept transfer (−0.017). Both paired intervals included zero. The external ranking sensitivity estimates discrimination for one alternative set of slopes. The 95% percentile intervals use the original seed (20261123) and 5,000 paired outcome-stratified resamples of the same 14 responders and 41 nonresponders. These intervals condition on the fitted slopes and observed outcome counts.

`external_cohort_intercept_rank_auc.csv` contains the two AUCs and their paired difference. `external_rank_verification.json` compares pairwise AUC with the Mann–Whitney U statistic for the original patients and every bootstrap resample. `external_rank_R_verification.json` verifies the two point AUCs and their difference after reconstructing ranks independently in R.

## Reproduction

From the unpacked Additional file 3, prepare the four-cohort workspace following the root README. Pass its path in place of `../v19_work` if needed. R requires jsonlite. Python requires NumPy and pandas. Set the PUCH input directory to the folder with the two pinned source files downloaded by the original PUCH reproduction instructions. The first two commands read only the training cohorts, and the last two apply the resulting slopes to PUCH.

```sh
Rscript --vanilla external_PUCH/training_diagnostics/training_diagnostics.R . ../v19_work ../training_diagnostic_results
python external_PUCH/training_diagnostics/verify_training_diagnostics.py --package-root . --workspace ../v19_work --results ../training_diagnostic_results
python external_PUCH/training_diagnostics/external_rank_sensitivity.py --package-root . --input-dir ../v19_PUCH_raw --results ../training_diagnostic_results
Rscript --vanilla external_PUCH/training_diagnostics/verify_external_rank.R . ../v19_PUCH_raw ../training_diagnostic_results
```

R computes the Mann–Whitney rank AUC and fits the weighted ridge models using nlminb. The Python verification uses the archived feature construction, counts pairwise wins with half credit for ties, and solves the weighted score equations by Newton iteration. It checks the SHA-256 of the original eight training files, HGNC mapping and rank universe. Each program also checks the frozen training scaler. The R script checks the reproduced original intercepts and slopes against the archived full-precision model. The output consists entirely of aggregate estimates.
