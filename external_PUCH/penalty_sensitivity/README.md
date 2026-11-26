# PUCH ridge penalty sensitivity

This exploratory sensitivity analysis repeats the transfer experiment with fixed ridge penalties of 0.1, 1 and 10. The original analysis fixed λ = 1. The response labels, 249 training patients (104 responses), single-sample ranks over the archived 733 genes, training means and sample standard deviations were identical across fits. Each fit minimized the **sum** of individual logistic losses plus λ/2 times the sum of squared slopes. Patients had equal weights, and the common intercept was unpenalized. All six models were fitted before evaluation on the 55 PUCH patients (14 responders). The comparison retained the original λ = 1 model throughout.

At λ = 1, the IFN-4-only model has intercept −0.354869 and standardized IFN-4 slope 0.541680. The joint model has intercept −0.354228, IFN-4 slope 0.591873 and APM-6 slope −0.188319. The APM coefficient is conditional on IFN-4 in the pooled training samples. The full-precision scaling and coefficients are in `../extension_v21/results/frozen_models.json`. The common APM genes are B2M, HLA-A, HLA-B, HLA-C, TAP1 and TAPBP. TAP2 was unavailable in the Gide analysis matrix and present in the PUCH matrix.

`ridge_penalty_sensitivity.csv` contains external AUC and Brier estimates and 95% percentile intervals for each model and their paired differences (joint minus IFN-4). The script uses the same 5,000 stratified PUCH patient bootstrap draws as the original transfer analysis. Intervals condition on the fitted models and observed outcome counts. At λ = 0.1, 1 and 10 respectively, the paired AUC differences were −0.0174, −0.0174 and −0.0192. The Brier differences were +0.000707, +0.000604 and −0.000025. Each interval for a difference included zero.

## Reproduction

From the unpacked Additional file 3, prepare the four-cohort workspace and download the pinned PUCH expression and clinical CSVs as described in `external_PUCH/extension_v21/README.md`. With NumPy and pandas installed, run the following command, adjusting the two input directory paths if needed.

```sh
python external_PUCH/penalty_sensitivity/ridge_sensitivity.py --package-root . --workspace ../v19_work --input-dir ../v19_PUCH_raw --output ../ridge_penalty_sensitivity.csv
```

The script checks the training file hashes, external input hashes, scaler and λ = 1 coefficients. It reproduces the archived λ = 1 estimates and all interval endpoints to absolute tolerance 1e−12 and verifies the gradient of the summed objective for each fit. Reported results are aggregate estimates.

## Primary melanoma site in PUCH

The original Table 3 reports 24 acral, eight mucosal, 18 cutaneous and five unknown sites across the 55 patients. Its Supplementary Information (one gene table, six figures) and Reporting Summary contain no patient identifier paired with a primary melanoma subtype. The index of the analysis-ready clinical CSV is the sample identifier. Its three columns specify dataset, response and sampling time. The separate survival CSV lists status and overall survival. The expression matrix provides gene expression by sample. [Cui et al., Table 3 and supplementary information](https://www.nature.com/articles/s41525-021-00169-w) and [pinned public analysis repository](https://github.com/xmuyulab/ims_gene_signature/tree/275be3c340e4cba698d085bc1a875b92dd45602c/data).

The public HRA000524 metadata workbook lists 59 individual/sample records and contains two sheets, “Individuals & samples” and “Files.” The sample sheet lists individual and sample identifiers, gender, sample name and sample description. The source also has sample detail pages with a tissue field (for example, “Foot”). That field identifies the sampled tissue. Primary-site classification requires the site of origin, including for biopsies taken from metastases. The individual detail page lists identifiers, gender and release date. The raw archive is under controlled access. [GSA-Human HRA000524](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000524) and [example sample record](https://ngdc.cncb.ac.cn/gsa-human/browse/sampleDetail/HRS123850).

Two further studies from the same center describe this 55-patient cohort. The 2022 study gives aggregate primary-site counts in Table 1. Its four supplementary tables concern gene pairs, training analyses, TCGA patients and published signatures. The 2024 study reports cohort-level baseline characteristics in Supplementary Table S5. None of these materials supplied a verified link between individual expression profiles and primary melanoma subtype. Subtype-specific performance was therefore not estimated. [2022 PUCH study](https://doi.org/10.3389/fimmu.2022.839901) and [2024 PUCH study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10976329/).
