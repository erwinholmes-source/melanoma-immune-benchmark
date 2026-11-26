PUCH: separate revision-stage external check

Purpose and boundary
The PUCH candidate was screened against a dated eligibility specification before inspecting its scores or performance. Its expression/clinical intersection must reproduce the original 55-patient population, pretreatment anti-PD-1 monotherapy setting, objective response labels and complete four-gene/IFN-4 coverage. It is a revision-stage external check, not an additional member of the original four-cohort retrospective benchmark. There is no data-dependent pooling, gene reselection, cut-point selection, or incremental-model fitting.

Primary sources
Cui C, Xu C, Yang W, et al. Ratio of the interferon-gamma signature to the immunosuppression signature predicts anti-PD-1 therapy response in melanoma. npj Genomic Medicine. 2021;6:7. DOI: 10.1038/s41525-021-00169-w.
https://www.nature.com/articles/s41525-021-00169-w
https://pmc.ncbi.nlm.nih.gov/articles/PMC7862369/
Author-linked repository: https://github.com/xmuyulab/ims_gene_signature
Pinned commit: 275be3c340e4cba698d085bc1a875b92dd45602c
Expression: https://github.com/xmuyulab/ims_gene_signature/blob/275be3c340e4cba698d085bc1a875b92dd45602c/data/mel_puch_exp_data.csv
Clinical: https://github.com/xmuyulab/ims_gene_signature/blob/275be3c340e4cba698d085bc1a875b92dd45602c/data/mel_puch_cli_data.csv
Loader: https://github.com/xmuyulab/ims_gene_signature/blob/275be3c340e4cba698d085bc1a875b92dd45602c/src/load_data.py
Response handling: https://github.com/xmuyulab/ims_gene_signature/blob/275be3c340e4cba698d085bc1a875b92dd45602c/src/calculate_performance.py

Scale and labels
The publication's Patients and specimens section describes 55 pretreatment FFPE tumor specimens from melanoma patients treated with anti-PD-1 monotherapy, assessed by RECIST 1.1. Its RNA-sequencing methods describe RSEM gene-level FPKM and a log2 transformation. The repository supplies processed log-scale expression; its loader adds log2(x+1) only to matrices whose maximum exceeds 500. The deposited PUCH matrix is therefore used directly, with no extra log transformation. The loader docstring also mentions quantile normalization; a raw-to-processed normalization reconstruction is not claimed. Housekeeping normalization in the original paper pertains to its own signature scores and is not added to the fixed z-score rule here.
The repository labels are 1 = CR/PR, -1 = SD and 0 = PD. The original authors excluded SD from their ROC plots; this check retains SD as nonresponse to match the benchmark's objective-response definition. Therefore these AUCs should not be directly compared with the original paper's SD-excluding AUC.

Sample reconciliation
The clinical table contains 60 unique pretreatment records (15 CR/PR, 6 SD, 39 PD). The expression matrix contains 55 unique sample columns. All 55 match clinical records, giving exactly 14 CR/PR, 6 SD and 35 PD, consistent with the published population. Five clinical records have no expression column (1 CR/PR, 4 PD); no manual sample removal is used. Matching is by the original identifiers in private execution, preserving deposited matrix column order. Identifiers and individual values are not included in this release.

Fixed computation
Four-gene score: CD7, CXCL9, LYZ, STAT1.
Nonoverlapping IFN-4: IFNG, IDO1, CXCL10, HLA-DRA.
For each gene, z-standardize over all 55 baseline patients without using response labels, using sample SD; average the four positively signed standardized values. Report Spearman correlation; AUC for CR/PR versus SD/PD; and paired four-gene minus IFN-4 AUC difference. AUC gives ties half credit. The 95% percentile intervals use 5,000 outcome-stratified bootstrap samples, seed 20261123, with the same resamples for both fixed scores. Bootstrap intervals are conditional on these cohort-standardized scores. No equivalence claim is made from a CI crossing zero.

Reproduction
Requires Python 3, numpy and pandas. From the supplement root folder:
python external_PUCH/puch_external_check.py --download --input-dir ../v19_PUCH_raw --output-dir ../v19_PUCH_results
The optional downloader stores authors' files only in the specified input directory. Do not redistribute that directory with this public aggregate release. If a network transfer fails, remove the unfinished .partial file and retry, or download the pinned source files directly. The script checks exact Git blob hashes and aborts if the expression file is incomplete or altered, if samples/labels fail reconciliation, or if target-gene coverage is incomplete. The JSON records SHA-256 hashes and software-independent computation settings. It contains no patient identifiers or individual predictions.

Interpretation limit
This is a small retrospective external dataset acquired during revision. It provides an additional same-population comparison of two fixed expression summaries. It does not establish prospective clinical utility or validate an incremental predictive model beyond IFN-4.
