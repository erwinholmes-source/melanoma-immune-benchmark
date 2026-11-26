# Melanoma immune-signature benchmark

This archive accompanies manuscript v29 (BioData Mining). It contains analysis code, gene lists, protocols and aggregate results. Riaz, Hugo, Liu and Gide contribute 253 pretreatment patients, of whom 249 had an evaluable response. Yale contributed to historical selection. PUCH provides a separate 55-patient transfer cohort. TCGA-SKCM (`skcm_tcga_pan_can_atlas_2018`, 367 metastatic melanomas) supplied expression profiles for defining an independent inflammatory reference gene pool.

## Evaluate a custom gene set

After preparing the four melanoma cohorts with the source commands below, run a fixed, equal-weight, positively signed panel. The optional PUCH directory enables external evaluation. For example, the common six-gene APM list is:

```sh
python run_benchmark.py --workspace ../v19_work --puch-dir ../v19_PUCH_raw \
    --genes B2M HLA-A HLA-B HLA-C TAP1 TAPBP --output ../APM6_benchmark.json
```

Use distinct HGNC symbols from the fixed 733-gene rank universe in `results/reference/universe_historical733.csv`, up to 242 genes per panel (the size of the TCGA inflammatory pool). The runner requires the same gene list in every cohort. For the manuscript's signed GEP analysis or its per-cohort available-gene analyses, use the original entry points below.

The reference output for the example above is in `examples/APM6_expected.json`. The JSON contains aggregate results for 10,000 size-matched TCGA draws, leave-one-cohort-out reference percentiles, paired repeated fivefold CV point estimates (AUC and Brier), and rank-model coefficients estimated with shared or cohort-specific intercepts. With PUCH supplied it also reports the paired external changes. For the cohort-specific intercept model the external output is AUC, because the external intercept is unknown. This custom-panel runner yields point estimates. The original study-specific bootstrap intervals and independent R checks remain available below. Keep the source matrices and generated personal workspaces outside the repository.

## Manuscript and archive guide

| v29 item | Where to find the underlying results |
|---|---|
| Fig. 1, design and cohort roles | `code/plot_figures_v24.py` |
| Fig. 2 and Table 3, five published panels measured with the same genes in five cohorts | `figures/v24/fig2_values.csv` gives the correlations, baseline counts, AUCs and intervals. `primary_panels/primary_panel_table.csv` supplies four-cohort summaries. The `common` and `available` rows in `external_PUCH/extension_v21/results/PUCH_panel_results.csv` distinguish shared and complete PUCH gene lists. |
| Fig. 3a, the four-gene comparator in four reference pools | `results/reference/` contains the original three pools. `external_TCGA_background/` contains the 242-gene TCGA-defined pool and the analyses excluding the four comparator genes. |
| Fig. 3b and Table 3, five published panels against two inflammatory gene pools | `tcga_panels/results/panel_reference_comparison.csv` contains both pool definitions, including sign-matched and all-positive GEP draws. `figures/v24/fig5_values.csv` lists the plotted values. |
| Post hoc leave-one-cohort-out and patient-weighted reference sensitivity (Table S37) | `code/leave_one_cohort_out_reference.py` and `results/reference/leave_one_cohort_out_reference.tsv` recompute reference percentiles across all 60 panel × pool × omission/weighting combinations using the fixed 10,000 matched draws. |
| Post hoc audit of Ayers expanded immune (18 genes, Table 2) and Davoli (7 genes) | `code/audit_additional_panels.py`, `results/reference/posthoc_additional_panels.csv` and `results/reference/posthoc_additional_panels_cohorts.csv` record the post hoc standalone concordance, cohort AUCs and size-matched inflammatory reference percentiles for the two compact panels omitted from the initial five-panel run. |
| Fig. 4a and Table 4, six additions to IFN-4 | `results/increment/increment_all_results.csv` gives paired AUC and Brier changes and intervals for the six gene sets in each primary cohort. |
| Fig. 4b, original frozen models in PUCH | `figures/v28/figure4/results/PUCH_roc_vertices.csv` gives empirical ROC curve coordinates at λ = 1. `external_PUCH/extension_v21/results/frozen_models.json` records the training coefficients and scaling. `PUCH_frozen_model_performance.csv` in the same folder gives external AUC and Brier estimates and paired changes. |
| Fig. 4c–d, external ΔAUC (ridge penalties and post hoc cohort-intercept sensitivity) and ΔBrier | The `joint_minus_IFN4` rows in `external_PUCH/penalty_sensitivity/ridge_penalty_sensitivity.csv` give the paired external ΔAUC and ΔBrier at λ = 0.1, 1 and 10; `external_PUCH/training_diagnostics/external_cohort_intercept_rank_auc.csv` supplies the post hoc cohort-intercept ΔAUC plotted as the open diamond in Fig. 4c. |
| Table 5, training intercepts and slopes | `external_PUCH/training_diagnostics/training_model_coefficients.csv` gives model coefficients for combinations of training-cohort intercepts, patient/cohort weights and predictors at λ = 1 (including both standalone and conditional APM fits). `external_PUCH/penalty_sensitivity/ridge_penalty_sensitivity.csv` gives penalty sensitivities, and `leave_one_cohort_out.csv` gives training-cohort omission fits. |
| Figs. 5–6, cellular sources and compartment comparisons | Carried forward from the single-cell module and Addendum 08 of the linked v16 Zenodo deposit and documented in Additional file 1. |
| Additional file 4 | Page 1 shows the earlier gene-set and IFNG-6 comparison (cross-referenced to main Fig. 3a–b and main Table 3). Page 2 shows precision and the four-gene incremental analysis. Page 3 shows the four-gene comparator and IFN-4 across cohorts. Page 4 contains the paired discrimination differences and exploratory equivalence margins for the four-gene score versus eight comparators. |

Gene lists in Fig. 2 and the size-matched results in Fig. 3b use 19 ICR, 17 GEP, 8 chemokine, 8 TLS and 6 APM genes shared across cohorts. The GEP comparator assigns one negative sign per random set in both inflammatory pools. The separately archived complete-PUCH fits and available-gene incremental models retain their original membership. Historic filenames inside older subfolders preserve their original manuscript numbering.

## Reproduce the manuscript figures from aggregate results

Run these commands from the extracted archive root. Write regenerated figures and tables to a sibling directory.

```sh
python -m pip install -r requirements_figures.txt
python code/plot_figures_v24.py --archive . --output ../v28_replay/figures_raw --tables ../v28_replay/figure_data
python figures/v28/figure4/plot_figure4.py --archive . --output ../v28_replay/figures
```

`code/plot_figures_v24.py` exports Figs. 1 and 2 (`Fig1.*`, `Fig2.*`) and v28 Fig. 3 (written under its historic v24 filename `Fig5.*`). `figures/v28/figure4/plot_figure4.py` exports v28 Fig. 4 (`Fig4.pdf` and `Fig4.png`, including the post hoc cohort-intercept sensitivity in panel c). Historic plotters under `figures/v25/figure6/` and `figures/v26/figure6/` are retained for archive provenance; `figures/v25/figure6/plot_figure6.py` omits the panel c cohort-intercept row added in v26–v28. Figs. 5 and 6 carry forward the cellular analyses from the v16 deposit.

## Reproduce the model diagnostics, external ROC curves and post hoc panel audit

First prepare the four-cohort training workspace and obtain the two PUCH source CSVs using the original-source instructions below. From the archive root, run the following commands. The training commands independently reconstruct the training features and fits in R and Python. The external ranking commands apply the cohort-intercept model slopes to the 55 PUCH patients. The ROC and figure commands replay v28 Fig. 4, and the audit command reproduces the post hoc Ayers expanded immune and Davoli summaries.

```sh
Rscript --vanilla external_PUCH/training_diagnostics/training_diagnostics.R . ../v19_work ../training_diagnostic_results
python external_PUCH/training_diagnostics/verify_training_diagnostics.py --package-root . --workspace ../v19_work --results ../training_diagnostic_results
python external_PUCH/training_diagnostics/external_rank_sensitivity.py --package-root . --input-dir ../v19_PUCH_raw --results ../training_diagnostic_results
Rscript --vanilla external_PUCH/training_diagnostics/verify_external_rank.R . ../v19_PUCH_raw ../training_diagnostic_results
python figures/v28/figure4/build_roc.py --archive . --input-dir ../v19_PUCH_raw --output ../v28_roc_replay
python figures/v28/figure4/plot_figure4.py --archive . --roc ../v28_roc_replay/PUCH_roc_vertices.csv --output ../v28_replay/figures
python code/audit_additional_panels.py --archive . --workspace ../v19_work --output ../v29_posthoc_audit
python code/leave_one_cohort_out_reference.py . ../v29_posthoc_audit/leave_one_cohort_out_reference.tsv
```

The diagnostic analysis followed the original PUCH results. The choice to test the equal-patient-weight, cohort-specific-intercept slopes was recorded after examining the training diagnostics in `external_PUCH/training_diagnostics/rank_transfer_design.txt`. `external_PUCH/training_diagnostics/README.md` documents the training coefficients and external rank estimates.

## Reproduce the inflammatory reference and original transfer analyses

First prepare the original-cohort workspace with the source reconstruction commands below. Follow `external_TCGA_background/README.md` to rebuild the TCGA gene-eligibility file from the pinned TCGA expression and sample tables. The two-pool runner reads the archived eligibility file. Obtain the two original clinical/expression CSVs for PUCH as documented below. Then run the following commands from the archive root.

```sh
python tcga_panels/run_panel_reference.py --archive . --workspace ../v19_work --output ../v24_panel_replay
Rscript --vanilla tcga_panels/verify_draws.R ../v19_work/prepared_v7 ../v24_panel_replay
python external_PUCH/extension_v21/puch_extension.py --package-root . --workspace ../v19_work --input-dir ../v19_PUCH_raw --output ../v24_puch_replay
Rscript --vanilla external_PUCH/extension_v21/verify_puch_extension.R . ../v19_work ../v19_PUCH_raw ../v24_puch_replay
python external_PUCH/penalty_sensitivity/ridge_sensitivity.py --package-root . --workspace ../v19_work --input-dir ../v19_PUCH_raw --output ../v24_puch_replay/ridge_penalty_sensitivity.csv
```

See `tcga_panels/README.md` for the draws, sign assignments, nominal tail probabilities and independent R calculations. `external_PUCH/penalty_sensitivity/README.md` documents the three fixed penalties. The original PUCH comparison was recorded before the expression profiles were accessed. The five-score external evaluation, APM transfer, cohort-intercept sensitivity and two-panel post hoc audit were specified after earlier study results were available. The paired PUCH intervals condition on the fitted models.

## Background analyses and initial PUCH screen

| Manuscript item | File |
|---|---|
| Table 4: five published gene-set adaptations added to IFN-4 | `results/increment/table4_incremental_prediction.csv` |
| Original four-gene comparison reproduced with the same code | `results/increment/four_gene_reproduction.csv` |
| All 24 cohort-by-panel comparisons, complete AUC/Brier intervals | `results/increment/increment_all_results.csv` |
| Three random backgrounds: 15,561 / 733 / 142 genes | `results/reference/background_comparison_summary.csv` |
| Scoring-rule sensitivity: available-gene z, common-gene z, common-gene rank | `results/reference/scoring_auc.csv` |
| Paired rank-minus-z sensitivity | `results/reference/rank_sensitivity_paired_differences.csv` |
| Separate PUCH check: 55 patients, 14 responders | `external_PUCH/PUCH_aggregate_results.csv` |
| Original PUCH eligibility and analysis record | `external_PUCH/PUCH_external_screen_spec.txt` |
| Common gene membership, omissions and negative signs | `results/reference/panel_coverage.csv` and `results/increment/gene_coverage.csv` |

The four primary RNA-seq cohorts contribute 49, 26, 102 and 72 response-evaluable patients (Riaz, Hugo, Liu and Gide). The 104 evaluable Yale patients are excluded from the new incremental analysis; their baseline expression is used only to reconstruct the historical 733- and 142-gene backgrounds.

All 20 published-set increments had approximate 95% AUC and Brier intervals containing zero. For APM, ΔAUC was 0.077 in Riaz, 0.256 in Hugo, −0.022 in Liu and −0.007 in Gide, with wide intervals in the first three cohorts. The APM random-set comparison used the six genes common to all four cohorts, whereas the incremental comparison used the available 7/7/7/6 genes in Riaz/Hugo/Liu/Gide. An extreme random-set percentile therefore did not establish a consistent conditional gain. The four-gene score ranked at the 99.72nd, 96.32nd, 86.95th and 90.60th percentiles against the broad, historical, original inflammatory and TCGA-defined inflammatory pools, respectively.

## Check the archive without patient data

Use Python 3.12 and the packages in `requirements.txt`:

```sh
python code/check_integrity.py
python code/check_aggregate_results.py
python code/check_v24_results.py
python code/check_v25_results.py
```

The integrity check covers every file listed in `FILE_MANIFEST.json`. The second command checks 346 interval and summary quantities. The third recomputes the published-gene-set comparisons against two inflammatory references and checks external-model penalty sensitivities. The fourth checks the areas under the two empirical ROC curves, the training summaries and coefficients, and the paired AUC contrast for the alternative external rankings. The numerical replay is documented in `validation/v25_replay.json`; earlier replays remain in `validation/`.

## Recalculate from original processed sources

Obtain and extract the v16 public code archive cited in the manuscript, [Zenodo record 24001782](https://doi.org/10.5281/zenodo.24001782). Its root folder is `Melanoma_immune_benchmark_v16_public`. The v19 reference runner reads its published A07/A10 aggregate references for reproduction checks. It does not expect private patient-score files or missing patient-level inputs inside that archive.

The nine required original processed expression/annotation files, URLs and exact source hashes are listed in `resources/source_files.json`. Download them locally under their original access/reuse conditions, either manually with those filenames or using the explicit download command below. Hash mismatches stop preparation; an updated file is not silently substituted. The frozen HGNC snapshot and 755-gene mapping are included as gene-level resources, so today's changing HGNC file is not required.

Run these commands from this supplement's root, using separate private work directories:

```sh
python -m pip install -r requirements.txt
python code/fetch_sources.py --destination ../v19_raw
python code/prepare_sources.py --raw-dir ../v19_raw --workspace ../v19_work
python code/run_increment_public.py --workspace ../v19_work --output ../v19_increment_results --workers 4
python code/run_reference_public.py --workspace ../v19_work --v7-root ../v19_work/prepared_v7 --v16-root ../Melanoma_immune_benchmark_v16_public --output ../v19_reference_results
```

Preparation is executable source parsing, not a request for unavailable author intermediates. It rebuilds the baseline populations and historical source layout, including score-free clinical tables and compact expression matrices. It preserves original baseline order, uses source CR/PR versus SD/PD labels, excludes on-treatment samples, averages repeated baseline lesions within patients, and applies each specified log transform once. It does not perform alignment or quantification of raw sequencing reads. The five rebuilt populations and original compact matrices were checked locally against the archived inputs; no previous patient scores were used for reconstruction. The recorded replay used cached copies of the original source files matching all nine hashes, rather than a new network-download availability test.

The incremental runner rebuilds its compact gene matrices and regenerates all 500 outcome-stratified patient bootstrap weight vectors from the documented seed sequence. There is no missing precomputed weight file to obtain. Patients are grouped across bootstrap copies; training-fold standardization, equal/signed gene means, ridge penalty 1, five-fold repeated validation and patient-averaged held-out probabilities match the original S11 procedure. It outputs aggregate metrics only. All 24 point estimates and interval endpoints were reproduced from the rebuilt sources to less than 1.2 × 10^-16.

The complete public reference entry point was also rerun from the reconstructed sources: all six core aggregate CSV files and all 99 archived NumPy arrays were reproduced exactly.

The reference runner uses the four-cohort intersection of uniquely measured, finite, nonconstant approved protein-coding genes for the broad background. The historical and inflamed backgrounds are reconstructed from all five baseline expression matrices without response labels. It draws 10,000 four-gene sets per background. The rank sensitivity uses a fixed common gene membership and within-sample centered percentile ranks over the broad universe, with 5,000 paired outcome-stratified bootstrap draws. Original published assay algorithms are not approximated by these gene-content comparisons.

## Separate PUCH check

`external_PUCH/PUCH_README.txt` documents the pinned author sources, sample reconciliation, expression scale and endpoint. The fixed four-gene and IFN-4 scores were correlated (rho 0.848); their AUCs were 0.779 and 0.744. The paired difference was 0.035 (95% interval −0.075 to 0.148). Stable disease was retained as nonresponse, matching this manuscript rather than the original source paper's SD-excluding ROC analysis. This comparison does not test a fitted incremental model.

```sh
python external_PUCH/puch_external_check.py --download --input-dir ../v19_PUCH_raw --output-dir ../v19_PUCH_results
```

## Definitions and analysis chronology

The original gene-pool definitions and numerical draws can be traced without patient data through `code/run_reference_public.py`, `results/reference/run_info.json`, `results/reference/gene_eligibility_audit.csv.gz`, the three `universe_*.csv` files and their `draws_*.csv.gz` files. The 15,561-gene pool is the four-cohort intersection of uniquely measured, finite, nonconstant HGNC Approved protein-coding genes, using the included frozen HGNC snapshot. It is not the complete human transcriptome. The 733-gene pool reproduces genes shared across the five original source-symbol matrices; the 142-gene subset requires Spearman rho >= 0.5 with CD8A in every one of those five baseline populations. Outcome labels do not define these pools. Reference-distribution 2.5th and 97.5th percentiles describe variation between random gene sets, not patient-bootstrap confidence intervals.

The historical execution specifications are in the cited v16 archive under the following paths relative to `Melanoma_immune_benchmark_v16_public/`. These are dated retrospective records, not prospective registrations.

| Record | Archived path and scope |
|---|---|
| Core v7, 2026-09-07 | `modules/Melanoma_immune_panel_Reproducibility_Data_v7/protocols/analysis_lock_v7.md`: fixed scoring and response definitions; states that Yale/Riaz outcomes and Hugo eligibility metadata had already been inspected. An earlier manuscript also reported outcome analyses for Hugo. |
| Addendum 07, 2026-11-02 | `modules/Melanoma_immune_panel_Reproducibility_Addendum07/protocols/analysis_lock_addendum_07.md`: IFN-4, pooling, equivalence and the 733/142 backgrounds; explicitly records earlier same-day exploratory pilots. |
| Addendum 08, 2026-11-02 | `modules/Melanoma_immune_panel_Reproducibility_Addendum08/protocols/analysis_lock_addendum_08.md`: single-cell response and bulk marker analyses. |
| Addendum 10, 2026-11-05 | `modules/Melanoma_immune_panel_Reproducibility_Addendum10/protocols/analysis_lock_addendum_10.md`: published gene-set adaptations, coverage and size-matched references. |
| Addendum 11, 2026-11-05 | `modules/Melanoma_immune_panel_Reproducibility_Addendum11/addendum/analysis_lock_addendum_11.md`: paired four-gene incremental models, full-pipeline bootstrap and signed GEP reference; explicitly post hoc. |

The broad background, new rank sensitivity and five additional incremental comparisons were revision-stage post hoc extensions on already available cohorts. Their definitions and implementation are documented in this supplement; they are not claimed as choices made before the original outcomes were known.

An additional source file included with v20 was `external_PUCH/PUCH_external_screen_spec.txt`, copied byte-for-byte from the original local PUCH eligibility/analysis record. It is identical to the later working copy. SHA-256: `697cc886cd3fe3aa1c5164a96ffa052ee9b3f1e6f52e9c72022ea5ace211a97f`. The original local file metadata records creation on 2026-11-23 at 20:32:05 Asia/Singapore, before the locally recorded clinical-file acquisition (20:32:38) and final aggregate results (20:53:22). These local records support the execution sequence; they do not constitute independent time-stamped registration. No date or other text was inserted into the original specification. The v19-labelled internal figures and numerical records retain their original names and content.

## Interpretation and reuse

The added analyses remain exploratory, and their intervals are not multiplicity adjusted. Internally validated increments are conditional on the specified models, genes and observed cohort populations. Wide intervals do not establish no increment. Random-reference percentiles are conditional on their gene pool; the inflamed pool was defined using expression correlation in the evaluated cohorts. The rank sensitivity tests a different scoring rule and does not demonstrate a transportable clinical assay.

Author-created code and aggregate results retain the project's MIT and CC BY 4.0 licenses, respectively (see `LICENSES/`). HGNC nomenclature data are CC0; their source is recorded in `resources/mapping_manifest.json`. These licenses do not grant new redistribution rights for separately acquired third-party clinical or expression data. New local preparations and any optional private audit outputs should remain outside the distributed supplement.
