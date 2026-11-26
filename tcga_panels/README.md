# Published gene sets against two inflammatory gene backgrounds

We asked whether the relative standing of five published gene sets changes when a separate melanoma cohort defines the inflammatory gene pool. This exploratory extension was specified during revision, after the earlier results were known. The decision record is in `analysis_decisions.txt`.

The 142-gene pool required Spearman correlation with CD8A of at least 0.50 in each of the five historical baseline populations. The 242-gene pool applied the same threshold in 367 metastatic TCGA-SKCM specimens, where 731 of 733 historical candidates could be evaluated. The gene-eligibility table and code for rebuilding the TCGA screen are in `../external_TCGA_background/`. Both comparisons evaluate response in the same four treated melanoma cohorts. TCGA supplied the expression data for selecting reference genes.

## Results

| Gene set | Common genes | Mean AUC | Original pool percentile | TCGA pool percentile | TCGA draws at or above observed / 10,000 | TCGA upper-tail P |
|---|---:|---:|---:|---:|---:|---:|
| ICR | 19 | 0.621 | 97.19 | 98.37 | 163 | 0.0164 |
| GEP, one negative gene per random set | 17 | 0.616 | 83.28 | 91.05 | 895 | 0.0896 |
| Chemokine | 8 | 0.602 | 27.78 | 55.39 | 4,461 | 0.4462 |
| TLS | 8 | 0.590 | 4.14 | 23.96 | 7,604 | 0.7604 |
| APM | 6 | 0.648 | 100.00 | 99.96 | 4 | 0.00050 |
| GEP sensitivity, all-positive random sets | 17 | 0.616 | 85.80 | 93.63 | 637 | 0.0638 |

The observed scores use the published gene content with this study's signed, equally weighted gene means. For both GEP comparisons, CD276 enters the observed score negatively. The sensitivity changes only the signs assigned to the random lists. The APM comparison uses B2M, HLA-A, HLA-B, HLA-C, TAP1 and TAPBP; TAP2 is unavailable in Gide. The four cohort AUCs have equal weight. The percentiles describe discrimination relative to random inflammatory lists of the same size. Additional information conditional on IFN-4 is assessed separately by the incremental modeling analyses.

## Calculation

Each gene was standardized across all pretreatment patients in its cohort using the sample standard deviation. Response was evaluated among 49/51 Riaz, 26/26 Hugo, 102/104 Liu and 72/72 Gide patients. A single sampled list was scored in all four cohorts. Each list was sampled uniformly without replacement from its pool. Chemokine and TLS share the same eight-gene random lists. Within each pool and size, 10,000 draws used NumPy PCG64 with seed `20261500 + size`. The sign-matched GEP draws used an independent seed sequence `[20261105, 11, 17, 1]` to choose one negative gene per list.

The percentile equals 100 times the proportion of draws below the observed mean AUC, with half weight for ties. The nominal upper-tail probability equals (1 + draws at or above observed) / 10,001. The 2.5th and 97.5th percentiles in the results file describe variation across random lists. Upper-tail probabilities are unadjusted for multiple comparisons.

## Reproduction

From the Additional file 3 root, prepare the RNA-seq workspace following the main README. The package's archived gene-eligibility table supplies both pools. To regenerate that table from the TCGA source matrices, follow `external_TCGA_background/README.md`.

```sh
python tcga_panels/run_panel_reference.py --archive . --workspace ../v19_work --output ../panel_results
Rscript tcga_panels/verify_draws.R ../v19_work/prepared_v7 ../panel_results
```

The Python run checks all 20 observed cohort AUCs, reproduces the six original pool percentile/tail comparisons, and replays all 10,000 archived external four-gene draws in every response cohort. Every 250th draw is also checked using explicit responder/nonresponder pair comparisons. The separate base R script reads the prepared source matrices, applies sample-standard-deviation scaling and computes AUC by direct pair counting for 40 lists per cohort, including both GEP reference sign patterns.

## Output files

- `results/panel_reference_comparison.csv` contains the outcomes of all 12 pool, panel and sign comparisons, including the number of ties.
- `results/common_panel_observed_auc.csv` contains the 20 observed cohort-by-panel AUCs and cohort-level sample counts.
- `results/reference_draws.npz` contains the two sorted gene lists, zero-based sampled gene indices, negative positions for sign-matched GEP, cohort AUCs and per-draw arithmetic means. Open with `numpy.load(path, allow_pickle=False)`.
- `results/spot_check_sets.csv` records the gene symbols and cohort AUCs for draws 1, 2,501, 7,501 and 10,000 across the ten pool/size/sign combinations.
- `results/verification.json` and `results/independent_R_checks.csv` record the replay checks and arithmetic differences.
