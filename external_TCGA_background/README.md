# External expression reference from metastatic TCGA-SKCM

This exploratory analysis asks how a random four-gene reference changes when CD8A coexpression is measured in external melanoma biopsies. The earlier study results were available when we specified this check. The original 733-gene candidate list controls assay coverage across all comparisons. The response cohorts supply expression and outcomes for AUC evaluation. TCGA supplies the gene correlations only. The external definition uses one cohort, whereas the original definition required the threshold in each of five cohorts.

## Inputs and reproduction

Use the package's existing requirements and follow the top-level README through source preparation to create the score-free `prepared_v7/` directory. Download the two pinned [cBioPortal TCGA-SKCM PanCancer Atlas](https://www.cbioportal.org/study/summary?id=skcm_tcga_pan_can_atlas_2018) files into a separate local directory. The script verifies each file's SHA-256 against the original source used in the analysis.

| File | Pinned source | SHA-256 |
| --- | --- | --- |
| data_mrna_seq_v2_rsem.txt | [RNA-seq RSEM estimates](https://media.githubusercontent.com/media/cBioPortal/datahub/86690e1ed9752b1dcd50b5657f5f05eafa4b6b78/public/skcm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt) | de12a5483d43195976e022b5c81c41ed4f5daa6fe4a09613e7db3f0d7c7bec4c |
| data_clinical_sample.txt | [Sample annotations](https://media.githubusercontent.com/media/cBioPortal/datahub/438e7e4169200607196a8833a07a076224bc6d9d/public/skcm_tcga_pan_can_atlas_2018/data_clinical_sample.txt) | 1da0e866ba04db38a3efc98da78bff1ffb6204558859c7546683f791995fc90d |

Run the command below from this analysis folder. Replace the paths with the actual locations. The archive root contains `resources/fixed_rank_universe.csv` and `results/reference/`.

```sh
python run_tcga_background.py --archive /path/to/additional-file-3 --prepared-v7 /path/to/prepared_v7 --tcga /path/to/TCGA_inputs --output /path/to/new_results
```

## Cohorts and scoring

All 367 metastatic specimens with profiled expression had unique patient identifiers. Among the 733 candidates, archived Entrez identifiers uniquely matched 731 finite, varying genes in TCGA. The identifier for FGF13 had two rows, and the identifier for MICA was absent. We used Spearman correlations of RSEM estimates across all 367 specimens, with mean ranks for ties. No survival or treatment outcome was used for gene selection. The CD8A correlation threshold of 0.50 identified 242 genes, including 141 of the original 142 and all four index genes. The correlation for the remaining original gene, TNFSF13, was 0.438.

For every draw, the same four distinct genes were assessed in Riaz, Hugo, Liu and Gide. We averaged the within-cohort gene z-scores for each individual after standardizing in the baseline population, then computed each AUC in the response-evaluable population. Higher expression indicated response throughout. We compared each draw's arithmetic mean of the four cohort AUCs with the analogous observed four-gene value of 0.621140. Each new pool had 10,000 draws with seed 20261124. Percentiles give ties half weight. The upper-tail empirical P uses a plus-one correction. These are exploratory reference comparisons.

| Random pool | Genes | Median random mean AUC | Central 95% of draws | Observed percentile | Empirical upper-tail P |
| --- | ---: | ---: | --- | ---: | ---: |
| Original CD8A-correlated genes | 142 | 0.607 | 0.576–0.630 | 86.95 | 0.131 |
| Original, four tested genes removed | 138 | 0.606 | 0.575–0.630 | 87.54 | 0.125 |
| TCGA-defined CD8A-correlated genes | 242 | 0.597 | 0.550–0.632 | 90.60 | 0.094 |
| TCGA-defined, four tested genes removed | 238 | 0.596 | 0.549–0.631 | 91.49 | 0.085 |

The observed four-gene score lay below the 95th percentile in each of these CD8A-correlated pools, compared with the 99.72nd percentile in the broad protein-coding background. The central 95% range describes variation between gene draws. The external comparison retains the original four cohorts for response evaluation.

## Outputs and numerical checks

- `results/tcga_gene_eligibility.csv` contains the Entrez match status, TCGA CD8A correlation and membership for each of the 733 candidates.
- `results/draws_auc_*.csv.gz` lists the genes and the four cohort AUCs for each draw. These files contain gene-level and cohort-level results.
- `results/tcga_random_reference_summary.csv` contains cohort-specific and equal-weight summaries. `results/tcga_run_info.json` records software versions, input hashes and check results.
- The script reproduces all four original observed AUCs and the archived cohort AUCs for every one of the 10,000 original inflamed-pool draws. It also checks AUCs for the first 20 draws in every pool against direct response/nonresponse pair counts.
