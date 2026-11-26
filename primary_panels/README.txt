Primary-cohort comparison of five published gene sets

From this directory, run:
  python3 build_primary_panel_table.py

Python 3 standard library only. By default the runner reads the two aggregate
files scoring_auc.csv and scoring_correlations.csv from ../results/reference.
Use --reference-dir if the directory layout differs. The runner overwrites
primary_panel_table.csv or a file specified by --output.

Data sources
* scoring_auc.csv and scoring_correlations.csv use the fixed lists of genes
  shared across Riaz, Hugo, Liu and Gide (method z_common). Gene expression was
  standardized within each cohort. The GEP adaptation reverses the CD276 sign.
* a10_random_sets.csv is the unaltered aggregate output of addendum 10,
  comparing mean AUC across the four cohorts with 10,000 random gene sets
  of the same size drawn from the 142-gene CD8A-correlated pool or the wider
  733-gene pool. This table uses only the 142-gene pool. Random sets used
  all-positive signs and were sampled without replacement within each set.
* a11_gep_signed_null.csv is the unaltered aggregate sign-matched GEP
  sensitivity from addendum 11, with one gene given a negative sign per
  set. The two null-table CSV files are also included in the published source
  archive (Zenodo DOI 10.5281/zenodo.24001782, see the main manuscript
  reference to the historical code archive for version provenance).

Each median rho is the unweighted median of the four within-cohort Spearman
correlations with IFN-4 (253 pretreatment patients in all). Each mean AUC
gives equal weight to the four cohort AUCs (249 response-evaluable patients,
104 responders). Ranges are observed cohort ranges. The percentiles use the
same fixed gene lists and give ties half weight. The upper-tail probability
equals (1 + number of draws at or above the observed mean AUC)/10,001.
These nominal probabilities describe gene-set sampling. They do not quantify
uncertainty from sampling patients. No patient-level data are included here.
