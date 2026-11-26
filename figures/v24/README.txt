Figures 1, 2 and 5

From an extracted Additional file 3 root, run

python code/plot_figures_v24.py --archive . --output ../figure_replay --tables ../table_replay

--panel-results may specify another location for panel_reference_comparison.csv if the external five-panel reanalysis is installed outside tcga_panels/results/. Dependencies are given in requirements_figures.txt.

Figure 2 reads scoring_auc.csv, scoring_correlations.csv and panel_coverage.csv in results/reference/, primary_panels/primary_panel_table.csv and external_PUCH/extension_v21/results/PUCH_panel_results.csv. It compares exactly the same constituents in all five cohorts. The runner verifies all gene lists against the external release and reproduces the median primary-cohort correlations and equal-weight mean AUC for each published score. The companion fig2_values.csv includes the original 95% AUC intervals.

Figure 5 checks each distribution of 10,000 four-gene draws against the archived mean-cohort AUC, reference median and midrank percentile. The original cohorts' distributions are from results/reference/random_reference_auc.npz. The TCGA-defined four-gene reference is from external_TCGA_background/results/. The two size-matched reference results for each published set are read from tcga_panels/results/panel_reference_comparison.csv. The GEP plot uses sign-matched draws from both pools. The companion fig5_values.csv contains all figure estimates.

Figures are vector PDF and 450-dpi PNG. The script relies on aggregate files and does not require access to participant-level data.

Figure 6 comes from results/increment/increment_all_results.csv and external_PUCH/extension_v21/results/PUCH_frozen_model_performance.csv. Run python code/plot_increment.py --archive . --output ../figure_replay.
