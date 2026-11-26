#!/usr/bin/env python3
"""Post hoc audit of Ayers refined expanded immune (18 genes, Table 2) and Davoli (7 genes) under the common scoring rule."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PANELS = {
    "Ayers_expanded_18": (
        "CD3D IDO1 CIITA CD3E CCL5 GZMK CD2 HLA-DRA CXCL13 IL2RG NKG7 "
        "HLA-E CXCR6 LAG3 TAGAP CXCL10 STAT1 GZMB"
    ).split(),
    "Davoli_7": "CD247 CD2 CD3E GZMH NKG7 PRF1 GZMK".split(),
}


def score_common(matrix, genes, negative_gene=None):
    values = matrix.loc[genes]
    standardized = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1, ddof=1), axis=0)
    if negative_gene:
        standardized.loc[negative_gene] *= -1
    return standardized.mean(axis=0)


def calc_auc(values, response):
    pos = values.loc[response.eq(1)].to_numpy()
    neg = values.loc[response.eq(0)].to_numpy()
    diff = pos[:, None] - neg
    return float(((diff > 0) + 0.5 * (diff == 0)).mean())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="Additional file 3 root")
    parser.add_argument("--workspace", type=Path, required=True, help="Prepared four-cohort workspace")
    parser.add_argument("--output", type=Path, required=True, help="Output directory for CSV summaries")
    args = parser.parse_args()

    archive = args.archive.resolve()
    sys.path[:0] = [str(archive / s) for s in ("code", "tcga_panels", "external_TCGA_background")]
    from run_increment_public import PATHS
    from run_panel_reference import simulation, summarize
    import run_tcga_background as background

    reference = archive / "results/reference"
    eligible = set(pd.read_csv(reference / "universe_broad_protein_coding.csv").gene)
    shared = {name: [g for g in genes if g in eligible] for name, genes in PANELS.items()}
    coverage = pd.read_csv(reference / "panel_coverage.csv")
    published = pd.read_csv(reference / "scoring_auc.csv")

    hgnc = pd.read_csv(archive / "resources/hgnc_complete_set.txt", sep="\t", dtype=str)
    hgnc = hgnc[hgnc.status.eq("Approved") & hgnc.entrez_id.notna()]
    unique = ~hgnc.symbol.duplicated(keep=False) & ~hgnc.entrez_id.duplicated(keep=False)
    symbols = dict(zip(hgnc.loc[unique, "entrez_id"].astype(int), hgnc.loc[unique, "symbol"]))

    cohort_rows = []
    for cohort, rel_path in PATHS.items():
        matrix = pd.read_csv(args.workspace / rel_path, index_col=0)
        clinical = pd.read_csv(args.workspace / "prepared_v7/inputs/bulk" / cohort / "clinical.csv", index_col=0)
        matrix = matrix[clinical.index]
        if cohort == "GSE91061":
            matrix.index = matrix.index.map(symbols)
            matrix = matrix[matrix.index.notna()]
        assert matrix.index.is_unique and clinical.index.is_unique

        std_genes = coverage[coverage.cohort.eq(cohort)].set_index("panel").fixed_common_genes
        ifn4 = score_common(matrix, std_genes["IFN-4"].split(";"))
        gep = score_common(matrix, std_genes["GEP-18"].split(";"), negative_gene="CD276")
        rep = published.query('cohort == @cohort and panel == "GEP-18" and method == "z_common"')
        np.testing.assert_allclose(calc_auc(gep, clinical.response), rep.iloc[0].auc, rtol=0, atol=1e-12)

        for name, genes in shared.items():
            vals = score_common(matrix, genes)
            cohort_rows.append({
                "cohort": cohort,
                "panel": name,
                "baseline_n": len(clinical),
                "evaluable_n": int(clinical.response.notna().sum()),
                "common_genes_n": len(genes),
                "spearman_rho_IFN4": float(vals.rank().corr(ifn4.rank())),
                "standalone_auc": calc_auc(vals, clinical.response),
            })

    cohort_df = pd.DataFrame(cohort_rows)
    eligibility = pd.read_csv(archive / "external_TCGA_background/results/tcga_gene_eligibility.csv")
    pools = {
        "original142": sorted(eligibility.loc[eligibility.in_historical142, "gene"]),
        "tcga242": sorted(eligibility.loc[eligibility.eligible_external, "gene"]),
    }
    matrices, outcomes = background.read_evaluation_inputs(args.workspace / "prepared_v7", sorted(eligibility.gene))

    summary_rows = []
    with np.load(archive / "tcga_panels/results/reference_draws.npz") as draws:
        for name, genes in shared.items():
            sub = cohort_df[cohort_df.panel.eq(name)]
            row = {
                "panel": name,
                "source_genes_n": len(PANELS[name]),
                "common_genes_n": len(genes),
                "omitted_genes": ";".join(sorted(set(PANELS[name]) - set(genes))),
                "primary_median_rho_IFN4": float(sub.spearman_rho_IFN4.median()),
                "primary_mean_auc": float(sub.standalone_auc.mean()),
            }
            for pool_name, pool_genes in pools.items():
                key = f"{pool_name}_n{len(genes)}_all_positive_mean"
                if key in draws:
                    dist = draws[key]
                else:
                    dist = simulation(pool_genes, len(genes), matrices, outcomes, background)["all_positive_mean"]
                res = summarize(dist, row["primary_mean_auc"])
                row[f"{pool_name}_percentile"] = float(res["percentile"])
                row[f"{pool_name}_nominal_p_upper"] = float(res["p_upper"])
            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    args.output.mkdir(parents=True, exist_ok=True)
    cohort_df.to_csv(args.output / "posthoc_additional_panels_cohorts.csv", index=False)
    summary_df.to_csv(args.output / "posthoc_additional_panels.csv", index=False)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
