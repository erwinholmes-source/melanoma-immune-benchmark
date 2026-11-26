"""Compare five published melanoma gene sets with two CD8A-correlated reference pools."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scipy

DRAW_COUNT = 10_000
PANEL_NAMES = ("ICR-20", "GEP-18", "Chemokine-12", "TLS-9", "APM-7")
POOL_NAMES = ("original142", "tcga242")


def direct_auc(scores, response):
    """Use all responder/nonresponder pairs, counting a tie as one half."""
    difference = scores[:, response, None] - scores[:, None, ~response]
    return ((difference > 0).mean(axis=(1, 2)) +
            0.5 * (difference == 0).mean(axis=(1, 2)))


def read_observations(workspace, archive, core, outcomes):
    # Use the same original baseline matrices and Riaz mapping as the public preparation.
    sys.path.insert(0, str(archive / "code"))
    from run_increment_public import PATHS
    hgnc = pd.read_csv(archive / "resources/hgnc_complete_set.txt", sep="\t", dtype=str)
    hgnc = hgnc.loc[hgnc.status.eq("Approved") & hgnc.entrez_id.notna()]
    valid = ~hgnc.symbol.duplicated(keep=False) & ~hgnc.entrez_id.duplicated(keep=False)
    entrez_to_symbol = dict(zip(hgnc.loc[valid, "entrez_id"].astype(int), hgnc.loc[valid, "symbol"]))
    definitions = pd.read_csv(archive / "results/reference/panel_coverage.csv", keep_default_na=False)
    target = pd.read_csv(archive / "results/reference/scoring_auc.csv").query("method == 'z_common'")
    output = []
    for cohort in core.COHORTS:
        clinical = pd.read_csv(workspace / "prepared_v7/inputs/bulk" / cohort / "clinical.csv", index_col=0)
        expr = pd.read_csv(workspace / PATHS[cohort], index_col=0)
        if cohort == "GSE91061":
            expr.index = expr.index.map(entrez_to_symbol)
            expr = expr.loc[expr.index.notna()]
        assert expr.index.is_unique and clinical.index.is_unique
        eligible = definitions[definitions.panel.isin(PANEL_NAMES) & definitions.cohort.eq(cohort)]
        genes_needed = sorted({g for gl in eligible.fixed_common_genes for g in gl.split(";")})
        expr = expr.loc[genes_needed, clinical.index]
        without_outcome = clinical.response.isna()
        z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=1), axis=0)
        y = clinical.loc[~without_outcome, "response"].to_numpy(dtype=bool)
        np.testing.assert_array_equal(y, outcomes[cohort])
        for panel in PANEL_NAMES:
            info = definitions.query("panel == @panel and cohort == @cohort").iloc[0]
            genes = info.fixed_common_genes.split(";")
            negative = set(filter(None, info.negative_genes.split(";")))
            signs = np.array([-1 if gene in negative else 1 for gene in genes])
            score = (z.loc[genes, ~without_outcome].to_numpy() * signs[:, None]).mean(axis=0)[None, :]
            auc = direct_auc(score, y)[0]
            old = target.query("cohort == @cohort and panel == @panel").iloc[0]
            np.testing.assert_allclose(auc, old.auc, rtol=0, atol=1e-14)
            output.append(dict(panel=panel, cohort=cohort, genes=len(genes),
                               observed_auc=auc, patients=len(y), responders=y.sum()))
    return pd.DataFrame(output)


def simulation(pool, count, matrices, outcomes, core):
    seed = 20261500 + count
    rng = np.random.default_rng(seed)
    indices = np.array([rng.choice(len(pool), count, replace=False) for _ in range(DRAW_COUNT)])
    signs = ["all_positive"]
    if count == 17:
        neg_rng = np.random.default_rng(np.random.SeedSequence([20261105, 11, 17, 1]))
        negative_position = neg_rng.integers(0, count, size=DRAW_COUNT)
        signs.append("one_negative")
    stored = {"gene_index": indices.astype(np.int16)}
    if count == 17:
        stored["negative_column"] = negative_position.astype(np.int8)
    for mode in signs:
        cohort_values = []
        for cohort, y in outcomes.items():
            z = matrices[cohort].loc[pool].to_numpy()
            block_values = []
            for start in range(0, DRAW_COUNT, 250):
                selected = z[indices[start:start+250]]
                score = selected.mean(axis=1)
                if mode == "one_negative":
                    positions = negative_position[start:start+len(score)]
                    score -= 2 * selected[np.arange(len(score)), positions] / count
                values = core.rank_auc(score, y)
                # Independent direct pair counts in the first draw of each block.
                np.testing.assert_allclose(values[0], direct_auc(score[:1], y)[0],
                                           rtol=0, atol=1e-14)
                block_values.append(values)
            cohort_values.append(np.concatenate(block_values))
        for cohort, values in zip(outcomes, cohort_values):
            stored[f"{mode}_{cohort}"] = values
        stored[f"{mode}_mean"] = np.mean(cohort_values, axis=0)
    return stored


def summarize(values, observed):
    lower, median, upper = np.quantile(values, [0.025, 0.50, 0.975])
    ties = np.count_nonzero(values == observed)
    exceed = np.count_nonzero(values >= observed)
    return dict(random_median=median, random_q025=lower, random_q975=upper,
                percentile=100 * ((values < observed).mean() + 0.5 * (values == observed).mean()),
                p_upper=(1 + exceed) / (1 + DRAW_COUNT), draws_ge_observed=exceed, ties=ties)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--archive", type=Path, required=True, help="Additional file 3 root directory")
    cli.add_argument("--workspace", type=Path, required=True, help="Output of code/prepare_sources.py")
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    archive = args.archive.resolve()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    ext_dir = archive / "external_TCGA_background"
    sys.path.insert(0, str(ext_dir))
    import run_tcga_background as core
    eligibility = pd.read_csv(ext_dir / "results/tcga_gene_eligibility.csv")
    pools = {"original142": sorted(eligibility.loc[eligibility.in_historical142, "gene"]),
             "tcga242": sorted(eligibility.loc[eligibility.eligible_external, "gene"])}
    assert len(pools["original142"]) == 142 and len(pools["tcga242"]) == 242
    matrices, outcomes = core.read_evaluation_inputs(args.workspace / "prepared_v7",
                                                     sorted(eligibility.gene))
    # Check the pooled four-gene comparator against every archived external draw.
    four = pd.read_csv(ext_dir / "results/draws_auc_tcga_cd8a.csv.gz")
    ix = pd.Index(pools["tcga242"]).get_indexer(four[[f"gene{j}" for j in range(1, 5)]]
                                               .to_numpy().ravel()).reshape(-1, 4)
    assert ix.min() >= 0 and len(ix) == DRAW_COUNT
    replay = core.compare_draws(pools["tcga242"], ix, matrices, outcomes)
    for c, values in replay.items():
        np.testing.assert_allclose(four[c], values, rtol=0, atol=1e-14)

    observations = read_observations(args.workspace, archive, core, outcomes)
    observations.to_csv(out / "common_panel_observed_auc.csv", index=False)
    internal = pd.read_csv(archive / "primary_panels/primary_panel_table.csv").set_index("panel")
    results = []
    spot_checks = []
    arrays = {}
    for pool_name, pool in pools.items():
        arrays[pool_name + "_symbols"] = np.asarray(pool)
        for size in sorted(set(observations.genes)):
            simulated = simulation(pool, size, matrices, outcomes, core)
            for sign in (["all_positive", "one_negative"] if size == 17 else ["all_positive"]):
                for idx in (0, 2500, 7500, 9999):
                    symbols = np.asarray(pool)[simulated["gene_index"][idx]]
                    negative = symbols[simulated["negative_column"][idx]] if sign == "one_negative" else ""
                    values = {c: simulated[f"{sign}_{c}"][idx] for c in outcomes}
                    spot_checks.append(dict(pool=pool_name, size=size, reference_signs=sign,
                                            draw=idx+1, genes=";".join(symbols),
                                            negative_gene=negative, **values))
            for key, arr in simulated.items():
                arrays[f"{pool_name}_n{size}_{key}"] = arr
            for panel, rows in observations.query("genes == @size").groupby("panel", sort=False):
                observed = np.mean(rows.observed_auc.to_numpy())
                for sign in (["all_positive", "one_negative"] if panel == "GEP-18" else ["all_positive"]):
                    summary = summarize(simulated[sign + "_mean"], observed)
                    if pool_name == "original142":
                        orig = internal.loc[panel]
                        expected = ((orig.gep_sign_matched_percentile, orig.gep_sign_matched_upper_p)
                                    if sign == "one_negative" else
                                    (orig.percentile_in_size_matched_random_sets, orig.nominal_upper_p))
                        np.testing.assert_allclose([summary["percentile"], summary["p_upper"]],
                                                   expected, rtol=0, atol=1e-12)
                    results.append(dict(panel=panel, pool=pool_name, reference_signs=sign,
                                        panel_genes=size, observed_mean_auc=observed,
                                        **summary, n_draws=DRAW_COUNT))
    result = pd.DataFrame(results)
    panel_order = {panel: i for i, panel in enumerate(PANEL_NAMES)}
    result["panel_order"] = result.panel.map(panel_order)
    result = result.sort_values(["pool", "panel_order", "reference_signs"]).drop(columns="panel_order")
    result.to_csv(out / "panel_reference_comparison.csv", index=False)
    np.savez_compressed(out / "reference_draws.npz", **arrays)
    pd.DataFrame(spot_checks).to_csv(out / "spot_check_sets.csv", index=False)
    check = {
        "observed_cohort_panel_AUCs_checked": len(observations),
        "internal_percentiles_and_tail_probabilities_reproduced": 6,
        "external_four_gene_draws_reproduced_in_each_primary_cohort": DRAW_COUNT,
        "direct_pair_aucs_per_size_sign_pool_and_cohort": DRAW_COUNT // 250,
        "seed_by_size": "20261500 + size",
        "negative_position_seed_sequence": [20261105, 11, 17, 1],
        "pool_audit_sha256": hashlib.sha256((ext_dir / "results/tcga_gene_eligibility.csv").read_bytes()).hexdigest(),
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
    }
    (out / "verification.json").write_text(json.dumps(check, indent=2) + "\n")
    print(result[["panel", "pool", "reference_signs", "observed_mean_auc", "percentile", "p_upper"]].to_string(index=False))


if __name__ == "__main__":
    main()
