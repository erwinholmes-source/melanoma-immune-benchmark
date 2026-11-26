"""Assess response AUCs for random sets selected using external TCGA-SKCM coexpression."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

COHORTS = {
    'GSE91061': ('Sample', 51, 49, 10),
    'GSE78220': ('sample_id', 26, 26, 14),
    'Liu2019': ('patient_id', 104, 102, 41),
    'Gide2019': ('patient_id', 72, 72, 39),
}
PANEL = ['CD7', 'CXCL9', 'LYZ', 'STAT1']
TCGA_HASHES = {
    'data_mrna_seq_v2_rsem.txt': 'de12a5483d43195976e022b5c81c41ed4f5daa6fe4a09613e7db3f0d7c7bec4c',
    'data_clinical_sample.txt': '1da0e866ba04db38a3efc98da78bff1ffb6204558859c7546683f791995fc90d',
}
DRAWS = 10000
SEED = 20261124


def external_screen(tcga, archive, candidates):
    for name, expected in TCGA_HASHES.items():
        actual = hashlib.sha256((tcga / name).read_bytes()).hexdigest()
        if expected != actual:
            raise ValueError(f'Wrong TCGA source file checksum for {name}: {actual}')
    samples = pd.read_csv(tcga / 'data_clinical_sample.txt', sep='\t', comment='#',
                          usecols=['SAMPLE_ID', 'SAMPLE_TYPE', 'PATIENT_ID'])
    source = pd.read_csv(tcga / 'data_mrna_seq_v2_rsem.txt', sep='\t', low_memory=False)
    metastases = samples[samples.SAMPLE_TYPE.eq('Metastasis') & samples.SAMPLE_ID.isin(source.columns)]
    assert metastases.PATIENT_ID.is_unique and metastases.SAMPLE_ID.is_unique and len(metastases) == 367
    mapping = pd.read_csv(archive / 'resources/fixed_rank_universe.csv').set_index('symbol')
    mapping = mapping.loc[candidates]
    assert mapping.index.is_unique and mapping.entrez_id.is_unique and mapping.entrez_id.notna().all()
    counts = source.Entrez_Gene_Id.value_counts()
    by_id = source.loc[~source.Entrez_Gene_Id.duplicated(keep=False)].set_index('Entrez_Gene_Id')
    matrix = by_id.reindex(mapping.entrez_id)[metastases.SAMPLE_ID].set_axis(candidates, axis=0)
    allfinite = np.isfinite(matrix).all(axis=1)
    varies = matrix.max(axis=1) > matrix.min(axis=1)
    valid = allfinite & varies
    assert valid['CD8A']
    ranks = matrix.loc[valid].rank(axis=1, method='average')
    correlation = ranks.corrwith(ranks.loc['CD8A'], axis=1)
    audit = mapping.copy()
    audit.index.name = 'gene'
    audit['n_tcga_rows_for_entrez'] = mapping.entrez_id.map(counts).fillna(0).astype(int)
    audit['tcga_source_symbol'] = mapping.entrez_id.map(by_id.Hugo_Symbol)
    audit['finite_all_metastases'] = allfinite
    audit['variable'] = varies
    audit['rho_with_CD8A'] = correlation
    audit['eligible_external'] = correlation.reindex(candidates).ge(0.5)
    audit['index_gene'] = audit.index.isin(PANEL)
    return audit, len(metastases)


def read_evaluation_inputs(prepared, candidates):
    matrices, outcomes = {}, {}
    for cohort, (id_col, n_pre, n_eval, n_resp) in COHORTS.items():
        clinical = pd.read_csv(prepared / 'inputs/bulk' / cohort / 'clinical.csv').set_index(id_col)
        expr = pd.read_csv(prepared / 'inputs/rank' / cohort / 'source_symbol_expression.csv.gz', index_col=0)
        expr = expr.loc[candidates, clinical.index]
        assert clinical.index.is_unique and expr.index.is_unique and np.isfinite(expr).all().all()
        z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=1), axis=0)
        y = clinical.response.dropna()
        assert (len(clinical), len(y), y.sum()) == (n_pre, n_eval, n_resp)
        assert y.isin([0, 1]).all() and np.isfinite(z).all().all()
        matrices[cohort] = z[y.index]
        outcomes[cohort] = y.to_numpy(dtype=bool)
    return matrices, outcomes


def rank_auc(scores, responders):
    n1 = int(responders.sum())
    n0 = len(responders) - n1
    r = rankdata(scores, axis=1, method='average')
    return (r[:, responders].sum(axis=1) - n1 * (n1 + 1) / 2) / (n1 * n0)


def compare_draws(genes, draws, matrices, outcomes):
    vals = {}
    for cohort in COHORTS:
        z = matrices[cohort].loc[genes].to_numpy()
        scores = z[draws].mean(axis=1)
        y = outcomes[cohort]
        vals[cohort] = rank_auc(scores, y)
        # Independently count all response/nonresponse pairs on the first 20 draws.
        checked = scores[:20]
        diff = checked[:, y, None] - checked[:, None, ~y]
        direct = (np.sum(diff > 0, axis=(1, 2)) + 0.5 * np.sum(diff == 0, axis=(1, 2))) / (y.sum() * (~y).sum())
        np.testing.assert_allclose(vals[cohort][:len(checked)], direct, rtol=0, atol=1e-14)
    vals['mean_independent4'] = np.mean(list(vals.values()), axis=0)
    return vals


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--archive', type=Path, required=True, help='Additional file 3 package root')
    cli.add_argument('--prepared-v7', type=Path, required=True, help='prepared_v7 from prepare_sources.py')
    cli.add_argument('--tcga', type=Path, required=True, help='Folder with the two original, hash-checked PanCancer Atlas source files')
    cli.add_argument('--output', type=Path, required=True)
    args = cli.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    refs = args.archive / 'results/reference'
    candidates = sorted(pd.read_csv(refs / 'universe_historical733.csv').gene)
    historical = sorted(pd.read_csv(refs / 'universe_inflamed142.csv').gene)
    assert len(candidates) == 733 and len(historical) == 142 and set(PANEL).issubset(historical)
    audit, n = external_screen(args.tcga, args.archive, candidates)
    audit['in_historical142'] = audit.index.isin(historical)
    audit.to_csv(out / 'tcga_gene_eligibility.csv')
    external = sorted(audit.index[audit.eligible_external])
    matrices, outcomes = read_evaluation_inputs(args.prepared_v7, candidates)
    panel = {c: rank_auc(matrices[c].loc[PANEL].mean(axis=0).to_numpy()[None, :], outcomes[c])[0] for c in COHORTS}
    panel['mean_independent4'] = np.mean(list(panel.values()))

    # Replay every archived random draw before adding external/exclusion backgrounds.
    original_draws = pd.read_csv(refs / 'draws_inflamed142.csv.gz')
    locations = pd.Index(historical).get_indexer(original_draws.iloc[:, 1:].to_numpy().ravel()).reshape(-1, 4)
    assert np.all(locations >= 0) and locations.shape == (DRAWS, 4)
    old_auc = compare_draws(historical, locations, matrices, outcomes)
    target = pd.read_csv(refs / 'null_auc_inflamed142.csv.gz')
    for c, values in old_auc.items():
        np.testing.assert_allclose(target[c], values, rtol=0, atol=1e-14)
    observed_target = pd.read_csv(refs / 'random_reference_summary.csv').query("background == 'inflamed142'").set_index('cohort')
    for c, obs in panel.items():
        np.testing.assert_allclose(obs, observed_target.loc[c, 'panel_auc'], rtol=0, atol=1e-14)

    pools = {'tcga_cd8a': external, 'tcga_cd8a_without_panel': sorted(set(external) - set(PANEL)),
             'historical142_without_panel': sorted(set(historical) - set(PANEL))}
    summaries = []
    gzip = {'method': 'gzip', 'mtime': 0}
    for name, genes in pools.items():
        rng = np.random.default_rng(SEED)
        draws = np.stack([rng.choice(len(genes), size=4, replace=False) for _ in range(DRAWS)])
        vals = compare_draws(genes, draws, matrices, outcomes)
        draws_frame = pd.DataFrame(np.asarray(genes)[draws], columns=['gene1', 'gene2', 'gene3', 'gene4'])
        draws_frame.insert(0, 'draw', np.arange(1, DRAWS + 1))
        for c, values in vals.items():
            obs = panel[c]
            q025, median, q975 = np.quantile(values, [0.025, 0.5, 0.975])
            pct = 100 * ((values < obs).mean() + 0.5 * (values == obs).mean())
            p = (1 + np.sum(values >= obs)) / (1 + DRAWS)
            summaries.append(dict(background=name, cohort=c, genes=len(genes), draws=DRAWS, seed=SEED,
                                  panel_auc=obs, random_median=median, random_q025=q025, random_q975=q975,
                                  panel_percentile=pct, p_upper_nominal=p))
            draws_frame[c] = values
        draws_frame.to_csv(out / f'draws_auc_{name}.csv.gz', index=False, compression=gzip)
    pd.DataFrame(summaries).to_csv(out / 'tcga_random_reference_summary.csv', index=False)
    info = {
        'sample_selection': 'TCGA-SKCM RNA-profiled metastatic specimens, one per patient',
        'tcga_n': n, 'candidate_genes': len(candidates), 'tcga_genes_with_correlation': int(audit.rho_with_CD8A.notna().sum()),
        'genes_per_pool': {k: len(v) for k, v in pools.items()},
        'n_shared_between_original142_and_tcga': len(set(external) & set(historical)),
        'panel_genes_in_tcga_pool': sorted(set(external) & set(PANEL)),
        'tcga_source_sha256': TCGA_HASHES, 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__,
        'validation': {'original_observed_auc_reproduced_all_cohorts': True,
                       'all_10000_historical_inflamed_draws_reproduced_all_cohorts': True,
                       'first_20_draws_each_pool_directly_checked_in_all_cohorts': True},
    }
    (out / 'tcga_run_info.json').write_text(json.dumps(info, indent=2) + '\n')
    print(json.dumps(info, indent=2))
    print(pd.DataFrame(summaries).query("cohort == 'mean_independent4'").to_string(index=False))


if __name__ == '__main__':
    main()
