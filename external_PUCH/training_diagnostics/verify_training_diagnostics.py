#!/usr/bin/env python3
"""Independently replay the R training-only summaries using the archived Python feature pipeline."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def read_scores(package, workspace, frozen):
    sys.path[:0] = [str(package / 'code'), str(package / 'external_PUCH/extension_v21')]
    import puch_extension as transport
    import run_increment_public as original

    reference = package / 'results/reference/universe_historical733.csv'
    assert transport.sha256(reference) == frozen['rank_universe_sha256']
    universe = pd.read_csv(reference).gene.tolist()
    hgnc_file = package / 'resources/hgnc_complete_set.txt'
    assert transport.sha256(hgnc_file) == frozen['hgnc_sha256']
    hgnc = pd.read_csv(hgnc_file, sep='\t', dtype=str)
    approved = hgnc.loc[hgnc.status.eq('Approved') & hgnc.entrez_id.notna()]
    unique = approved.loc[~approved.symbol.duplicated(keep=False) & ~approved.entrez_id.duplicated(keep=False)]
    symbols = dict(zip(unique.entrez_id.astype(int), unique.symbol))
    frames = []
    for cohort in original.COHORTS:
        clinical_file = workspace / f'prepared_v7/inputs/bulk/{cohort}/clinical.csv'
        expression_file = workspace / original.PATHS[cohort]
        assert transport.sha256(clinical_file) == frozen['training_files'][cohort]['clinical_sha256']
        assert transport.sha256(expression_file) == frozen['training_files'][cohort]['expression_sha256']
        labels = pd.read_csv(clinical_file, index_col=0).response.dropna().astype(int)
        labels.index = labels.index.astype(str)
        expression = pd.read_csv(expression_file, index_col=0)
        expression.columns = expression.columns.astype(str)
        if cohort == 'GSE91061':
            expression.index = expression.index.astype(int)
            expression = expression.loc[expression.index.isin(symbols)]
            expression.index = expression.index.map(symbols)
        ranks = transport.fixed_rank_features(expression.loc[:, labels.index], universe, frozen['gene_sets'])
        ranks.columns = ['IFN4', 'APM6']
        ranks['cohort'], ranks['y'] = cohort, labels
        frames.append(ranks)
    data = pd.concat(frames, ignore_index=True)
    features = ['IFN4', 'APM6']
    np.testing.assert_allclose(data[features].mean(), frozen['predictor_mean'], rtol=0, atol=1e-12)
    np.testing.assert_allclose(data[features].std(), frozen['predictor_sample_sd'], rtol=0, atol=1e-12)
    data[features] = (data[features] - frozen['predictor_mean']) / frozen['predictor_sample_sd']
    return data, original.COHORTS


def check_coefficients(row, data, cohort_order):
    if hasattr(row, 'omitted'):
        data = data.loc[data.cohort.ne(row.omitted)]
    if row.intercept == 'common':
        intercepts = np.ones((len(data), 1))
        names = ['intercept_common']
    else:
        intercepts = np.column_stack([data.cohort.eq(c) for c in cohort_order])
        names = [f'intercept_{c}' for c in cohort_order]
    slopes = row.model.split('+')
    design = np.column_stack((intercepts, data[slopes]))
    names += ['beta_' + feature for feature in slopes]
    penalty = np.r_[np.zeros(intercepts.shape[1]), np.ones(len(slopes))]
    y = data.y.to_numpy()
    if row.weighting == 'equal_patients':
        w = np.full(len(data), 249/len(data))
    else:
        sizes = data.cohort.value_counts()
        w = 249/(len(sizes) * data.cohort.map(sizes).to_numpy())
    assert abs(w.sum() - row.total_weight) < 1e-10 and row.n == len(y) and row.responders == y.sum()
    b = np.zeros(design.shape[1])
    for unused_iteration in range(50):
        p = np.exp(-np.logaddexp(0, -design @ b))
        g = design.T @ (w*(p-y)) + penalty*b
        if np.max(np.abs(g)) < 1e-11:
            break
        hessian = design.T @ ((w*p*(1-p))[:, None] * design) + np.diag(penalty)
        b -= np.linalg.solve(hessian, g)
    else:
        raise ArithmeticError('Newton solver failed to converge')
    expected = np.array([getattr(row, name) for name in names])
    return np.max(np.abs(b - expected)), np.max(np.abs(g))


def check_summary(row, data):
    group = data if row.cohort == 'Pooled' else data.loc[data.cohort.eq(row.cohort)]
    pos = group.y.eq(1).to_numpy()
    def auc(name):
        scores = group[name].to_numpy()
        comparisons = scores[pos, None] - scores[~pos]
        return np.mean((comparisons > 0) + 0.5*(comparisons == 0))
    assert len(group) == row.n and pos.sum() == row.responders
    left, right = group.IFN4, group.APM6
    values = [auc('IFN4'), auc('APM6'), left.rank().corr(right.rank()), left.corr(right), left.mean(), right.mean()]
    expected = [row.auc_IFN4, row.auc_APM6, row.rho_IFN4_APM6, row.pearson_IFN4_APM6, row.mean_zIFN4, row.mean_zAPM6]
    return np.max(np.abs(np.array(values) - expected))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('package-root', 'workspace', 'results'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    package = args.package_root.resolve()
    frozen = json.loads((package/'external_PUCH/extension_v21/results/frozen_models.json').read_text())
    data, cohorts = read_scores(package, args.workspace, frozen)
    assert len(data) == frozen['training_patients'] and data.y.sum() == frozen['responders']
    summaries = pd.read_csv(args.results/'cohort_rank_summaries.csv')
    fits = pd.read_csv(args.results/'training_model_coefficients.csv')
    omissions = pd.read_csv(args.results/'leave_one_cohort_out.csv')
    assert len(summaries) == 5 and len(fits) == 12 and len(omissions) == 4
    observed_configs = set(fits[['intercept','weighting','model']].itertuples(index=False, name=None))
    expected_configs = {
        (i,w,m) for i in ('common','cohort') for w in ('equal_patients','equal_cohorts')
        for m in ('IFN4','APM6','IFN4+APM6')}
    assert observed_configs == expected_configs and set(omissions.omitted) == set(cohorts)
    rows = [*fits.itertuples(), *omissions.itertuples()]
    coef_errors, grad_errors = zip(*(check_coefficients(row, data, cohorts) for row in rows))
    summary_error = max(check_summary(row, data) for row in summaries.itertuples())
    assert max(coef_errors) < 1e-6 and summary_error < 1e-12 and max(grad_errors) < 1e-10
    report = {
        'status':'passed', 'n_training_patients':len(data), 'responders':int(data.y.sum()), 'external_rows_loaded':0,
        'input_checks':'Frozen rank universe, HGNC file and all eight training expression/clinical file hashes matched.',
        'scope':'Python average ranks and pairwise AUC, direct weighted Newton ridge fits versus base-R ranks, U-statistic and nlminb.',
        'summary_rows':len(summaries), 'coefficient_rows':len(fits), 'omission_rows':len(omissions),
        'max_abs_summary_error':float(summary_error), 'max_abs_coefficient_error':float(max(coef_errors)),
        'max_abs_python_gradient':float(max(grad_errors)), 'numpy':np.__version__, 'pandas':pd.__version__,
    }
    (args.results/'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
