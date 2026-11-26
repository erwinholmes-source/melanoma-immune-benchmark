#!/usr/bin/env python3
"""Post hoc PUCH ranking experiment under one fixed set of alternate training assumptions."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def verify_rank_auc(scores, pos, neg, positive_draws, negative_draws):
    # Rank-based U statistic, independently of the pairwise win matrix used for the main estimates.
    npos, nneg = len(pos), len(neg)
    offset = npos * (npos + 1) / 2
    points = pd.Series(scores).rank(method='average').to_numpy()
    auc = (points[pos].sum() - offset) / (npos * nneg)
    draws = np.column_stack((pos[positive_draws], neg[negative_draws]))
    ranked = pd.DataFrame(scores[draws]).rank(axis=1, method='average').to_numpy()
    boot = (ranked[:, :npos].sum(axis=1) - offset) / (npos * nneg)
    return auc, boot


def run(package, inputs, results):
    sys.path[:0] = [str(package / 'external_PUCH'), str(package / 'external_PUCH/extension_v21')]
    import puch_external_check as puch
    import puch_extension as transport

    frozen = json.loads((package/'external_PUCH/extension_v21/results/frozen_models.json').read_text())
    universe_path = package/'results/reference/universe_historical733.csv'
    assert transport.sha256(universe_path) == frozen['rank_universe_sha256']
    universe = pd.read_csv(universe_path).gene.tolist()
    assert len(universe) == len(set(universe)) == 733
    fits = pd.read_csv(results/'training_model_coefficients.csv')
    selected = fits.loc[fits.intercept.eq('cohort') & fits.weighting.eq('equal_patients')].set_index('model')
    assert selected.index.is_unique and len(selected) == 3
    slopes = {'IFN4': np.array([selected.at['IFN4','beta_IFN4']]),
              'IFN4+APM6': selected.loc['IFN4+APM6', ['beta_IFN4','beta_APM6']].to_numpy(dtype=float)}
    assert np.isfinite(np.concatenate(list(slopes.values()))).all()

    expr_path = inputs/'mel_puch_exp_data.csv'
    clin_path = inputs/'mel_puch_cli_data.csv'
    assert puch.digest(expr_path)['git_blob_sha1'] == puch.EXP_BLOB
    assert puch.digest(clin_path)['git_blob_sha1'] == puch.CLI_BLOB
    expr = pd.read_csv(expr_path, index_col=0)
    assert expr.shape == (20293,55) and expr.columns.is_unique and expr.index.is_unique
    features = transport.fixed_rank_features(expr, universe, frozen['gene_sets'])
    assert list(features.columns) == ['IFN4_rank','APM6_rank']
    standardized = (features.to_numpy() - frozen['predictor_mean']) / frozen['predictor_sample_sd']
    scores = {name: standardized[:, :len(beta)] @ beta for name, beta in slopes.items()}
    clinical = pd.read_csv(clin_path, index_col=0)
    assert clinical.index.is_unique and len(clinical) == 60 and features.index.equals(expr.columns)
    clinical = clinical.loc[features.index]
    assert clinical.response.value_counts().to_dict() == {0:35, 1:14, -1:6}
    assert set(clinical.treatment) == {'pre'} and set(clinical.data_set) == {'PUCH'}
    pos = np.flatnonzero(clinical.response.eq(1)); neg = np.flatnonzero(clinical.response.ne(1))
    rng = np.random.default_rng(puch.SEED)
    ip = rng.integers(len(pos), size=(puch.N_BOOT,len(pos)))
    inn = rng.integers(len(neg), size=(puch.N_BOOT,len(neg)))
    estimates = {}
    max_point_error, max_draw_error = 0., 0.
    for name, score in scores.items():
        point, draws = transport.auc(score, pos, neg, ip, inn)
        point_u, draws_u = verify_rank_auc(score, pos, neg, ip, inn)
        max_point_error = max(max_point_error, abs(point-point_u))
        max_draw_error = max(max_draw_error, float(np.max(abs(draws-draws_u))))
        estimates[name] = (point, draws)
    assert max_point_error < 1e-12 and max_draw_error < 1e-12
    estimates['joint_minus_IFN4'] = tuple(a-b for a,b in zip(estimates['IFN4+APM6'], estimates['IFN4']))
    rows = []
    for name, (point, draws) in estimates.items():
        lo, hi = transport.interval(draws)
        rows.append(dict(model_or_contrast=name, n=len(clinical), responders=len(pos), auc=point, auc_lo=lo, auc_hi=hi))
    table = pd.DataFrame(rows)
    original = pd.read_csv(package/'external_PUCH/extension_v21/results/PUCH_frozen_model_performance.csv')
    original_ifn = original.loc[original.model_or_contrast.eq('IFN4_rank'), ['auc','auc_lo','auc_hi']]
    np.testing.assert_allclose(table.loc[table.model_or_contrast.eq('IFN4'), ['auc','auc_lo','auc_hi']].to_numpy(),
                               original_ifn.to_numpy(), rtol=0, atol=1e-12)
    table.to_csv(results/'external_cohort_intercept_rank_auc.csv', index=False, float_format='%.15g')
    verification = {
        'status':'passed', 'patients':len(clinical), 'responders':len(pos), 'bootstrap_draws':puch.N_BOOT,
        'outcome_stratified_bootstrap_seed':puch.SEED, 'scope':'Pairwise AUC vs Mann–Whitney U on raw linear score ranks and each fixed bootstrap resample.',
        'max_abs_point_error':max_point_error,'max_abs_bootstrap_draw_error':max_draw_error,
        'IFN4_auc_and_interval_reproduce_original_common_intercept_experiment':True,
        'external_parameters_estimated':0}
    (results/'external_rank_verification.json').write_text(json.dumps(verification, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('package-root','input-dir','results'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    run(args.package_root.resolve(), args.input_dir, args.results)
