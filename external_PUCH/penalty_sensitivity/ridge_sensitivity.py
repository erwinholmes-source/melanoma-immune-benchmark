#!/usr/bin/env python3
"""Check three fixed ridge penalties using the archived frozen-model pipeline."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


class FixedPenalty:
    """Reuse the archived unit-penalty fitter through an exact change of variables."""

    def __init__(self, core, penalty):
        self.core = core
        self.penalty = penalty

    def ridge_fit(self, x, y, w, allow_fallback):
        scale = np.sqrt(self.penalty)
        # gamma = sqrt(lambda) * beta converts lambda * ||beta||^2 to ||gamma||^2.
        coefficients, diag = self.core.ridge_fit(x / scale, y, w, allow_fallback)
        coefficients[1:] /= scale
        residuals = w * (self.core.expit(coefficients[0] + x @ coefficients[1:]) - y)
        gradients = np.r_[residuals.sum(), x.T @ residuals + self.penalty * coefficients[1:]]
        diag['grad_max'] = np.abs(gradients).max()
        assert diag['grad_max'] < 1e-6
        return coefficients, diag


def run(package, workspace, input_dir, output):
    extension = package / 'external_PUCH/extension_v21'
    sys.path[:0] = [str(package / 'code'), str(package / 'external_PUCH'), str(extension)]
    import a11_core_original as core
    import puch_external_check as puch
    import puch_extension as transport
    import run_increment_public as prior

    reference = json.loads((extension / 'results/frozen_models.json').read_text())
    universe_path = package / 'results/reference/universe_historical733.csv'
    assert transport.sha256(universe_path) == reference['rank_universe_sha256']
    universe = pd.read_csv(universe_path).gene.tolist()
    models = {penalty: transport.fit_original(workspace, package, universe, reference['gene_sets'],
                                               prior, FixedPenalty(core, penalty))
              for penalty in (0.1, 1.0, 10.0)}
    for key in ('training_files', 'hgnc_sha256', 'cohorts', 'training_patients', 'responders'):
        assert reference[key] == models[1.0][key]
    for key in ('predictor_mean', 'predictor_sample_sd'):
        np.testing.assert_allclose(models[1.0][key], reference[key], rtol=0, atol=1e-12)
    for name, model in reference['models'].items():
        for key in ('slopes', 'intercept'):
            np.testing.assert_allclose(models[1.0]['models'][name][key], model[key], rtol=0, atol=1e-12)

    expr_path = input_dir / 'mel_puch_exp_data.csv'
    clinical_path = input_dir / 'mel_puch_cli_data.csv'
    assert puch.digest(expr_path)['git_blob_sha1'] == puch.EXP_BLOB
    assert puch.digest(clinical_path)['git_blob_sha1'] == puch.CLI_BLOB
    expression = pd.read_csv(expr_path, index_col=0)
    clinical = pd.read_csv(clinical_path, index_col=0).loc[expression.columns]
    assert clinical.response.value_counts().to_dict() == {0: 35, 1: 14, -1: 6}
    y = clinical.response.eq(1).to_numpy(dtype=int)
    pos, neg = np.flatnonzero(y), np.flatnonzero(1-y)
    rng = np.random.default_rng(puch.SEED)
    ip = rng.integers(len(pos), size=(puch.N_BOOT, len(pos)))
    inn = rng.integers(len(neg), size=(puch.N_BOOT, len(neg)))
    indices = np.c_[pos[ip], neg[inn]]
    rows = []
    for penalty, fitted in models.items():
        predictions = transport.external_predictions(fitted, expression, universe, core)
        metrics = {}
        for name, probabilities in predictions.items():
            losses = (y - probabilities) ** 2
            auc, auc_draws = transport.auc(probabilities, pos, neg, ip, inn)
            metrics[name] = (np.array([auc, losses.mean()]), np.array([auc_draws, losses[indices].mean(axis=1)]))
        first, joint = predictions.keys()
        metrics['joint_minus_IFN4'] = tuple(a-b for a, b in zip(metrics[joint], metrics[first]))
        for name, (point, draws) in metrics.items():
            row = dict(penalty=penalty, model_or_contrast=name, n=len(y), responders=len(pos))
            for column, estimate, samples in zip(('auc', 'brier'), point, draws):
                row[column] = estimate
                row[column+'_lo'], row[column+'_hi'] = transport.interval(samples)
            if name in fitted['models']:
                spec = fitted['models'][name]
                row['intercept'], row['beta_IFN4'] = spec['intercept'], spec['slopes'][0]
                if len(spec['slopes']) > 1:
                    row['beta_APM6'] = spec['slopes'][1]
                row['max_abs_gradient'] = spec['max_abs_gradient']
            rows.append(row)
    results = pd.DataFrame(rows)
    original = pd.read_csv(extension / 'results/PUCH_frozen_model_performance.csv').set_index('model_or_contrast')
    rebuilt = results[results.penalty.eq(1)].set_index('model_or_contrast').loc[original.index, original.columns]
    np.testing.assert_allclose(rebuilt, original, atol=1e-12, rtol=0)
    results.to_csv(output, index=False, float_format='%.15g')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('package-root', 'workspace', 'input-dir', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    run(args.package_root.resolve(), args.workspace, args.input_dir, args.output)
