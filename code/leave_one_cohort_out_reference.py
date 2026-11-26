"""Recompute five primary signatures' reference percentiles with one cohort omitted.

Usage: python code/leave_one_cohort_out_reference.py ARCHIVE_OR_EXTRACTED_ROOT OUTPUT.tsv
Inputs are the existing 10,000 matched draws by cohort. No model is refitted.
"""
import argparse
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

COHORT_NAMES = {'GSE91061': 'Riaz', 'GSE78220': 'Hugo', 'Liu2019': 'Liu', 'Gide2019': 'Gide'}


def summarize(observed, simulated, weights=None):
    threshold = np.average(observed, weights=weights)
    draws = np.average(simulated, axis=0, weights=weights)
    assert len(draws) == 10_000
    exceeded = np.count_nonzero(draws >= threshold)
    tied = np.count_nonzero(draws == threshold)
    return {
        'mean_auc': threshold,
        'median_random_auc': np.median(draws),
        'percentile': 100 * (len(draws) - exceeded + 0.5 * tied) / len(draws),
        'draws_ge': exceeded,
        'ties': tied,
        'nominal_p_upper': (1 + exceeded) / (1 + len(draws)),
    }


def main(archive, output):
    if archive.is_dir():
        tables = archive / 'tcga_panels/results'
        observations = pd.read_csv(tables / 'common_panel_observed_auc.csv')
        reference = pd.read_csv(tables / 'panel_reference_comparison.csv')
        packed_draws = (tables / 'reference_draws.npz').read_bytes()
    else:
        with ZipFile(archive) as files:
            root, = [s.removesuffix('results/reference_draws.npz') for s in files.namelist()
                     if s.endswith('tcga_panels/results/reference_draws.npz')]
            observations = pd.read_csv(files.open(root + 'results/common_panel_observed_auc.csv'))
            reference = pd.read_csv(files.open(root + 'results/panel_reference_comparison.csv'))
            packed_draws = files.read(root + 'results/reference_draws.npz')
    rows = []
    with np.load(BytesIO(packed_draws)) as draws:
        for panel, frame in observations.groupby('panel', sort=False):
            frame = frame.set_index('cohort').loc[list(COHORT_NAMES)]
            size, = frame.genes.unique()
            mode = 'one_negative' if panel == 'GEP-18' else 'all_positive'
            for pool in ['original142', 'tcga242']:
                scores = np.stack([draws[f'{pool}_n{size}_{mode}_{c}'] for c in frame.index])
                for omit in [None, *frame.index]:
                    keep = frame.index != omit
                    values = summarize(frame.observed_auc.to_numpy()[keep], scores[keep])
                    if omit is None:
                        original = reference.query('panel == @panel and pool == @pool and reference_signs == @mode').iloc[0]
                        np.testing.assert_allclose(values['percentile'], original.percentile, rtol=0, atol=1e-10)
                        assert values['draws_ge'] == original.draws_ge_observed
                    rows.append(dict(panel=panel, common_genes=size, pool=pool,
                                     cohort_weighting='equal', omitted=COHORT_NAMES.get(omit, 'none'),
                                     cohorts_used=sum(keep), **values))
                weighted = summarize(frame.observed_auc.to_numpy(), scores, weights=frame.patients.to_numpy())
                rows.append(dict(panel=panel, common_genes=size, pool=pool,
                                 cohort_weighting='response-evaluable patient count', omitted='none',
                                 cohorts_used=len(frame), **weighted))
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False, sep='\t', float_format='%.9f')


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('archive', type=Path)
    cli.add_argument('output', type=Path)
    args = cli.parse_args()
    main(args.archive, args.output)
