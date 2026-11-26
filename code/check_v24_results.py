#!/usr/bin/env python3
"""Recompute v24 reference summaries and model contrasts from aggregate files."""
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COHORTS = ("GSE91061", "GSE78220", "Liu2019", "Gide2019")


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def check_random_references():
    result_dir = ROOT / "tcga_panels/results"
    comparisons = read_rows(result_dir / "panel_reference_comparison.csv")
    observed_rows = read_rows(result_dir / "common_panel_observed_auc.csv")
    observed = {(r["panel"], r["cohort"]): float(r["observed_auc"]) for r in observed_rows}
    assert len(comparisons) == 12 and len(observed) == 20
    with np.load(result_dir / "reference_draws.npz", allow_pickle=False) as draws:
        for r in comparisons:
            prefix = f"{r['pool']}_n{r['panel_genes']}_{r['reference_signs']}_"
            observed_auc = float(r["observed_mean_auc"])
            four_aucs = [observed[r["panel"], cohort] for cohort in COHORTS]
            np.testing.assert_allclose(np.mean(four_aucs), observed_auc, atol=1e-14, rtol=0)
            pooled = np.mean([draws[prefix + c] for c in COHORTS], axis=0)
            np.testing.assert_array_equal(pooled, draws[prefix + "mean"])
            number = int(r["n_draws"])
            assert len(pooled) == number == 10000
            greater_equal = np.count_nonzero(pooled >= observed_auc)
            ties = np.count_nonzero(pooled == observed_auc)
            assert greater_equal == int(r["draws_ge_observed"]) and ties == int(r["ties"])
            estimates = [*np.quantile(pooled, [.025, .5, .975]),
                         100 * (number - greater_equal + .5 * ties) / number,
                         (1 + greater_equal) / (1 + number)]
            columns = ["random_q025", "random_median", "random_q975", "percentile", "p_upper"]
            np.testing.assert_allclose(estimates, [float(r[k]) for k in columns], atol=1e-12, rtol=0)
    return len(comparisons)


def check_ridge_contrasts():
    source = ROOT / "external_PUCH"
    reference = read_rows(source / "extension_v21/results/PUCH_frozen_model_performance.csv")
    rows = read_rows(source / "penalty_sensitivity/ridge_penalty_sensitivity.csv")
    data = {(float(r["penalty"]), r["model_or_contrast"]): r for r in rows}
    assert len(data) == len(rows) == 9 and len(reference) == 3
    for r in reference:
        comparable = data[1., r["model_or_contrast"]]
        for key in ("n", "responders", "auc", "auc_lo", "auc_hi", "brier", "brier_lo", "brier_hi"):
            np.testing.assert_allclose(float(comparable[key]), float(r[key]), atol=1e-12, rtol=0)
    for penalty in (.1, 1., 10.):
        single = data[penalty, "IFN4_rank"]
        joint = data[penalty, "IFN4_APM6_rank"]
        change = data[penalty, "joint_minus_IFN4"]
        for metric in ("auc", "brier"):
            difference = float(joint[metric]) - float(single[metric])
            np.testing.assert_allclose(difference, float(change[metric]), atol=1e-12, rtol=0)
    return len(rows)


if __name__ == "__main__":
    print(json.dumps({"two_pool_summary_rows_checked": check_random_references(),
                      "external_model_sensitivity_rows_checked": check_ridge_contrasts(),
                      "passed": True}, indent=2))
