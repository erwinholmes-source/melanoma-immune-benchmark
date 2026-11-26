#!/usr/bin/env python3
"""Check published ROC coordinates, cohort fit summaries and the external ranking contrast."""
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TRAINING = ROOT / "external_PUCH/training_diagnostics"
ORIGINAL = ROOT / "external_PUCH/extension_v21/results"


def read(name):
    with name.open(newline="") as fp:
        return list(csv.DictReader(fp))


def check_curves():
    models = read(ORIGINAL / "PUCH_frozen_model_performance.csv")
    rows = read(ROOT / "figures/v25/figure6/results/PUCH_roc_vertices.csv")
    assert set(rows[0]) == {"model", "false_positive_rate", "true_positive_rate"}
    areas = {}
    for model in models:
        name = model["model_or_contrast"]
        if name == "joint_minus_IFN4":
            continue
        points = [r for r in rows if r["model"] == name]
        x, y = (np.array([float(r[col]) for r in points])
                for col in ("false_positive_rate", "true_positive_rate"))
        np.testing.assert_array_equal([x[0], y[0], x[-1], y[-1]], [0, 0, 1, 1])
        assert np.all(np.diff(x) >= 0) and np.all(np.diff(y) >= 0)
        area = np.trapezoid(y, x)
        np.testing.assert_allclose(area, float(model["auc"]), atol=1e-12, rtol=0)
        areas[name] = float(area)
    assert len(areas) == 2 and {r["model"] for r in rows} == set(areas)
    return areas


def check_training():
    frozen = json.loads((ORIGINAL / "frozen_models.json").read_text())
    summaries = {r["cohort"]: r for r in read(TRAINING / "cohort_rank_summaries.csv")}
    pooled = summaries.pop("Pooled")
    assert set(summaries) == set(frozen["cohorts"])
    for field in ("n", "responders"):
        assert sum(int(r[field]) for r in summaries.values()) == int(pooled[field])
        assert all(int(r[field]) == frozen["cohorts"][c][field] for c, r in summaries.items())
    fitted = read(TRAINING / "training_model_coefficients.csv")
    assert len({(r["intercept"], r["weighting"], r["model"]) for r in fitted}) == len(fitted) == 12
    common = {r["model"]: r for r in fitted
              if r["intercept"] == "common" and r["weighting"] == "equal_patients"}
    for key, name in [("IFN4", "IFN4_rank"), ("IFN4+APM6", "IFN4_APM6_rank")]:
        fit = frozen["models"][name]
        coeffs = [float(common[key]["beta_IFN4"])]
        if key == "IFN4+APM6":
            coeffs.append(float(common[key]["beta_APM6"]))
        np.testing.assert_allclose([float(common[key]["intercept_common"]), *coeffs],
                                   [fit["intercept"], *fit["slopes"]], atol=1e-8, rtol=0)
    omit = read(TRAINING / "leave_one_cohort_out.csv")
    assert len(omit) == 4 and {r["omitted"] for r in omit} == set(summaries)
    for row in omit:
        for field in ("n", "responders"):
            assert int(row[field]) == int(pooled[field]) - int(summaries[row["omitted"]][field])
    return len(fitted), len(omit)


def check_external_ranking(areas):
    rows = {r["model_or_contrast"]: r for r in read(TRAINING / "external_cohort_intercept_rank_auc.csv")}
    assert set(rows) == {"IFN4", "IFN4+APM6", "joint_minus_IFN4"}
    assert all(int(r["n"]) == 55 and int(r["responders"]) == 14 for r in rows.values())
    np.testing.assert_allclose(float(rows["IFN4"]["auc"]), areas["IFN4_rank"], atol=1e-12, rtol=0)
    delta = float(rows["IFN4+APM6"]["auc"]) - float(rows["IFN4"]["auc"])
    np.testing.assert_allclose(delta, float(rows["joint_minus_IFN4"]["auc"]), atol=1e-12, rtol=0)
    return delta


if __name__ == "__main__":
    areas = check_curves()
    fits, omissions = check_training()
    diff = check_external_ranking(areas)
    print(json.dumps({"passed": True, "primary_external_ROC_AUCs": areas,
                      "training_fits": fits, "training_omissions": omissions,
                      "post_hoc_external_rank_AUC_difference": diff}, indent=2))
