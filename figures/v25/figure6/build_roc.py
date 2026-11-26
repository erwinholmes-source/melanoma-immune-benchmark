"""Reconstruct and verify the two PUCH ROC curves from the archived frozen models."""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def roc_vertices(probabilities, response):
    """Retain endpoints and corners, grouping ties at each distinct probability."""
    # Sort privately in memory. Neither scores nor decision thresholds are exported.
    order = np.argsort(-probabilities, kind="stable")
    scores = probabilities[order]
    labels = response[order]
    last_in_group = np.r_[scores[:-1] != scores[1:], True]
    tp = np.r_[0, labels.cumsum()[last_in_group]]
    fp = np.r_[0, (1 - labels).cumsum()[last_in_group]]
    # Removing interior collinear points preserves the full empirical curve and its area.
    dx, dy = np.diff(fp), np.diff(tp)
    turns = dx[:-1] * dy[1:] != dy[:-1] * dx[1:]
    keep = np.r_[True, turns, True]
    return fp[keep] / fp[-1], tp[keep] / tp[-1]


def build(archive, input_dir, output):
    root = archive.resolve()
    sys.path[:0] = [str(root / "code"), str(root / "external_PUCH"),
                     str(root / "external_PUCH/extension_v21")]
    import a11_core_original as core
    import puch_extension as transfer
    import puch_external_check as sources

    registry = root / "external_PUCH/extension_v21/results"
    fit = json.loads((registry / "frozen_models.json").read_text())
    assert fit["training_patients"] == 249 and float(fit["ridge_lambda"]) == 1.
    universe_file = root / "results/reference/universe_historical733.csv"
    assert transfer.sha256(universe_file) == fit["rank_universe_sha256"]
    universe = pd.read_csv(universe_file).gene.tolist()
    expr_file = input_dir / "mel_puch_exp_data.csv"
    clin_file = input_dir / "mel_puch_cli_data.csv"
    assert sources.digest(expr_file)["git_blob_sha1"] == sources.EXP_BLOB
    assert sources.digest(clin_file)["git_blob_sha1"] == sources.CLI_BLOB
    profiles = pd.read_csv(expr_file, index_col=0)
    clinical = pd.read_csv(clin_file, index_col=0)
    assert profiles.index.is_unique and profiles.columns.is_unique and clinical.index.is_unique
    probabilities = transfer.external_predictions(fit, profiles, universe, core)
    clinical = clinical.loc[profiles.columns]
    assert set(clinical.treatment) == {"pre"} and set(clinical.data_set) == {"PUCH"}
    assert clinical.response.value_counts().to_dict() == {0: 35, 1: 14, -1: 6}
    response = clinical.response.eq(1).to_numpy(dtype=int)
    performance = pd.read_csv(registry / "PUCH_frozen_model_performance.csv").set_index("model_or_contrast")
    rows = []
    checks = {}
    for name, scores in probabilities.items():
        x, y = roc_vertices(scores, response)
        area = np.trapezoid(y, x)
        # Independent Mann–Whitney comparison of each responder/nonresponder pair.
        pair_diff = scores[response == 1, None] - scores[response == 0]
        pair_auc = np.mean((pair_diff > 0) + 0.5 * (pair_diff == 0))
        brier = np.mean((response - scores) ** 2)
        np.testing.assert_allclose([area, pair_auc, brier],
                                   performance.loc[name, ["auc", "auc", "brier"]].to_numpy(float),
                                   atol=1e-12, rtol=0)
        rows.extend((name, fpr, tpr) for fpr, tpr in zip(x, y))
        checks[name] = {"vertices": len(x), "ROC_trapezoidal_AUC": float(area),
                        "pairwise_AUC": float(pair_auc), "Brier": float(brier)}

    output.mkdir(parents=True, exist_ok=True)
    with (output / "PUCH_roc_vertices.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["model", "false_positive_rate", "true_positive_rate"])
        writer.writerows(rows)
    report = {"N": len(response), "responders": int(response.sum()), "ridge_penalty": 1,
              "fixed_rank_universe_genes": len(universe), "models": checks}
    (output / "PUCH_roc_verification.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="Additional file 3 root")
    parser.add_argument("--input-dir", type=Path, required=True, help="Pinned PUCH clinical and expression files")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.archive, args.input_dir, args.output)
