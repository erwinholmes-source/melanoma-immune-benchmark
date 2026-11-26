#!/usr/bin/env python3
"""Evaluate a fixed, positively signed gene list on the prepared melanoma cohorts."""
import argparse
import json
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
for module_dir in ("code", "external_PUCH", "external_PUCH/extension_v21",
                   "external_TCGA_background", "tcga_panels"):
    sys.path.insert(0, str(ROOT / module_dir))

import a11_core_original as core
import puch_extension as puch_models
import puch_external_check as puch_sources
import run_increment_public as increment
import run_tcga_background as background
from run_panel_reference import simulation, summarize


def training_inputs(workspace, universe, genes):
    features, outcomes, comparison = [], [], {}
    scoring_genes = increment.IFN4 + genes
    preprocess = increment.prep_factory(np.ones(len(genes)))
    for index, cohort in enumerate(increment.COHORTS):
        prepared = workspace / "prepared_v7/inputs"
        clinical = pd.read_csv(prepared / "bulk" / cohort / "clinical.csv", index_col=0)
        y = clinical.response.dropna().astype(int)
        matrix = pd.read_csv(prepared / "rank" / cohort / "source_symbol_expression.csv.gz", index_col=0)
        if not clinical.index.is_unique or not y.isin([0, 1]).all():
            raise ValueError(f"Invalid patient labels in {cohort}")
        observed = matrix.loc[scoring_genes, y.index].T.to_numpy()
        scores = puch_models.fixed_rank_features(matrix.loc[:, y.index], universe,
                                                 {"IFN4": increment.IFN4, "candidate": genes})
        metrics, *_ = core.paired_cv(observed, y.to_numpy(), np.ones(len(y)),
                                     np.random.SeedSequence([20261105, 11, index, 0]),
                                     preprocess_fn=preprocess)
        comparison[cohort] = {key: metrics[key] for key in increment.METRICS}
        features.append(scores.to_numpy())
        outcomes.append(y.to_numpy())
    return np.concatenate(features), np.concatenate(outcomes), [len(y) for y in outcomes], comparison


def reference_test(workspace, genes):
    audit = pd.read_csv(ROOT / "external_TCGA_background/results/tcga_gene_eligibility.csv")
    pool = sorted(audit.loc[audit.eligible_external, "gene"])
    if len(genes) > len(pool):
        raise ValueError("The candidate exceeds the size of the TCGA inflammatory reference pool")
    matrices, outcomes = background.read_evaluation_inputs(workspace / "prepared_v7",
                                                           sorted(set(pool) | set(genes)))
    generated = simulation(pool, len(genes), matrices, outcomes, background)
    observed = np.asarray([
        background.rank_auc(frame.loc[genes].mean(axis=0).to_numpy()[None, :], outcomes[cohort])[0]
        for cohort, frame in matrices.items()
    ])
    reference = np.stack([generated["all_positive_" + cohort] for cohort in matrices])
    results = {}
    for omitted in [None, *matrices]:
        keep = np.array([cohort != omitted for cohort in matrices])
        value = float(observed[keep].mean())
        result = summarize(reference[keep].mean(axis=0), value)
        results[omitted or "all_cohorts"] = {"observed_mean_auc": value, **result}
    return results


def cohort_intercept_fit(z, y, sizes):
    indicators = np.repeat(np.eye(len(sizes)), sizes, axis=0)
    design = np.column_stack([indicators, z])
    penalty = np.r_[np.zeros(len(sizes)), np.ones(z.shape[1])]
    def objective(b):
        linear = design @ b
        loss = np.sum(np.logaddexp(0, linear) - y * linear) + 0.5 * np.dot(penalty * b, b)
        gradient = design.T @ (core.expit(linear) - y) + penalty * b
        return loss, gradient
    fitted = minimize(objective, np.zeros(design.shape[1]), method="BFGS", jac=True,
                      options={"gtol": 1e-9, "maxiter": 500})
    if not np.isfinite(fitted.x).all() or np.max(np.abs(fitted.jac)) > 1e-6:
        raise ArithmeticError("Cohort-intercept ridge fit did not converge")
    return fitted.x


def puch_inputs(directory, universe, genes):
    paths = [directory / "mel_puch_exp_data.csv", directory / "mel_puch_cli_data.csv"]
    for path, expected in zip(paths, (puch_sources.EXP_BLOB, puch_sources.CLI_BLOB)):
        if puch_sources.digest(path)["git_blob_sha1"] != expected:
            raise ValueError(f"Unexpected source file for {path.name}")
    matrix, clinical = [pd.read_csv(path, index_col=0) for path in paths]
    scores = puch_models.fixed_rank_features(matrix, universe, {"IFN4": increment.IFN4, "candidate": genes})
    y = clinical.loc[scores.index, "response"].eq(1).to_numpy(dtype=int)
    return scores.to_numpy(), y


def model_diagnostics(x, y, sizes, external):
    center, scale = x.mean(axis=0), x.std(axis=0, ddof=1)
    if not np.isfinite(scale).all() or min(scale) == 0:
        raise ValueError("A training rank score is constant or nonfinite")
    train = (x - center) / scale
    result = {"predictor_mean": center.tolist(), "predictor_sample_sd": scale.tolist()}
    for specification in ("common_intercept", "cohort_intercepts"):
        training, performance = {}, {}
        for width, model_name in enumerate(("IFN4", "IFN4_plus_candidate"), 1):
            z = train[:, :width]
            if specification == "common_intercept":
                beta, _ = core.ridge_fit(z, y, np.ones(len(y)), allow_fallback=False)
            else:
                beta = cohort_intercept_fit(z, y, sizes)
            intercepts, slopes = beta[:-width], beta[-width:]
            training[model_name] = {"intercepts": intercepts.tolist(), "slopes": slopes.tolist()}
            if external is not None:
                scores, labels = external
                risk = ((scores - center) / scale)[:, :width] @ slopes
                if specification == "common_intercept":
                    probability = core.expit(risk + intercepts[0])
                    auc, brier = core.weighted_metrics(labels, probability, np.ones(len(labels)))
                    performance[model_name] = {"auc": auc, "brier": brier}
                else:
                    auc = background.rank_auc(risk[None, :], labels.astype(bool))[0]
                    performance[model_name] = {"auc": auc}
        result[specification] = {"training": training}
        if external is not None:
            performance["joint_minus_IFN4"] = {
                key: value - performance["IFN4"][key]
                for key, value in performance["IFN4_plus_candidate"].items()
            }
            result[specification]["PUCH"] = performance
    return result


def run(genes, workspace, puch_dir):
    universe = pd.read_csv(ROOT / "results/reference/universe_historical733.csv").gene.tolist()
    if not genes or len(set(genes)) != len(genes):
        raise ValueError("Supply at least one gene, with no duplicate symbols")
    missing = sorted(set(genes) - set(universe))
    if missing:
        raise ValueError("Genes absent from the fixed 733-gene rank universe: " + ", ".join(missing))
    x, y, sizes, paired = training_inputs(workspace, universe, genes)
    external = puch_inputs(puch_dir, universe, genes) if puch_dir else None
    return {"candidate": genes, "training_cohorts": increment.COHORTS, "evaluable_counts": sizes,
            "tcga_matched_reference": reference_test(workspace, genes), "paired_cv": paired,
            "rank_transfer": model_diagnostics(x, y, sizes, external)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", nargs="+", required=True, help="Distinct positive HGNC symbols in the fixed rank universe")
    parser.add_argument("--workspace", type=Path, required=True, help="Private output of code/prepare_sources.py")
    parser.add_argument("--puch-dir", type=Path, help="Private directory holding the two original PUCH CSVs")
    parser.add_argument("--output", type=Path, required=True, help="Destination for one aggregate JSON file")
    args = parser.parse_args()
    summary = run(args.genes, args.workspace, args.puch_dir)
    serialized = json.dumps(summary, indent=2, allow_nan=False,
                            default=lambda value: value.item() if isinstance(value, np.generic) else value)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", dir=args.output.parent, delete=False, encoding="utf-8") as stream:
        stream.write(serialized + "\n")
    Path(stream.name).replace(args.output)
