#!/usr/bin/env python3
"""Rebuild the exploratory PUCH panel and fixed-coefficient results from public inputs."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixed_rank_features(expression, universe, membership):
    if not expression.index.is_unique:
        raise ValueError("Duplicate gene symbols in model matrix")
    resolved = expression.loc[universe].astype(float)
    if not np.isfinite(resolved.to_numpy()).all():
        raise ValueError("Incomplete 733-gene expression vector")
    ranks = resolved.rank(axis=0, method="average")
    # Sum integer/half-integer midranks before scaling so equal rank sums remain exactly tied.
    return pd.DataFrame({key: (ranks.loc[genes].sum(axis=0)/len(genes) - 1)/(len(universe) - 1) - 0.5
                         for key, genes in membership.items()})


def fit_original(workspace, package, universe, membership, prior, core):
    counts = {}
    fingerprints = {}
    xs, ys = [], []
    hgnc = package / "resources/hgnc_complete_set.txt"
    hg = pd.read_csv(hgnc, sep="\t", dtype=str)
    approved = hg.loc[hg.status.eq("Approved") & hg.entrez_id.notna()]
    unique = approved.loc[~approved.symbol.duplicated(keep=False) & ~approved.entrez_id.duplicated(keep=False)]
    symbols = dict(zip(unique.entrez_id.astype(int), unique.symbol))
    for cohort in prior.COHORTS:
        clinical_path = workspace / f"prepared_v7/inputs/bulk/{cohort}/clinical.csv"
        profile_path = workspace / prior.PATHS[cohort]
        clinical = pd.read_csv(clinical_path, index_col=0)
        clinical.index = clinical.index.astype(str)
        assert clinical.index.is_unique
        response = clinical.response.dropna().astype(int)
        profile = pd.read_csv(profile_path, index_col=0)
        if cohort == "GSE91061":
            profile.index = profile.index.astype(int)
            profile = profile.loc[profile.index.isin(symbols)]
            profile.index = profile.index.map(symbols)
        profile.columns = profile.columns.astype(str)
        assert profile.columns.is_unique and set(response.unique()) == {0, 1}
        xs.append(fixed_rank_features(profile.loc[:, response.index], universe, membership))
        ys.append(response.to_numpy())
        counts[cohort] = {"n": len(response), "responders": int(response.sum())}
        fingerprints[cohort] = {"expression_sha256": sha256(profile_path), "clinical_sha256": sha256(clinical_path)}
    assert [item["n"] for item in counts.values()] == [49, 26, 102, 72]
    assert [item["responders"] for item in counts.values()] == [10, 14, 41, 39]
    x = pd.concat(xs, axis=0).to_numpy()
    y = np.concatenate(ys)
    means, sds = x.mean(axis=0), x.std(axis=0, ddof=1)
    z = (x - means) / sds
    models = {}
    for name, width in (("IFN4_rank", 1), ("IFN4_APM6_rank", 2)):
        coefficients, diag = core.ridge_fit(z[:, :width], y, np.ones(len(y)), allow_fallback=False)
        models[name] = {"predictors": list(membership)[:width], "intercept": float(coefficients[0]),
                        "slopes": coefficients[1:].tolist(), "max_abs_gradient": float(diag["grad_max"])}
    return {
        "cohorts": counts, "training_patients": len(y), "responders": int(y.sum()),
        "gene_sets": membership, "rank_universe_n": len(universe), "predictor_mean": means.tolist(),
        "predictor_sample_sd": sds.tolist(), "models": models, "ridge_lambda": 1.0,
        "training_files": fingerprints, "hgnc_sha256": sha256(hgnc),
        "estimation": "One pooled intercept, unit patient weights, fixed ridge slopes, no study indicators",
        "feature": "Within-specimen mean centered percentile rank across the archived 733 genes; ties averaged.",
    }


def external_predictions(model, expression, universe, core):
    x = fixed_rank_features(expression, universe, model["gene_sets"])
    z = (x - np.asarray(model["predictor_mean"])) / np.asarray(model["predictor_sample_sd"])
    return {name: core.expit(spec["intercept"] + z[spec["predictors"]].to_numpy() @ np.array(spec["slopes"]))
            for name, spec in model["models"].items()}


def interval(values):
    return list(map(float, np.quantile(values, [0.025, 0.975])))


def auc(values, pos, neg, ip, inn):
    x, z = values[pos], values[neg]
    comparisons = (x[:, None] > z).astype(float) + 0.5 * (x[:, None] == z)
    return float(comparisons.mean()), comparisons[ip[:, :, None], inn[:, None, :]].mean(axis=(1, 2))


def paired_spearman(x, y, ind):
    # Ordinary paired patient resampling; rank within each replicate to retain ties.
    x_rank = pd.DataFrame(x[ind]).rank(axis=1).to_numpy()
    y_rank = pd.DataFrame(y[ind]).rank(axis=1).to_numpy()
    x_rank -= x_rank.mean(axis=1, keepdims=True)
    y_rank -= y_rank.mean(axis=1, keepdims=True)
    boot = (x_rank * y_rank).sum(axis=1) / np.sqrt((x_rank**2).sum(axis=1) * (y_rank**2).sum(axis=1))
    rho = float(np.corrcoef(pd.Series(x).rank(), pd.Series(y).rank())[0, 1])
    return rho, interval(boot)


def run(package, workspace, input_dir, output_dir):
    sys.path.insert(0, str(package / "code"))
    sys.path.insert(0, str(package / "external_PUCH"))
    import run_increment_public as prior
    import a11_core_original as core
    import puch_external_check as original
    out = output_dir
    out.mkdir(parents=True, exist_ok=True)
    ref = package / "results/reference"
    coverage = pd.read_csv(ref / "panel_coverage.csv")
    common = {name: rows.fixed_common_genes.drop_duplicates().item().split(";")
              for name, rows in coverage.groupby("panel", sort=False)}
    universe_file = ref / "universe_historical733.csv"
    universe = pd.read_csv(universe_file).gene.tolist()
    assert len(set(universe)) == len(universe) == 733
    model_genes = {"IFN4_rank": list(prior.IFN4), "APM6_rank": common["APM-7"]}
    assert len(model_genes["APM6_rank"]) == 6
    model = fit_original(workspace, package, universe, model_genes, prior, core)
    model["rank_universe_sha256"] = sha256(universe_file)
    model_path = out / "frozen_models.json"
    model_path.write_text(json.dumps(model, indent=2) + "\n")

    # Evaluate the fixed models after saving their predictor scaling and coefficients.
    expr_path = input_dir / "mel_puch_exp_data.csv"
    clinical_path = input_dir / "mel_puch_cli_data.csv"
    assert original.digest(expr_path)["git_blob_sha1"] == original.EXP_BLOB
    assert original.digest(clinical_path)["git_blob_sha1"] == original.CLI_BLOB
    expression = pd.read_csv(expr_path, index_col=0)
    assert expression.index.is_unique and expression.columns.is_unique and expression.shape == (20293, 55)
    predictors = external_predictions(json.loads(model_path.read_text()), expression, universe, core)
    clinical = pd.read_csv(clinical_path, index_col=0)
    assert clinical.index.is_unique and len(clinical) == 60
    clinical = clinical.loc[expression.columns]
    assert set(clinical.treatment) == {"pre"} and set(clinical.data_set) == {"PUCH"}
    assert clinical.response.value_counts().to_dict() == {0:35, 1:14, -1:6}
    y = clinical.response.eq(1).to_numpy().astype(int)
    pos, neg = np.flatnonzero(y), np.flatnonzero(1-y)
    n, reps = len(y), original.N_BOOT
    rng = np.random.default_rng(original.SEED)
    ip = rng.integers(0, len(pos), size=(reps, len(pos)))
    inn = rng.integers(0, len(neg), size=(reps, len(neg)))
    all_indices = np.c_[pos[ip], neg[inn]]
    corr_indices = np.random.default_rng(20261124).integers(0, n, size=(reps, n))

    genesets = {"IFN-4": prior.IFN4, **prior.PANELS}
    pool = sorted({g for genes in genesets.values() for g in genes if g in expression.index})
    expr = expression.loc[pool].T
    assert np.isfinite(expr.to_numpy()).all() and (expr.std(ddof=1)>0).all()
    z = (expr - expr.mean()) / expr.std(ddof=1)
    baseline = z[prior.IFN4].mean(axis=1).to_numpy()
    ref_auc, ref_boot = auc(baseline, pos, neg, ip, inn)
    rows = []
    for name, nominal in genesets.items():
        for version, targets in (("available", [g for g in nominal if g in z.columns]), ("common", common[name])):
            sgn = np.array([-1 if name == "GEP-18" and g == "CD276" else 1 for g in targets])
            scores = (z[targets].to_numpy() * sgn).mean(axis=1)
            val, bs = auc(scores, pos, neg, ip, inn)
            lo, hi = interval(bs)
            dl, dh = interval(bs - ref_boot)
            rho, (rl, rh) = paired_spearman(scores, baseline, corr_indices)
            rows.append(dict(panel=name, version=version, n=n, responders=len(pos), nominal_n=len(nominal),
                             n_used=len(targets), genes_used=";".join(targets),
                             omitted_nominal_genes=";".join(g for g in nominal if g not in targets),
                             auc=val, auc_lo=lo, auc_hi=hi, rho_IFN4=rho, rho_lo=rl, rho_hi=rh,
                             delta_auc_vs_IFN4=val-ref_auc, delta_lo=dl, delta_hi=dh))
    panel = pd.DataFrame(rows)
    panel.to_csv(out / "PUCH_panel_results.csv", index=False, float_format="%.15g")
    archived = pd.read_csv(package / "external_PUCH/PUCH_aggregate_results.csv").set_index("quantity")
    d = panel.set_index(["panel", "version"])
    for name, source in [("IFN-4", "IFN4_AUC"), ("Four-gene", "four_gene_AUC")]:
        assert np.allclose(d.loc[(name,"available"), ["auc","auc_lo","auc_hi"]].to_numpy(float),
                           archived.loc[source].to_numpy(float), atol=1e-12, rtol=0)
    four = d.loc[("Four-gene","available")]
    assert abs(four.rho_IFN4 - archived.loc["spearman_four_gene_IFN4", "estimate"]) < 1e-12
    assert np.allclose(four[["delta_auc_vs_IFN4","delta_lo","delta_hi"]].to_numpy(float),
                       archived.loc["delta_AUC_four_gene_minus_IFN4"].to_numpy(float), atol=1e-12, rtol=0)

    metrics = {}
    boot = {}
    for name, probabilities in predictors.items():
        point_auc, auc_boot = auc(probabilities, pos, neg, ip, inn)
        losses = (y-probabilities)**2
        metrics[name] = {"auc":point_auc, "brier":float(losses.mean())}
        boot[name] = {"auc":auc_boot, "brier":losses[all_indices].mean(axis=1)}
    baseline_name, combined_name = model["models"]
    metrics["joint_minus_IFN4"] = {k: metrics[combined_name][k] - metrics[baseline_name][k] for k in ("auc", "brier")}
    boot["joint_minus_IFN4"] = {k: boot[combined_name][k] - boot[baseline_name][k] for k in ("auc","brier")}
    transported = []
    for name, pair in metrics.items():
        row = {"model_or_contrast":name, "n":n, "responders":len(pos), **pair}
        for k in ("auc", "brier"):
            row[k+"_lo"], row[k+"_hi"] = interval(boot[name][k])
        transported.append(row)
    pd.DataFrame(transported).to_csv(out / "PUCH_frozen_model_performance.csv", index=False, float_format="%.15g")
    checks = {"external_inputs": {p.name:sha256(p) for p in (expr_path, clinical_path)},
              "available_genes":{k:len(set(v) & set(pool)) for k,v in genesets.items()},
              "reproduced_prior_PUCH_fixed_scores_and_95CI":True, "fixed_rank_reference_coverage":len(universe),
              "numpy":np.__version__, "pandas":pd.__version__, "auc_boot_seed":original.SEED, "rho_boot_seed":20261124,
              "bootstrap_draws":reps, "training_complete_before_loading_external_labels":True,
              "conditional_inference":"AUC/Brier intervals use 5000 paired, outcome-stratified external patient resamples, conditional on scores or fitted models. Spearman intervals use 5000 unstratified paired patient resamples, conditional on scores."}
    (out / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True, help="Directory populated by the supplement's prepare_sources.py")
    parser.add_argument("--input-dir", type=Path, required=True, help="Pinned PUCH expression and clinical CSVs")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.package_root.resolve(), args.workspace.resolve(), args.input_dir.resolve(), args.output.resolve())
