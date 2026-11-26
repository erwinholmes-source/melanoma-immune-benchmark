"""Summarize the four original RNA-seq cohorts for the published gene-set table."""
import argparse
import csv
import math
import statistics
from pathlib import Path


HERE = Path(__file__).resolve().parent
COHORTS = ("GSE91061", "GSE78220", "Liu2019", "Gide2019")
PANELS = ("ICR-20", "GEP-18", "Chemokine-12", "TLS-9", "APM-7")


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def fixed_common_rows(rows, panel):
    selected = [r for r in rows if r["panel"] == panel and r["method"] == "z_common"]
    if len(selected) != len(COHORTS) or {r["cohort"] for r in selected} != set(COHORTS):
        raise ValueError(f"Expected one common-genes record per primary cohort for {panel}")
    return selected


def build(reference, nulls, signed_nulls):
    correlation_rows = read_rows(reference / "scoring_correlations.csv")
    auc_rows = read_rows(reference / "scoring_auc.csv")
    matched = {r["measure"]: r for r in read_rows(nulls) if r["null"] == "inflamed_genes"}
    signed = {r["reference_signs"]: r for r in read_rows(signed_nulls) if r["pool"] == "inflamed_genes"}
    results = []
    for panel in PANELS:
        ref = matched[f"{panel} (common)"]
        auc_records = fixed_common_rows(auc_rows, panel)
        rho_records = fixed_common_rows(correlation_rows, panel)
        auc = [float(r["auc"]) for r in auc_records]
        rho = [float(r["rho_with_IFN4_same_method"]) for r in rho_records]
        if {int(r["genes_used"]) for r in auc_records} != {int(ref["genes"])}:
            raise ValueError(f"Inconsistent fixed gene set size for {panel}")
        if not math.isclose(statistics.mean(auc), float(ref["observed_mean_auc"]), abs_tol=1e-10):
            raise ValueError(f"Inconsistent observed mean AUC for {panel}")
        if int(ref["null_size"]) != 10000:
            raise ValueError(f"Unexpected random draw count for {panel}")
        sign_matched = {}
        if panel == "GEP-18":
            if not math.isclose(float(ref["percentile"]), float(signed["all_positive"]["percentile"])):
                raise ValueError("GEP reference percentile does not match its original all-positive null")
            sign_matched = signed["one_negative"]
        results.append({
            "panel": panel,
            "genes_common": int(ref["genes"]),
            "cohorts": len(auc),
            "baseline_n_for_correlations": sum(int(r["n_baseline"]) for r in rho_records),
            "response_evaluable_n_for_auc": sum(int(r["n"]) for r in auc_records),
            "responders": sum(int(r["responders"]) for r in auc_records),
            "median_cohort_rho_with_IFN4": statistics.median(rho),
            "min_cohort_rho_with_IFN4": min(rho),
            "max_cohort_rho_with_IFN4": max(rho),
            "mean_cohort_auc": statistics.mean(auc),
            "min_cohort_auc": min(auc),
            "max_cohort_auc": max(auc),
            "percentile_in_size_matched_random_sets": ref["percentile"],
            "nominal_upper_p": ref["p_upper"],
            "gep_sign_matched_percentile": sign_matched.get("percentile", ""),
            "gep_sign_matched_upper_p": sign_matched.get("p_upper", ""),
        })
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=HERE.parent / "results/reference")
    parser.add_argument("--null-table", type=Path, default=HERE / "a10_random_sets.csv")
    parser.add_argument("--gep-sign-table", type=Path, default=HERE / "a11_gep_signed_null.csv")
    parser.add_argument("--output", type=Path, default=HERE / "primary_panel_table.csv")
    args = parser.parse_args()
    rows = build(args.reference_dir, args.null_table, args.gep_sign_table)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
