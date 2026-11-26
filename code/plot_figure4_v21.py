"""Plot the four gene-reference distributions and APM increments in Figure 4."""
import argparse
import csv
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--increment", type=Path, required=True)
    parser.add_argument("--tcga", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    archive = np.load(args.reference / "random_reference_auc.npz")
    summaries = {
        r["background"]: r
        for path in [args.reference / "background_comparison_summary.csv",
                     args.tcga / "tcga_random_reference_summary.csv"]
        for r in rows(path) if r["cohort"] == "mean_independent4"
    }
    tcga_draws = np.array([
        float(r["mean_independent4"])
        for r in rows(args.tcga / "draws_auc_tcga_cd8a.csv.gz")
    ])
    apm = {
        r["cohort"]: r for r in rows(args.increment / "table4_incremental_prediction.csv")
        if r["panel"] == "APM-7"
    }
    pools = [
        ("broad_protein_coding", "Shared protein-coding genes", "#7796A8"),
        ("historical733", "Shared immune-panel genes", "#5295A6"),
        ("inflamed142", "CD8A-correlated in study cohorts", "#287F77"),
        ("tcga_cd8a", "CD8A-correlated in TCGA", "#5A6489"),
    ]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": False, "axes.linewidth": .65})
    fig = plt.figure(figsize=(7.2, 4.65), facecolor="white")
    ink = "#263740"
    muted = "#41515A"
    observed_color = "#AA3A32"
    fig.text(.020, .944, "a", size=12, weight="bold", color=ink)
    fig.text(.058, .944, "Reference pools for four-gene sets", size=9.5, weight="bold", color=ink)
    fig.text(.544, .944, "b", size=12, weight="bold", color=ink)
    fig.text(.585, .944, "Added prediction from APM", size=9.5, weight="bold", color=ink)
    fig.legend(
        handles=[Line2D([], [], color=observed_color, lw=1.55, label="Four-gene score (0.621)"),
                 Line2D([], [], color=ink, lw=.95, ls="--", label="Pool median")],
        loc="upper left", bbox_to_anchor=(.055, .911), borderaxespad=0,
        frameon=False, ncol=2, handlelength=1.4, columnspacing=1.2, fontsize=7.2
    )
    bins = np.linspace(.30, .72, 71)
    for index, (key, label, color) in enumerate(pools):
        s = summaries[key]
        values = tcga_draws if key == "tcga_cd8a" else archive[f"{key}__mean_independent4"]
        obs, median, percentile = map(float, (s["panel_auc"], s["random_median"], s["panel_percentile"]))
        assert len(values) == 10000 and bins[0] <= values.min() <= values.max() <= bins[-1]
        assert abs(np.median(values) - median) < 1e-12 and abs(obs - .6211395250979418) < 1e-12
        ties = np.isclose(values, obs, atol=1e-12, rtol=0)
        rank = 100 * (np.sum((values < obs) & ~ties) + .5 * ties.sum()) / len(values)
        assert abs(rank - percentile) < 1e-8
        bottom = .677 - index * .175
        ax = fig.add_axes([.058, bottom, .316, .108])
        ax.hist(values, bins=bins, color=color, edgecolor="white", linewidth=.23)
        ax.axvline(median, color=ink, lw=.97, ls="--")
        ax.axvline(obs, color=observed_color, lw=1.55)
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0, ax.get_ylim()[1] * 1.075)
        ax.set_yticks([])
        ax.set_xticks([.3, .4, .5, .6, .7])
        ax.tick_params(axis="x", labelsize=7.4, pad=2, length=2, width=.6,
                       labelbottom=(index == len(pools) - 1))
        if index == len(pools) - 1:
            ax.set_xlabel("Mean response AUC in four primary cohorts", fontsize=7.7, labelpad=4)
        fig.text(.058, bottom + .117, f"{label} ({int(s['genes']):,})",
                 size=7.9, weight="bold", color=ink)
        fig.text(.392, bottom + .084, f"Median  {median:.3f}", color=muted, size=7.5, va="top")
        fig.text(.392, bottom + .049, f"Rank  {percentile:.2f}%", color=muted, size=7.5, va="top")
    ax = fig.add_axes([.633, .152, .343, .673])
    cohorts = [("GSE91061", "Riaz"), ("GSE78220", "Hugo"),
               ("Liu2019", "Liu"), ("Gide2019", "Gide")]
    for index, (key, label) in enumerate(cohorts):
        rec = apm[key]
        x, lo, hi = map(float, [rec["delta_auc"], rec["delta_auc_low"], rec["delta_auc_high"]])
        boot = [float(r["delta_auc"]) for r in rows(args.increment / "bootstrap" / key / "APM-7_bootstrap.csv")]
        assert len(boot) == 500 and np.max(np.abs(np.quantile(boot, [.025, .975]) - [lo, hi])) < 1e-12
        assert abs(float(rec["auc_M1"]) - float(rec["auc_M0"]) - x) < 1e-10
        y = 3 - index
        ax.errorbar(x, y, xerr=[[x-lo], [hi-x]], fmt="o", color="#795E37", ms=4.7, lw=1.4,
                    capsize=2.6, zorder=3)
        ax.text(.755, y + .20, f"{x:+.3f} [{lo:+.3f}, {hi:+.3f}]", fontsize=7.15,
                va="bottom", ha="right", color=muted)
        if index < 3:
            ax.axhline(y - .48, color="#E0E6E7", lw=.65)
    ax.axvline(0, lw=.9, ls="--", color="#79878B", zorder=1)
    ax.set_xlim(-.20, .78)
    ax.set_xticks([-.2, 0, .2, .4, .6])
    ax.set_ylim(-.40, 3.61)
    ax.set_yticks([3, 2, 1, 0], [f"{name}\n(n = {apm[key]['n']})" for key, name in cohorts])
    ax.tick_params(axis="y", length=0, labelsize=7.7, pad=3)
    ax.tick_params(axis="x", labelsize=7.4, length=2, pad=2)
    ax.set_xlabel("Held-out ΔAUC after adding APM to IFN-4", fontsize=7.75, labelpad=4)
    fig.text(.058, .031, "TCGA supplied expression for gene selection only.", color=muted, size=7.2)
    fig.text(.585, .031, "Points and 95% patient bootstrap intervals", color=muted, size=7.2)
    args.output.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output / "Fig4.pdf", metadata={"Title": "Reference pools and incremental prediction"})
    fig.savefig(args.output / "Fig4.png", dpi=450)
    plt.close(fig)


if __name__ == "__main__":
    main()
