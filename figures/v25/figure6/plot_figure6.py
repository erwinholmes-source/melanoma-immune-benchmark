"""Generate Figure 6 from archived cohort summaries and aggregate PUCH ROC vertices."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COHORTS = [("GSE91061", "Riaz", 49), ("GSE78220", "Hugo", 26),
           ("Liu2019", "Liu", 102), ("Gide2019", "Gide", 72)]
PANELS = [("ICR-20", "ICR"), ("GEP-18", "GEP"), ("Chemokine-12", "Chemokine"),
          ("TLS-9", "TLS"), ("APM-7", "APM"), ("Four-gene", "Four-gene")]
INK, TEAL, BROWN, GRAY = "#243b46", "#246d7b", "#9c591f", "#63727a"
FAINT = "#dce2e5"


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def caption(fig, x, y, letter, title):
    fig.text(x, y, letter, fontsize=12, fontweight="bold", color=INK)
    fig.text(x + .031, y, title, fontsize=9.5, fontweight="bold", color=INK)


def axes_style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for edge in ("bottom", "left"):
        ax.spines[edge].set_color(GRAY)
        ax.spines[edge].set_linewidth(.65)
    ax.tick_params(labelsize=7.4, color=GRAY, labelcolor=INK, length=2.6, width=.65)


def interval(ax, low, high, estimate, y, color, filled=True):
    assert low <= high
    ax.hlines(y, low, high, color=color, lw=1.25, zorder=2)
    ax.vlines([low, high], y - .087, y + .087, color=color, lw=.85)
    ax.plot(estimate, y, "o", ms=4.2, markeredgewidth=.9, color=color,
            markerfacecolor=color if filled else "white", zorder=3)


def plot_within(fig, rows):
    caption(fig, .024, .955, "a", "Within-cohort change in AUC after adding a gene set to IFN-4")
    lookup = {(row["cohort"], row["panel"]): row for row in rows}
    assert len(lookup) == len(rows) == 24
    for j, (key, label, n) in enumerate(COHORTS):
        ax = fig.add_axes([.142 + j * .211, .651, .19, .227])
        ax.axvline(0, color=GRAY, lw=.8, ls=(0, (3, 2)))
        for i, (panel, _) in enumerate(PANELS):
            row = lookup[key, panel]
            assert int(row["n"]) == n
            lo, hi, est = (float(row[k]) for k in ("delta_auc_low", "delta_auc_high", "delta_auc"))
            color = BROWN if panel == "APM-7" else GRAY if panel == "Four-gene" else TEAL
            interval(ax, lo, hi, est, 5 - i, color)
        ax.set_ylim(-.55, 5.55)
        ax.set_xlim(-.25, .75)
        ax.set_xticks([-.2, 0, .3, .6], ["−0.2", "0", "0.3", "0.6"])
        ax.set_yticks(range(5, -1, -1), [s for _, s in PANELS] if j == 0 else [""] * 6)
        ax.set_title(f"{label} (n = {n})", fontsize=8.3, fontweight="bold", color=INK, pad=5)
        axes_style(ax)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=4)
    fig.text(.555, .589, "Held-out AUC difference (augmented minus IFN-4)",
             fontsize=8, ha="center", color=INK)


def plot_roc(fig, rows, baseline):
    caption(fig, .024, .526, "b", "External ROC curves in PUCH")
    ax = fig.add_axes([.110, .143, .333, .314])
    ax.set_aspect("equal", adjustable="box")
    ax.plot([0, 1], [0, 1], color="#adb7bb", ls=(0, (2, 2)), lw=.8, zorder=1)
    models = [("IFN4_rank", "IFN-4", TEAL, "-", 1.8),
              ("IFN4_APM6_rank", "IFN-4 + APM6", BROWN, (0, (3, 1.8)), 1.35)]
    for key, label, color, style, width in models:
        points = [r for r in rows if r["model"] == key]
        x = np.array([float(r["false_positive_rate"]) for r in points])
        y = np.array([float(r["true_positive_rate"]) for r in points])
        assert np.all(np.diff(x) >= 0) and np.all(np.diff(y) >= 0)
        assert (x[0], y[0], x[-1], y[-1]) == (0, 0, 1, 1)
        area = float(np.trapezoid(y, x))
        np.testing.assert_allclose(area, float(baseline[key]["auc"]), atol=1e-12, rtol=0)
        ax.plot(x, y, lw=width, ls=style, color=color, label=f"{label}  (AUC {area:.3f})",
                solid_capstyle="butt", dash_capstyle="butt", zorder=3)
    ax.set(xlim=(-.02, 1.02), ylim=(-.02, 1.02), xlabel="False-positive rate", ylabel="True-positive rate")
    ax.set_xticks([0, .5, 1], ["0", "0.5", "1"])
    ax.set_yticks([0, .5, 1], ["0", "0.5", "1"])
    axes_style(ax)
    ax.legend(loc="upper left", bbox_to_anchor=(-.11, -.255), fontsize=7.5, frameon=False,
              handlelength=2.1, borderpad=0, labelspacing=.48, borderaxespad=0, labelcolor=INK)


def plot_sensitivity(fig, rows, key, letter, bottom, xticks, ticklabels, limits):
    caption(fig, .482, bottom + .159, letter,
            "External ΔAUC by ridge penalty" if key == "auc" else "External ΔBrier by ridge penalty")
    ax = fig.add_axes([.576, bottom, .379, .127])
    ax.axvline(0, color=GRAY, lw=.85, ls=(0, (3, 2)), zorder=1)
    indexed = {float(r["penalty"]): r for r in rows if r["model_or_contrast"] == "joint_minus_IFN4"}
    assert set(indexed) == {.1, 1., 10.}
    for y, penalty in zip([2, 1, 0], [.1, 1., 10.]):
        row = indexed[penalty]
        assert int(row["n"]) == 55 and int(row["responders"]) == 14
        interval(ax, float(row[key + "_lo"]), float(row[key + "_hi"]), float(row[key]),
                 y, TEAL, filled=penalty == 1.)
    ax.set(ylim=(-.5, 2.5), xlim=limits)
    ax.set_xticks(xticks, ticklabels)
    ax.set_yticks([2, 1, 0], ["λ = 0.1", "λ = 1", "λ = 10"])
    axes_style(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=5)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--archive", type=Path, required=True, help="Additional file 3 root")
    cli.add_argument("--roc", type=Path, default=Path(__file__).resolve().parent / "results/PUCH_roc_vertices.csv")
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    extension = args.archive / "external_PUCH"
    inc = read_rows(args.archive / "results/increment/increment_all_results.csv")
    sensitivity = read_rows(extension / "penalty_sensitivity/ridge_penalty_sensitivity.csv")
    baseline = {r["model_or_contrast"]: r for r in
                read_rows(extension / "extension_v21/results/PUCH_frozen_model_performance.csv")}
    rc = {"font.family": "DejaVu Sans", "font.size": 8.2, "text.color": INK,
          "axes.labelcolor": INK, "pdf.fonttype": 42, "savefig.facecolor": "white"}
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(7.2, 5.64))
        plot_within(fig, inc)
        plot_roc(fig, read_rows(args.roc), baseline)
        plot_sensitivity(fig, sensitivity, "auc", "c", .367,
                         [-.06, -.04, -.02, 0, .02], ["−0.06", "−0.04", "−0.02", "0", "0.02"], (-.061, .023))
        plot_sensitivity(fig, sensitivity, "brier", "d", .107,
                         [-.002, 0, .002, .004], ["−0.002", "0", "0.002", "0.004"], (-.0031, .0043))
        args.output.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output / "Fig6.pdf", metadata={"Title": "Incremental prediction in melanoma cohorts and external transfer to PUCH"})
        fig.savefig(args.output / "Fig6.png", dpi=450)
        plt.close(fig)


if __name__ == "__main__":
    main()
