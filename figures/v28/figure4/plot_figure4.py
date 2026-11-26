"""Generate v28 Figure 4 from archived cohort summaries, aggregate PUCH ROC vertices, and cohort-intercept sensitivity."""
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
DIVIDER = "#c9d3d8"


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


def interval(ax, low, high, estimate, y, color, filled=True, marker="o", ms=4.2):
    assert low <= high
    ax.hlines(y, low, high, color=color, lw=1.25, zorder=2)
    ax.vlines([low, high], y - .10, y + .10, color=color, lw=.85)
    ax.plot(estimate, y, marker, ms=ms, markeredgewidth=.95, color=color,
            markerfacecolor=color if filled else "white", zorder=3)


def plot_within(fig, rows):
    caption(fig, .024, .957, "a", "Within-cohort change in AUC after adding a gene set to IFN-4")
    lookup = {(row["cohort"], row["panel"]): row for row in rows}
    assert len(lookup) == len(rows) == 24
    for j, (key, label, n) in enumerate(COHORTS):
        ax = fig.add_axes([.142 + j * .211, .664, .19, .220])
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
    fig.text(.555, .604, "Held-out AUC difference (augmented minus IFN-4)",
             fontsize=8, ha="center", color=INK)


def plot_roc(fig, rows, baseline):
    caption(fig, .024, .542, "b", "External ROC curves in PUCH")
    ax = fig.add_axes([.088, .178, .320, .308])
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
    ax.legend(loc="upper left", bbox_to_anchor=(-.10, -.245), fontsize=7.5, frameon=False,
              handlelength=2.1, borderpad=0, labelspacing=.48, borderaxespad=0, labelcolor=INK)


def plot_auc_sensitivity(fig, rows, cohort_row):
    caption(fig, .462, .542, "c", "External ΔAUC in PUCH")
    ax = fig.add_axes([.635, .348, .332, .165])
    ax.axvline(0, color=GRAY, lw=.85, ls=(0, (3, 2)), zorder=1)
    indexed = {float(r["penalty"]): r for r in rows if r["model_or_contrast"] == "joint_minus_IFN4"}
    assert set(indexed) == {.1, 1., 10.}
    for y, penalty in zip([3.1, 2.1, 1.1], [.1, 1., 10.]):
        row = indexed[penalty]
        assert int(row["n"]) == 55 and int(row["responders"]) == 14
        interval(ax, float(row["auc_lo"]), float(row["auc_hi"]), float(row["auc"]),
                 y, TEAL, filled=penalty == 1.)
    ax.axhline(.48, color=DIVIDER, lw=.75, ls=(0, (2, 2)), zorder=1)
    assert int(cohort_row["n"]) == 55 and int(cohort_row["responders"]) == 14
    interval(ax, float(cohort_row["auc_lo"]), float(cohort_row["auc_hi"]), float(cohort_row["auc"]),
             -.22, BROWN, filled=False, marker="D", ms=4.3)
    ax.set(ylim=(-.85, 3.65), xlim=(-.065, .073))
    ax.set_xticks([-.06, -.03, 0, .03, .06], ["−0.06", "−0.03", "0", "0.03", "0.06"])
    ax.set_yticks([3.1, 2.1, 1.1, -.22],
                  ["λ = 0.1", "λ = 1", "λ = 10", "Cohort intercepts\n(post hoc, λ = 1)"])
    axes_style(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=5)


def plot_brier_sensitivity(fig, rows):
    caption(fig, .462, .258, "d", "External ΔBrier by ridge penalty")
    ax = fig.add_axes([.635, .095, .332, .132])
    ax.axvline(0, color=GRAY, lw=.85, ls=(0, (3, 2)), zorder=1)
    indexed = {float(r["penalty"]): r for r in rows if r["model_or_contrast"] == "joint_minus_IFN4"}
    assert set(indexed) == {.1, 1., 10.}
    for y, penalty in zip([2, 1, 0], [.1, 1., 10.]):
        row = indexed[penalty]
        assert int(row["n"]) == 55 and int(row["responders"]) == 14
        interval(ax, float(row["brier_lo"]), float(row["brier_hi"]), float(row["brier"]),
                 y, TEAL, filled=penalty == 1.)
    ax.set(ylim=(-.5, 2.5), xlim=(-.0031, .0043))
    ax.set_xticks([-.002, 0, .002, .004], ["−0.002", "0", "0.002", "0.004"])
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
    cohort_rows = read_rows(extension / "training_diagnostics/external_cohort_intercept_rank_auc.csv")
    cohort_diff = next(r for r in cohort_rows if r["model_or_contrast"] == "joint_minus_IFN4")
    baseline = {r["model_or_contrast"]: r for r in
                read_rows(extension / "extension_v21/results/PUCH_frozen_model_performance.csv")}
    rc = {"font.family": "DejaVu Sans", "font.size": 8.2, "text.color": INK,
          "axes.labelcolor": INK, "pdf.fonttype": 42, "savefig.facecolor": "white"}
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(7.2, 5.88))
        plot_within(fig, inc)
        plot_roc(fig, read_rows(args.roc), baseline)
        plot_auc_sensitivity(fig, sensitivity, cohort_diff)
        plot_brier_sensitivity(fig, sensitivity)
        args.output.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output / "Fig4.pdf", metadata={"Title": "Incremental prediction in melanoma cohorts and external transfer to PUCH"})
        fig.savefig(args.output / "Fig4.png", dpi=450)
        plt.close(fig)


if __name__ == "__main__":
    main()
