"""Rebuild manuscript Figures 1, 2 and 5 using published aggregate results."""
import argparse
import csv
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
import matplotlib.pyplot as plt
import numpy as np

COHORTS = [
    ("GSE91061", "Riaz"), ("GSE78220", "Hugo"),
    ("Liu2019", "Liu"), ("Gide2019", "Gide"), ("PUCH", "PUCH"),
]
PANELS = [
    ("ICR-20", "ICR", 19), ("GEP-18", "GEP", 17),
    ("Chemokine-12", "Chemokine", 8), ("TLS-9", "TLS", 8),
    ("APM-7", "APM", 6), ("Four-gene", "Four-gene", 4), ("IFN-4", "IFN-4", 4),
]
INK = "#20343D"
MUTED = "#465A63"
TEAL = "#246E78"
PURPLE = "#65548F"


def read_csv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if path is not None:
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def save(fig, outdir, number, title):
    fig.savefig(outdir / f"Fig{number}.pdf", metadata={"Title": title})
    fig.savefig(outdir / f"Fig{number}.png", dpi=450)
    plt.close(fig)


def figure1(outdir):
    fig = plt.figure(figsize=(7.2, 4.82), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()

    def panel(x, y, w, h, title, body, face, edge, title_size=9.1):
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0,rounding_size=0.009",
            fc=face, ec=edge, lw=.7,
        ))
        ax.text(x + .015, y + h - .028, title, weight="bold", size=title_size,
                color=edge, va="top", linespacing=1.2)
        title_lines = title.count("\n") + 1
        yy = y + h - .028 - .039 * title_lines - .009
        ax.text(x + .015, yy, body, fontsize=8.05, color=INK, va="top",
                linespacing=1.38)

    def label(y, txt):
        ax.text(.025, y, txt, fontsize=11, color=INK, weight="bold", va="top")

    def arrow(x1, x2, y):
        ax.annotate("", (x2, y), (x1, y), arrowprops={
            "arrowstyle": "-|>", "color": "#6B7E86", "lw": .8,
            "mutation_scale": 8, "shrinkA": 0, "shrinkB": 0})

    label(.973, "Pretreatment melanoma cohorts")
    panel(.025, .729, .407, .175, "Four primary RNA-seq cohorts",
          "253 baseline patients, 249 response-evaluable\n"
          "Riaz 49 · Hugo 26 · Liu 102 · Gide 72", "#F1F7F7", TEAL)
    panel(.448, .729, .237, .175, "Yale series",
          "104 evaluable patients\nHistorical selection cohort", "#F6F7F8", "#566973")
    panel(.701, .729, .274, .175, "External PUCH series",
          "55 evaluable patients\nScore and model evaluation", "#F6F4F9", PURPLE)
    label(.686, "Gene-set concordance and cellular context")
    panel(.025, .409, .464, .213, "Bulk RNA measurements",
          "Five published lists: ICR, GEP, chemokine, TLS, APM\n"
          "Four-gene comparator and IFN-4 index\n"
          "Biopsy rankings and cohort-specific response AUCs", "#F5F8F8", TEAL)
    panel(.505, .409, .470, .213, "Cellular sources of the four-gene score",
          "Tirosh atlas: transcript detection by lineage\n"
          "Jerby-Arnon atlas: bulk compartment markers\n"
          "Sade-Feldman series: sorted CD45-positive cells", "#FAF7F1", "#78623E")
    label(.365, "Response prediction benchmarks")
    panel(.025, .038, .293, .264, "1   Random references",
          "Compare the mean response AUC\n"
          "with size-matched draws. Define\n"
          "inflammatory genes in the study\n"
          "cohorts and separately in 367\n"
          "TCGA metastatic tumors.", "#F1F7F7", TEAL)
    panel(.354, .038, .292, .264, "2   Increment within cohorts",
          "Add each candidate to IFN-4.\n"
          "Compare held-out AUC and\n"
          "Brier score on the same patients\n"
          "in each primary cohort.", "#F1F7F7", TEAL, 8.8)
    panel(.682, .038, .293, .264, "3   Transfer to PUCH",
          "Fit IFN-4 and IFN-4 plus six\n"
          "APM genes in the primary\n"
          "cohorts. Apply the fixed rank\n"
          "reference, scaling and models\n"
          "to PUCH.", "#F6F4F9", PURPLE)
    arrow(.321, .351, .171)
    arrow(.649, .679, .171)
    save(fig, outdir, 1, "Study cohorts and comparisons")


def compile_cells(archive):
    ref = archive / "results/reference"
    rows_auc = {(r["panel"], r["cohort"]): r for r in read_csv(ref / "scoring_auc.csv")
                if r["method"] == "z_common"}
    rows_rho = {(r["panel"], r["cohort"]): r for r in read_csv(ref / "scoring_correlations.csv")
                if r["method"] == "z_common"}
    members = {(r["panel"], r["cohort"]): r for r in read_csv(ref / "panel_coverage.csv")}
    external = {r["panel"]: r for r in read_csv(
        archive / "external_PUCH/extension_v21/results/PUCH_panel_results.csv") if r["version"] == "common"}
    primary_table = {r["panel"]: r for r in read_csv(archive / "primary_panels/primary_panel_table.csv")}
    cells = []
    for panel, label, size in PANELS:
        puch = external[panel]
        gene_list = set(puch["genes_used"].split(";"))
        assert len(gene_list) == size == int(puch["n_used"])
        for cohort, short in COHORTS:
            if cohort == "PUCH":
                rho = float(puch["rho_IFN4"])
                auc, lo, hi = map(float, [puch["auc"], puch["auc_lo"], puch["auc_hi"]])
                n_base = n_auc = int(puch["n"])
                responders = int(puch["responders"])
            else:
                key = (panel, cohort)
                outcome = rows_auc[key]
                r = rows_rho[key]
                assert int(outcome["genes_used"]) == size
                assert set(members[key]["fixed_common_genes"].split(";")) == gene_list
                rho = float(r["rho_with_IFN4_same_method"])
                auc, lo, hi = map(float, [outcome["auc"],outcome["ci95_low"],outcome["ci95_high"]])
                n_base, n_auc, responders = int(r["n_baseline"]), int(outcome["n"]), int(outcome["responders"])
            assert 0 <= lo <= auc <= hi <= 1 and -1 <= rho <= 1 and n_base >= n_auc > responders > 0
            cells.append(dict(panel=panel, common_genes=size, cohort=short, baseline=n_base,
                              response_evaluable=n_auc, responders=responders, rho_with_IFN4=rho,
                              auc=auc, auc_low=lo, auc_high=hi))
        if panel in primary_table:
            primary = cells[-5:-1]
            assert abs(np.median([r["rho_with_IFN4"] for r in primary]) -
                       float(primary_table[panel]["median_cohort_rho_with_IFN4"])) < 1e-12
            assert abs(np.mean([r["auc"] for r in primary]) -
                       float(primary_table[panel]["mean_cohort_auc"])) < 1e-12
    return cells


def figure2(archive, outdir, source_csv):
    cells = compile_cells(archive)
    write_csv(source_csv, cells)
    fig = plt.figure(figsize=(7.2, 4.32), facecolor="white")
    concordance = colors.LinearSegmentedColormap.from_list(
        "coexpression", ["#F6FAFA", "#C4DCDD", "#71A8AA", "#205E68"])
    aucmap = colors.LinearSegmentedColormap.from_list(
        "response", ["#96432B", "#DFAC92", "#FAF7F0", "#85B8B6", "#226470"])
    cols = [0, 1, 2, 3, 4.18]
    ys = [0, 1, 2, 3, 4, 5.28, 6.28]
    for left, title, metric, n_key, cmap, letter in [
        (.165, "Concordance with IFN-4", "rho_with_IFN4", "baseline", concordance, "a"),
        (.595, "Response discrimination", "auc", "response_evaluable", aucmap, "b")]:
        w = .38
        ax = fig.add_axes([left, .202, w, .57])
        ax.set(xlim=(-.5, cols[-1]+.5), ylim=(ys[-1]+.5, -.5))
        ax.set_axis_off()
        title_start = .034 if letter == "a" else .567
        fig.text(title_start, .955, letter, color=INK, fontsize=12, weight="bold")
        fig.text(title_start + .032, .955, title, color=INK, fontsize=10.1, weight="bold")
        for i, (panel, panel_label, size) in enumerate(PANELS):
            for j, (_, short) in enumerate(COHORTS):
                v = cells[i*5+j][metric]
                self_comparison = metric=="rho_with_IFN4" and panel=="IFN-4"
                face = (.96,.965,.966,1) if self_comparison else cmap(v)
                ax.add_patch(Rectangle((cols[j]-.49, ys[i]-.477), .98, .954,
                                       linewidth=0, facecolor=face))
                luminance = .2126*face[0]+.7152*face[1]+.0722*face[2]
                ax.text(cols[j], ys[i], "–" if self_comparison else f"{v:.3f}",
                        ha="center", va="center", fontsize=8.25,
                        color=("#869397" if self_comparison else ("white" if luminance < .51 else INK)))
            if metric == "rho_with_IFN4":
                row_label = f"{panel_label} ({size})" if i < 5 else panel_label
                ax.text(-.63, ys[i], row_label, ha="right", va="center",
                        clip_on=False, fontsize=8.3, color=INK)
        # Count headings follow the population represented by each panel.
        for j, (_, short) in enumerate(COHORTS):
            ax.text(cols[j], -.72, f"{short}\n{cells[j][n_key]}", va="bottom", ha="center",
                    fontsize=7.65, linespacing=1.25, color=INK, clip_on=False)
        for midpoint, lo, hi, label in [(1.5,-.47,3.47,"Primary cohorts"),
                                       (4.18,3.73,4.63,"External")]:
            ax.text(midpoint,-1.91,label,ha="center",va="bottom",fontsize=7.45,
                    color=MUTED,clip_on=False)
            ax.plot([lo,hi],[-1.83,-1.83],color="#B6C2C6",lw=.65,clip_on=False)
        cax=fig.add_axes([left+.006, .100, w-.012, .025])
        cb=fig.colorbar(ScalarMappable(cmap=cmap,norm=colors.Normalize(0,1)),cax=cax,
                        orientation="horizontal",ticks=[0,.25,.5,.75,1])
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=7.2,length=2,width=.5,pad=2)
        cb.ax.set_xticklabels(["0","0.25","0.50","0.75","1"])
        label="Spearman correlation (ρ)" if metric=="rho_with_IFN4" else "Standalone AUC (chance = 0.50)"
        fig.text(left+w/2,.016,label,ha="center",color=INK,fontsize=7.8)
    save(fig,outdir,2,"Agreement with IFN-4 and response discrimination by cohort")


def four_gene_distributions(archive):
    ref=archive/"results/reference"
    tcga=archive/"external_TCGA_background/results"
    summaries={r["background"]:r for p in
               [ref/"background_comparison_summary.csv", tcga/"tcga_random_reference_summary.csv"]
               for r in read_csv(p) if r["cohort"]=="mean_independent4"}
    with np.load(ref/"random_reference_auc.npz") as npz:
        arrays={key:npz[f"{key}__mean_independent4"].copy() for key in
                ["broad_protein_coding","historical733","inflamed142"]}
    arrays["tcga_cd8a"]=np.array([float(r["mean_independent4"])
                                   for r in read_csv(tcga/"draws_auc_tcga_cd8a.csv.gz")])
    return summaries, arrays


def figure5(archive, panel_results, outdir, source_csv):
    summaries,arrays=four_gene_distributions(archive)
    comparison={(r["panel"],r["pool"],r["reference_signs"]):r for r in read_csv(panel_results)}
    fig=plt.figure(figsize=(7.2,4.65),facecolor="white")
    fig.text(.017,.945,"a",size=12,color=INK,weight="bold")
    fig.text(.052,.945,"Four-gene comparator",size=10.1,color=INK,weight="bold")
    fig.text(.51,.945,"b",size=12,color=INK,weight="bold")
    fig.text(.545,.945,"Published gene sets",size=10.1,color=INK,weight="bold")
    legend_style=[Line2D([],[],color="#A74130",lw=1.6,label="Observed AUC (0.621)"),
                  Line2D([],[],color=INK,lw=1.0,ls="--",label="Pool median")]
    fig.legend(handles=legend_style,loc="upper left",bbox_to_anchor=(.049,.912),
               borderaxespad=0,frameon=False,ncol=2,fontsize=7.25,
               handlelength=1.35,columnspacing=1.25)
    collected=[]
    edges=np.linspace(.30,.72,71)
    for j,(key,title,colour) in enumerate([
        ("broad_protein_coding","Shared protein-coding genes","#7899AA"),
        ("historical733","Shared immune-panel genes","#5295A6"),
        ("inflamed142","CD8A-correlated in study cohorts","#287F77"),
        ("tcga_cd8a","CD8A-correlated in TCGA",PURPLE)]):
        row=summaries[key]
        x=arrays[key]
        obs, median, pct=map(float,[row["panel_auc"],row["random_median"],row["panel_percentile"]])
        ties=np.isclose(x,obs,atol=1e-12,rtol=0)
        emp=100*(np.sum((x<obs)&~ties)+.5*ties.sum())/len(x)
        assert len(x)==10000 and abs(emp-pct)<1e-8 and abs(np.median(x)-median)<1e-12
        assert x.min()>=edges[0] and x.max()<=edges[-1]
        y=.682-j*.185
        ax=fig.add_axes([.052,y,.276,.098])
        ax.hist(x,bins=edges,color=colour,ec="white",lw=.2)
        ax.axvline(median,color=INK,lw=1.0,ls="--")
        ax.axvline(obs,color="#A74130",lw=1.55)
        ax.set(xlim=(edges[0],edges[-1]),yticks=[],xticks=[.3,.4,.5,.6,.7])
        ax.spines[['top','right','left']].set_visible(False)
        ax.tick_params(axis="x",labelsize=7.4,length=2,width=.55,labelbottom=(j==3),pad=2)
        ax.set_ylim(0,ax.get_ylim()[1]*1.04)
        fig.text(.052,y+.108,f"{title} ({int(row['genes']):,})",size=7.5,weight="bold",color=INK)
        fig.text(.344,y+.067,f"Median  {median:.3f}",size=7.4,color=MUTED)
        fig.text(.344,y+.030,f"Rank  {pct:.2f}%",size=7.4,color=MUTED)
        if j==3:ax.set_xlabel("Mean response AUC across four cohorts",labelpad=4,fontsize=7.45)
        collected.append(dict(panel="Four-gene",size=4,pool=key,signs="all_positive",
                              observed_mean_auc=obs,random_median=median,percentile=pct))
    pools=[("original142","Study (142)","#237D78","o"),("tcga242","TCGA (242)",PURPLE,"D")]
    legend=[Line2D([],[],marker=m,ls="",ms=5.0,color=c,label=txt,
                   markerfacecolor=("white" if k=="tcga242" else c),markeredgewidth=.95)
            for k,txt,c,m in pools]
    fig.legend(handles=legend,loc="upper left",bbox_to_anchor=(.544,.912),
               borderaxespad=0,frameon=False,ncol=2,fontsize=7.25,
               handlelength=1.32,columnspacing=1.15)
    ax=fig.add_axes([.63,.178,.207,.560])
    ax.set(xlim=(-2,102),ylim=(-.44,4.5),xticks=[0,25,50,75,100],yticks=[])
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis="x",length=2,width=.6,pad=3,labelsize=7.4)
    ax.set_xlabel("Random-set percentile",fontsize=7.7)
    ax.grid(axis="x",alpha=.39,color="#ACB8BC",lw=.5)
    ax.set_axisbelow(True)
    for k,col in enumerate([.883,.960]):
        fig.text(col,.785,("Study" if k==0 else "TCGA"),ha="center",fontsize=7.3,weight="bold",
                 color=pools[k][2])
    for i,(panel,txt,size) in enumerate(PANELS[:5]):
        y=4-i; ys=[y+.085,y-.085]
        sign="one_negative" if panel=="GEP-18" else "all_positive"
        pct=[]
        for (pool,_,color,marker),yv,col in zip(pools,ys,[.883,.960]):
            row=comparison[(panel,pool,sign)]
            assert int(row["panel_genes"])==size and int(row["n_draws"])==10000
            v=float(row["percentile"]);pct.append(v)
            ax.scatter(v,yv,s=22,marker=marker,facecolor="white" if pool=="tcga242" else color,
                       edgecolor=color,linewidth=.95,zorder=4)
            fig.text(col,.178 + .560*(y+.44)/4.94,f"{v:.2f}",ha="center",va="center",
                     fontsize=7.75,color=INK)
            collected.append(dict(panel=panel,size=size,pool=pool,signs=sign,
                              observed_mean_auc=row["observed_mean_auc"],
                              random_median=row["random_median"],percentile=v))
        ax.plot(pct,ys,color="#A6B1B4",lw=1,zorder=2)
        ax.text(-5,y,f"{txt} ({size})",ha="right",va="center",fontsize=7.55,clip_on=False)
    fig.text(.545,.052,"GEP draws use one negatively signed gene.",fontsize=7.4,color=MUTED)
    write_csv(source_csv,collected)
    save(fig,outdir,5,"Response discrimination against random gene sets")


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--archive",type=Path,required=True,help="Additional file 3 root")
    cli.add_argument("--output",type=Path,required=True)
    cli.add_argument("--tables",type=Path,help="Optional directory for figure source values")
    cli.add_argument("--panel-results",type=Path,
                     help="Override path to the published-panel random-set comparison CSV")
    cli.add_argument("--figures",default="1,2,5",help="Comma-separated figure numbers")
    args=cli.parse_args()
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"pdf.fonttype":42,
                         "axes.labelcolor":INK,"xtick.color":MUTED,"ytick.color":MUTED,
                         "axes.edgecolor":"#84959B","axes.linewidth":.7,
                         "savefig.facecolor":"white"})
    args.output.mkdir(parents=True,exist_ok=True)
    if args.tables:args.tables.mkdir(parents=True,exist_ok=True)
    selected={int(s) for s in args.figures.split(",")}
    if not selected<={1,2,5}:cli.error("Figures 1, 2 and 5 are available.")
    if 1 in selected:figure1(args.output)
    if 2 in selected:figure2(args.archive,args.output,args.tables/"fig2_values.csv" if args.tables else None)
    if 5 in selected:
        panel_results=args.panel_results or args.archive/"tcga_panels/results/panel_reference_comparison.csv"
        figure5(args.archive,panel_results,args.output,args.tables/"fig5_values.csv" if args.tables else None)


if __name__=="__main__":
    main()
