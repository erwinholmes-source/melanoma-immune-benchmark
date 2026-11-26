#!/usr/bin/env python3
"""One-page PDF/PNG from aggregate v19 output only."""
import argparse,sys
from pathlib import Path
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--legacy-library-dir',type=Path);a=p.parse_args()
if a.legacy_library_dir:sys.path.append(str(a.legacy_library_dir))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
R=a.results;O=a.output;O.mkdir(parents=True,exist_ok=True)
summary=pd.read_csv(R/'random_reference_summary.csv').query("cohort=='mean_independent4'").set_index('background')
paired=pd.read_csv(R/'rank_vs_z_paired_auc.csv')
null=np.load(R/'random_reference_auc.npz')
fig=plt.figure(figsize=(13.4,7.65),facecolor='white')
gs=fig.add_gridspec(3,2,left=.075,right=.965,bottom=.255,top=.805,width_ratios=[1.03,1],wspace=.34,hspace=.38)
fig.text(.075,.94,'Reference-background and rank-scoring sensitivity',fontsize=20,weight='bold')
fig.text(.075,.901,'Four pretreatment melanoma RNA-seq cohorts | 249 response-evaluable patients',fontsize=11,color='#44515e')
fig.text(.075,.844,'A  Conditional random-reference distributions',fontsize=12,weight='bold')
fig.text(.582,.844,'B  Paired change in response discrimination',fontsize=12,weight='bold')
colors=['#71899C','#558EAE','#237A73'];labels=['15,561 common protein-coding genes','733 assay-shared genes','142 CD8A-correlated genes']
for i,(bg,label,col) in enumerate(zip(['broad_protein_coding','historical733','inflamed142'],labels,colors)):
    ax=fig.add_subplot(gs[i,0]);v=null[f'{bg}__mean_independent4'];s=summary.loc[bg]
    ax.hist(v,bins=np.linspace(.25,.72,65),density=True,color=col,edgecolor='white',linewidth=.25,alpha=.90)
    ax.axvline(s.panel_auc,color='#B43F3F',lw=1.7,zorder=3)
    ax.axvline(s.random_median,color=col,lw=1.3,ls='--')
    ax.set_xlim(.25,.72);ax.set_ylim(0,ax.get_ylim()[1]*1.65)
    ax.text(.015,.97,label,transform=ax.transAxes,ha='left',va='top',fontsize=10,weight='bold')
    ax.text(.015,.72,f'Median {s.random_median:.3f} | panel percentile {s.panel_percentile:.2f}',transform=ax.transAxes,ha='left',va='top',fontsize=9,color='#35424e')
    ax.set_yticks([]);ax.spines['left'].set_visible(False)
    ax.set_xticks([.3,.4,.5,.6,.7]);ax.tick_params(axis='x',labelsize=9)
    if i<2:ax.set_xticklabels([])
    else:ax.set_xlabel('Mean cohort AUC (equal weight per cohort)',labelpad=7)
    if i==0:ax.text(s.panel_auc+.006,ax.get_ylim()[1]*.87,f'Observed\n{s.panel_auc:.3f}',color='#B43F3F',fontsize=9,ha='left',va='top')
ax=fig.add_subplot(gs[:,1]);cs=['GSE91061','GSE78220','Liu2019','Gide2019'];names=['Riaz (n = 49)','Hugo (n = 26)','Liu (n = 102)','Gide (n = 72)']
for i,(c,n) in enumerate(zip(cs,names)):
    y=3-i
    for pn,offset,col,mark in [('Four-gene',.13,'#26728C','o'),('IFN-4',-.13,'#61784E','s')]:
        row=paired[(paired.cohort==c)&(paired.panel==pn)].iloc[0]
        ax.errorbar(row.rank_minus_z_common_auc,y+offset,xerr=[[row.rank_minus_z_common_auc-row.paired_ci95_low],[row.paired_ci95_high-row.rank_minus_z_common_auc]],fmt=mark,color=col,markersize=6,capsize=3,lw=1.5)
    if i<3:ax.axhline(y-.5,color='#E8EBEE',lw=.8)
ax.axvline(0,color='#66717D',ls='--',lw=1,zorder=0)
ax.set_yticks([3,2,1,0],names);ax.set_ylim(-.55,3.6);ax.set_xlim(-.11,.13)
ax.set_xticks([-.10,-.05,0,.05,.10]);ax.set_xlabel('AUC change: rank minus common-gene z score',labelpad=10)
ax.legend(handles=[Line2D([0],[0],marker='o',color='#26728C',label='Four-gene',lw=1.5),Line2D([0],[0],marker='s',color='#61784E',label='IFN-4',lw=1.5)],loc='upper right',bbox_to_anchor=(1,1.04),ncol=2,frameon=False,fontsize=10)
fig.text(.075,.135,'A: 10,000 random four-gene sets per background; the same gene draws were evaluated in all four cohorts. Dashed lines mark null medians.',fontsize=9,color='#46525E')
fig.text(.075,.108,'B: Fixed common gene membership; within-sample ranks use the same 15,561-gene background. Bars: paired 95% bootstrap intervals (5,000 resamples).',fontsize=9,color='#46525E')
fig.text(.075,.081,'All analyses are exploratory and post hoc. The 142-gene pool is conditioned on CD8A correlation in five historical cohorts. AUC overlap does not establish equivalence.',fontsize=9,color='#46525E')
fig.text(.075,.054,'Source expression transformations were retained, including author-transformed Riaz values; rank scoring is a sensitivity analysis, not a locked assay validation.',fontsize=9,color='#46525E')
fig.savefig(O/'Reference_background_and_rank_sensitivity_v19.pdf',metadata={'Title':'Reference-background and rank-scoring sensitivity','Author':'Analysis of archived melanoma cohorts'})
fig.savefig(O/'Reference_background_and_rank_sensitivity_v19.png',dpi=180)
plt.close(fig)
print('Wrote one-page PDF and PNG')
