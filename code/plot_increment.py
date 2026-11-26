"""Plot the within-cohort increments and separate frozen-model assessment from aggregate results."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    inc=records(args.archive/'results/increment/increment_all_results.csv')
    ext=records(args.archive/'external_PUCH/extension_v21/results/PUCH_frozen_model_performance.csv')
    idx={(r['cohort'],r['panel']):r for r in inc}
    assert len(idx)==24
    cohorts=[('GSE91061','Riaz',49),('GSE78220','Hugo',26),('Liu2019','Liu',102),('Gide2019','Gide',72)]
    panels=[('ICR-20','ICR'),('GEP-18','GEP'),('Chemokine-12','Chemokine'),('TLS-9','TLS'),('APM-7','APM'),('Four-gene','Four-gene')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'pdf.fonttype':42})
    fig=plt.figure(figsize=(7.15,4.9))
    ink='#263B47';blue='#317B89';apm='#805826';grey='#6B767D'
    fig.text(.028,.956,'a',fontweight='bold',fontsize=12,color=ink)
    fig.text(.070,.956,'Within-cohort change in AUC after adding a gene set to IFN-4',fontsize=10.2,fontweight='bold',color=ink)
    for j,(cid,label,n) in enumerate(cohorts):
        ax=fig.add_axes([.145+j*.211,.475,.19,.37])
        ax.axvline(0,color=grey,lw=.8,ls='--')
        for i,(key,name) in enumerate(panels):
            r=idx[cid,key]; y=5-i
            assert int(r['n'])==n and float(r['delta_auc_low'])<=float(r['delta_auc_high'])
            lo,hi,est=[float(r[k]) for k in ('delta_auc_low','delta_auc_high','delta_auc')]
            color=apm if key=='APM-7' else blue if key!='Four-gene' else grey
            ax.hlines(y,lo,hi,color=color,lw=1.25)
            ax.vlines([lo,hi],y-.095,y+.095,lw=.8,color=color)
            ax.plot(est,y,'o',color=color,ms=3.7,clip_on=True)
        ax.set_yticks(range(5,-1,-1),[x[1] for x in panels] if j==0 else ['']*6)
        ax.set_ylim(-.55,5.55);ax.set_xlim(-.25,.75);ax.set_xticks([-.2,0,.3,.6],['−0.2','0','0.3','0.6'])
        ax.set_title(f'{label} (n = {n})',fontsize=8.5,fontweight='bold',pad=5,color=ink)
        ax.tick_params(axis='y',length=0,pad=4);ax.tick_params(axis='x',length=2.5,labelsize=7)
        for s in ['top','left','right']:ax.spines[s].set_visible(False)
        ax.spines['bottom'].set_color(grey)
    fig.text(.56,.393,'Held-out AUC difference (augmented minus baseline)',ha='center',fontsize=8.5,color=ink)
    fig.text(.028,.306,'b',fontweight='bold',fontsize=12,color=ink)
    fig.text(.070,.306,'External transfer to PUCH (n = 55)',fontsize=10.2,fontweight='bold',color=ink)
    ax=fig.add_axes([.073,.068,.89,.190]);ax.axis('off')
    ax.plot([0,1],[.74,.74],lw=.7,color='#79919C',transform=ax.transAxes)
    for x,t in zip([.0,.42,.69],['Fixed rank model','AUC','Brier score']):
        ax.text(x,.90,t,fontsize=8.3,fontweight='bold',transform=ax.transAxes,color=ink)
    for k,(r,name) in enumerate(zip(ext,['IFN-4','IFN-4 + APM (6 genes)','Paired change (95% CI)'])):
        y=.54-.29*k;fmt=lambda x,prec:f'{float(x):.{prec}f}'.replace('-','−')
        isdelta=k==2
        a=fmt(r['auc'],3); b=fmt(r['brier'],4)
        if isdelta:
            a+=f" ({fmt(r['auc_lo'],3)} to {fmt(r['auc_hi'],3)})"
            b+=f" ({fmt(r['brier_lo'],4)} to {fmt(r['brier_hi'],4)})"
        for x,s in zip([0,.42,.69],[name,a,b]):
            ax.text(x,y,s,fontsize=8.3,transform=ax.transAxes,color=ink,fontweight='bold' if isdelta and x==0 else 'normal')
    args.output.mkdir(parents=True,exist_ok=True)
    for extname in ['pdf','png']:fig.savefig(args.output/('Fig6.'+extname),dpi=450)
    plt.close(fig)
if __name__=='__main__':main()
