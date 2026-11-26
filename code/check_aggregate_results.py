#!/usr/bin/env python3
"""Check all archived summaries using aggregate-only bootstrap and null arrays."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[1];checks=[]
def add(name,error):
 assert error<1e-10,(name,error);checks.append({'check':name,'max_abs_difference':float(error)})
d=pd.read_csv(B/'results/increment/increment_all_results.csv')
for row in d.itertuples():
 b=pd.read_csv(B/f'results/increment/bootstrap/{row.cohort}/{row.panel}_bootstrap.csv')
 assert len(b)==500
 for m in ['auc_M0','auc_M1','delta_auc','brier_M0','brier_M1','delta_brier']:
  q=b[m].quantile([.025,.975]).values;target=[getattr(row,m+'_low'),getattr(row,m+'_high')];add('increment '+row.cohort+' '+row.panel+' '+m,np.max(np.abs(q-target)))
new=d[d.panel!='Four-gene'];assert len(new)==20
assert ((new.delta_auc_low<0)&(new.delta_auc_high>0)).all() and ((new.delta_brier_low<0)&(new.delta_brier_high>0)).all()
s=pd.read_csv(B/'results/reference/random_reference_summary.csv')
for row in s.itertuples():
 v=pd.read_csv(B/f'results/reference/null_auc_{row.background}.csv.gz')[row.cohort].values;assert len(v)==10000
 for k,q in [('random_q025',.025),('random_q25',.25),('random_median',.5),('random_q75',.75),('random_q975',.975)]:add('null '+row.background+' '+row.cohort+' '+k,abs(np.quantile(v,q)-getattr(row,k)))
 # Percentiles and Monte Carlo counts are audited from the unrounded NPZ values.
 arr=np.load(B/'results/reference/random_reference_auc.npz')[row.background+'__'+row.cohort]
 add('null mean '+row.background+' '+row.cohort,abs(arr.mean()-row.random_mean))
score=pd.read_csv(B/'results/reference/scoring_auc.csv');boots=np.load(B/'results/reference/scoring_auc_bootstrap.npz')
for row in score.itertuples():
 vals=boots[row.cohort+'__'+row.panel+'__'+row.method];add('rank AUC CI '+row.cohort+' '+row.panel+' '+row.method,np.max(np.abs(np.quantile(vals,[.025,.975])-[row.ci95_low,row.ci95_high])))
paired=pd.read_csv(B/'results/reference/rank_vs_z_paired_auc.csv')
for row in paired.itertuples():
 vals=boots[row.cohort+'__'+row.panel+'__rank_common']-boots[row.cohort+'__'+row.panel+'__z_common'];add('paired rank CI '+row.cohort+' '+row.panel,np.max(np.abs(np.quantile(vals,[.025,.975])-[row.paired_ci95_low,row.paired_ci95_high])))
print(json.dumps({'passed':True,'checks':len(checks),'maximum_error':max(x['max_abs_difference'] for x in checks),'scope':'Aggregate arithmetic only; no patient data were needed.'},indent=2))
