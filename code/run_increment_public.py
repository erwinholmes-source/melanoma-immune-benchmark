#!/usr/bin/env python3
"""Rebuild and run Table 4 incremental comparisons; outputs aggregate results only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import a11_core_original as core
PKG=Path(__file__).resolve().parents[1]
COHORTS=['GSE91061','GSE78220','Liu2019','Gide2019']
IFN4=['IFNG','IDO1','CXCL10','HLA-DRA']
PANELS={
'Four-gene':['CD7','CXCL9','LYZ','STAT1'],
'ICR-20':['IFNG','TBX21','CD8B','CD8A','IL12B','STAT1','IRF1','CXCL9','CXCL10','CCL5','GNLY','PRF1','GZMA','GZMB','GZMH','CD274','CTLA4','FOXP3','IDO1','PDCD1'],
'GEP-18':['CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6','HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT'],
'Chemokine-12':['CCL2','CCL3','CCL4','CCL5','CCL8','CCL18','CCL19','CCL21','CXCL9','CXCL10','CXCL11','CXCL13'],
'TLS-9':['CD79B','CD1D','CCR6','LAT','SKAP1','CETP','EIF1AY','RBP5','PTGDS'],
'APM-7':['B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','TAPBP']}
PATHS={'GSE91061':'melanoma_revision_work/external_audit/riaz_rld.csv','GSE78220':'v7_work/new_bulk/hugo2016/expression_gene_by_patient.csv','Liu2019':'v7_work/extra_bulk/liu2019/expression_baseline_log2TPM1.csv.gz','Gide2019':'v7_work/new_bulk/gide2019/expression_gene_by_patient.csv'}
METRICS=['auc_M0','auc_M1','delta_auc','brier_M0','brier_M1','delta_brier']
def prep_factory(signs):
 def prep(train,test,w,clinical=False):
  assert not clinical
  gz,gt,mu,sd,const=core.weighted_standardize(train,test,w)
  s=np.c_[gz[:,:4].mean(axis=1),(gz[:,4:]*signs).mean(axis=1)];st=np.c_[gt[:,:4].mean(axis=1),(gt[:,4:]*signs).mean(axis=1)]
  a,b,smu,ssd,sconst=core.weighted_standardize(s,st,w)
  return a,b,dict(gene_mean=mu.tolist(),gene_sd=sd.tolist(),score_mean=smu.tolist(),score_sd=ssd.tolist(),clinical_mean=[],clinical_sd=[],constant_genes=np.flatnonzero(const).tolist(),constant_scores=np.flatnonzero(sconst).tolist(),constant_clinical=[])
 return prep

def prepare(W):
 I=W/'increment_inputs';I.mkdir(exist_ok=True)
 h=pd.read_csv(PKG/'resources/hgnc_complete_set.txt',sep='\t',dtype=str);h=h[h.status.eq('Approved')&h.entrez_id.notna()];h=h[~h.symbol.duplicated(keep=False)&~h.entrez_id.duplicated(keep=False)];mapping=dict(zip(h.entrez_id.astype(int),h.symbol))
 genes=sorted(set(IFN4)|{g for gl in PANELS.values() for g in gl})
 for c in COHORTS:
  a=pd.read_csv(W/f'prepared_v7/inputs/bulk/{c}/clinical.csv',index_col=0);a.index=a.index.astype(str);m=pd.read_csv(W/PATHS[c],index_col=0);m.columns=m.columns.astype(str);m=m.loc[:,a.index]
  if c=='GSE91061':
   m.index=m.index.astype(int);m=m.loc[[i for i in m.index if i in mapping]].copy();m.index=[mapping[i] for i in m.index]
  assert m.index.is_unique and a.index.is_unique
  # Preserve the historical A10 compact export precision before fold-wise scaling.
  m.loc[[g for g in genes if g in m.index]].to_csv(I/f'{c}_genes.csv',float_format='%.12g')
  a.index.name='patient_id';a[['response']].to_csv(I/f'{c}_clinical.csv')
 return I

def run_one(arg):
 c,I,O=arg;I=Path(I);O=Path(O);idx=COHORTS.index(c);out=O/c;out.mkdir(parents=True,exist_ok=True)
 a=pd.read_csv(I/f'{c}_clinical.csv');a=a[a.response.notna()];ids=a.patient_id.astype(str).to_numpy();y=a.response.to_numpy(int);mat=pd.read_csv(I/f'{c}_genes.csv',index_col=0)
 rng=np.random.default_rng(np.random.SeedSequence([20261105,11,idx,1]));weights=[]
 for draw in range(500):
  w=np.zeros(len(y),int)
  for lab in [0,1]:
   ii=np.flatnonzero(y==lab);w+=np.bincount(rng.choice(ii,len(ii),replace=True),minlength=len(y))
  weights.append(w)
 rows=[]
 for pn,gl0 in PANELS.items():
  gl=[g for g in gl0 if g in mat.index];signs=np.array([-1 if pn=='GEP-18' and g=='CD276' else 1 for g in gl]);core.preprocess=prep_factory(signs);x=mat.loc[IFN4+gl,ids].T.to_numpy(float)
  primary,*_=core.paired_cv(x,y,np.ones(len(y)),np.random.SeedSequence([20261105,11,idx,0]))
  br=[]
  for draw,w in enumerate(weights):
   metrics,*_=core.paired_cv(x,y,w,np.random.SeedSequence([20261105,11,idx,2,draw]));br.append(dict(draw=draw,**metrics))
  b=pd.DataFrame(br);b.to_csv(out/f'{pn}_bootstrap.csv',index=False)
  s=dict(cohort=c,panel=pn,n=len(y),responders=int(y.sum()),n_genes=len(gl),missing_genes=';'.join(g for g in gl0 if g not in gl),bootstrap_requested=500,bootstrap_valid=len(b),bootstrap_reduced_folds=int((b.nfold<5).sum()),bootstrap_fallbacks=int(b.fallbacks.sum()),bootstrap_constant_folds=int(b.constant_folds.sum()),**primary)
  for key in METRICS:s[key+'_low'],s[key+'_high']=map(float,b[key].quantile([.025,.975]))
  rows.append(s)
 pd.DataFrame(rows).to_csv(out/'summary.csv',index=False);print(c,'finished',flush=True);return rows

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);I=prepare(a.workspace)
 with ProcessPoolExecutor(max_workers=a.workers) as pool:r=list(pool.map(run_one,[(c,str(I),str(a.output)) for c in COHORTS]))
 d=pd.DataFrame([x for rr in r for x in rr]);d.to_csv(a.output/'increment_all_results.csv',index=False);d[d.panel!='Four-gene'].to_csv(a.output/'table4_incremental_prediction.csv',index=False)
 reference=pd.read_csv(PKG/'results/increment/increment_all_results.csv').set_index(['cohort','panel']);new=d.set_index(['cohort','panel']);cols=METRICS+[k+s for k in METRICS for s in ['_low','_high']];error=float(np.max(np.abs(new.loc[reference.index,cols].values-reference[cols].values)))
 assert error<1e-8,error
 (a.output/'reproduction_receipt.json').write_text(json.dumps({'passed':True,'comparisons':len(d),'bootstrap_draws_per_comparison':500,'max_abs_aggregate_difference':error,'numpy':np.__version__,'pandas':pd.__version__,'scope':'All 24 point estimates and 95% interval endpoints rebuilt from prepared source inputs and regenerated patient bootstrap multiplicities. No private precomputed weights, predictions or patient-score tables were required.'},indent=2))
 print('All published aggregate results reproduced',error)
