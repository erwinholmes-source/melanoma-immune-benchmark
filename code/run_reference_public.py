#!/usr/bin/env python3
"""Post hoc v19 background/rank analysis; see analysis_plan.txt. No manuscript edits.

Portable invocation after installing requirements.txt:
  python run_reference.py --workspace /path/to/workspace --output /path/to/output
Supply prepared score-free v7 inputs and the original public v16 archive. No patient identifiers or scores are written unless --private-audit is explicitly supplied.
"""
import argparse, hashlib, json, platform, sys
from pathlib import Path
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--workspace', type=Path, required=True)
ap.add_argument('--v7-root', type=Path, required=True)
ap.add_argument('--v16-root', type=Path, required=True)
ap.add_argument('--private-audit', type=Path)
ap.add_argument('--output', type=Path, required=True)
ap.add_argument('--legacy-scipy-dir', type=Path)
args = ap.parse_args()
if args.legacy_scipy_dir:
    sys.path.append(str(args.legacy_scipy_dir))
import scipy
from scipy import stats

W=args.workspace.resolve(); O=args.output.resolve(); O.mkdir(parents=True, exist_ok=True)
R=args.v7_root.resolve(); V16=args.v16_root.resolve()
OLD=V16/'modules/Melanoma_immune_panel_Reproducibility_Addendum07/reference_results/addendum07'
A10=V16/'modules/Melanoma_immune_panel_Reproducibility_Addendum10/reference_results/addendum10'
if args.private_audit: args.private_audit.mkdir(parents=True,exist_ok=True)
IND=['GSE91061','GSE78220','Liu2019','Gide2019']; ALL=['GSE215868']+IND
ID={'GSE215868':'gsm','GSE91061':'Sample','GSE78220':'sample_id','Liu2019':'patient_id','Gide2019':'patient_id'}
SOURCES={
 'HGNC':('v7_work/mapping/hgnc_complete_set.txt','6f43d6ff43aa9fdfa5fb2f20a20a7cace66e6e02e2a0dcf19d9b726e2e248d20'),
 'GSE91061':('melanoma_revision_work/external_audit/riaz_rld.csv','ac69e754f36226a24c529fa8df6ae116fb0fbf30e365ed96757f9ee92d0dcec1'),
 'GSE78220':('v7_work/new_bulk/hugo2016/expression_gene_by_patient.csv','4b16bfa6f38340c06513dd22d38d5db8683eeac5fc4b3510510008f05f634a91'),
 'Liu2019':('v7_work/extra_bulk/liu2019/expression_baseline_log2TPM1.csv.gz','472567fe70f1cc36a77fb2696a2684fb7e6ee230aae1d402f51684f28dd63d9c'),
 'Gide2019':('v7_work/new_bulk/gide2019/expression_gene_by_patient.csv','eeb35954ce3b81452bb4c0689216941ff0218d8bd89e6d73a5498793cc011fcb')}
PANELS={
 'Four-gene':['CD7','CXCL9','LYZ','STAT1'],
 'IFN-4':['IFNG','IDO1','CXCL10','HLA-DRA'],
 'ICR-20':['IFNG','TBX21','CD8B','CD8A','IL12B','STAT1','IRF1','CXCL9','CXCL10','CCL5','GNLY','PRF1','GZMA','GZMB','GZMH','CD274','CTLA4','FOXP3','IDO1','PDCD1'],
 'GEP-18':['CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6','HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT'],
 'Chemokine-12':['CCL2','CCL3','CCL4','CCL5','CCL8','CCL18','CCL19','CCL21','CXCL9','CXCL10','CXCL11','CXCL13'],
 'TLS-9':['CD79B','CD1D','CCR6','LAT','SKAP1','CETP','EIF1AY','RBP5','PTGDS'],
 'APM-7':['B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','TAPBP']}
NS=10000; NB=5000
manifest=[]
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def register(p,role,expected=None):
    h=sha(p)
    if expected is not None: assert h==expected,(str(p),h)
    label=p.name
    for prefix,base in [('prepared-v7',R),('v16-public',V16),('source-workspace',W)]:
        try: label=prefix+'/'+str(p.relative_to(base)); break
        except ValueError: pass
    manifest.append({'path':label,'role':role,'bytes':p.stat().st_size,'sha256':h})
    return p
def zrow(m):return m.sub(m.mean(axis=1),axis=0).div(m.std(axis=1,ddof=1),axis=0)
def auc_rows(v,y):
    v=np.asarray(v,float); y=np.asarray(y,int); n1=int(y.sum()); n0=len(y)-n1
    ranks=stats.rankdata(v,axis=1,method='average')
    return (ranks[:,y==1].sum(axis=1)-n1*(n1+1)/2)/(n1*n0)
def bootstrap(v,y,wp,wn):
    a=v[y==1]; b=v[y==0]
    q=(a[:,None]>b[None,:])+0.5*(a[:,None]==b[None,:])
    return np.einsum('bi,ij,bj->b',wp,q,wn,optimize=True)/(len(a)*len(b))
def rho(a,b):return float(stats.spearmanr(a,b).statistic)
def dump(obj,path):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')

derived_source_comparison=[]
for k,(rel,h) in SOURCES.items():
    register(W/rel,'full_expression' if k!='HGNC' else 'gene_annotation',h if k=='HGNC' else None)
    derived_source_comparison.append({'source':k,'relative_path':rel,'historical_sha256':h,'reconstructed_sha256':sha(W/rel),'byte_identical_to_historical':sha(W/rel)==h})
if (W/'v7_work/mapping/mapping_manifest.json').exists(): register(W/'v7_work/mapping/mapping_manifest.json','HGNC_download_provenance')
hg=pd.read_csv(W/SOURCES['HGNC'][0],sep='\t',dtype=str)
approved=hg[hg.status.eq('Approved')].copy()
u=approved[approved.entrez_id.notna()]
u=u[~u.symbol.duplicated(keep=False)&~u.entrez_id.duplicated(keep=False)]
ez2sym=dict(zip(u.entrez_id.astype(int),u.symbol))
pc=set(approved.loc[approved.locus_group.eq('protein-coding gene'),'symbol'])
assert pc and len(pc)>15000
clin={}; compact={}; cz={}; full={}; fz={}; audits=[]; sourcechecks=[]
for c in ALL:
    p=register(R/f'inputs/bulk/{c}/clinical.csv','score_free_clinical')
    a=pd.read_csv(p); a.index=a[ID[c]].astype(str); assert a.index.is_unique
    a.index.name=None;clin[c]=a
    p=register(R/f'inputs/rank/{c}/source_symbol_expression.csv.gz','historical_background_expression')
    m=pd.read_csv(p,index_col=0);m.columns=m.columns.astype(str);m=m.loc[:,a.index].astype(float)
    compact[c]=m;cz[c]=zrow(m)
    if c not in IND:continue
    raw=pd.read_csv(W/SOURCES[c][0],index_col=0);raw.columns=raw.columns.astype(str)
    assert set(a.index)<=set(raw.columns)
    raw=raw.loc[:,a.index].astype(float)
    nraw=len(raw)
    if c=='GSE91061':
        raw.index=raw.index.astype(int)
        raw=raw.loc[[e for e in raw.index if e in ez2sym]].copy()
        raw.index=[ez2sym[e] for e in raw.index]
    else:raw.index=raw.index.astype(str)
    duplicates=set(raw.index[raw.index.duplicated(keep=False)])
    m=raw.loc[~raw.index.duplicated(keep=False)].copy()
    assert m.index.is_unique
    compare=[g for g in compact[c].index.unique() if g in m.index and (compact[c].index==g).sum()==1]
    d=float(np.nanmax(np.abs(m.loc[compare].values-compact[c].loc[compare].values)))
    assert d<1e-9,(c,d)
    sourcechecks.append({'cohort':c,'test':'full_vs_compact_expression','max_abs_difference':d,'n_compared_genes':len(compare),'passed':True})
    full[c]=m;fz[c]=zrow(m)
    finite=np.isfinite(m.values).all(axis=1); varying=m.std(axis=1,ddof=1).gt(0).values
    eligible=set(m.index[finite&varying])&pc
    audits.append({'cohort':c,'raw_rows':nraw,'unique_mapped_rows':len(m),'duplicate_symbols_removed':len(duplicates),'baseline_patients':len(a),'response_evaluable':int(a.response.notna().sum()),'responders':int(a.response.fillna(0).sum()),'finite_nonconstant_approved_protein_coding':len(eligible)})
    print(f'INPUT {c}: {len(m)} unique mapped genes, {len(eligible)} eligible PC, {len(a)} baseline',flush=True)

# The broad universe is defined without clinical outcomes.
broad=sorted(set.intersection(*[set(full[c].index[np.isfinite(full[c].values).all(axis=1)&full[c].std(axis=1,ddof=1).gt(0).values])&pc for c in IND]))
hist=sorted(set.intersection(*[set(cz[c].index[cz[c].notna().all(axis=1)]) for c in ALL]))
assert len(hist)==733,len(hist)
infl=[g for g in hist if all(rho(cz[c].loc[g],cz[c].loc['CD8A'])>=0.5 for c in ALL)]
assert len(infl)==142,len(infl)
assert set(PANELS['Four-gene']+PANELS['IFN-4'])<=set(broad)
universes={'broad_protein_coding':broad,'historical733':hist,'inflamed142':infl}
for name,genes in universes.items():pd.DataFrame({'gene':genes}).to_csv(O/f'universe_{name}.csv',index=False)
allgenes=sorted(set.union(*[set(m.index) for m in full.values()]))
broadset=set(broad);histset=set(hist);inflset=set(infl)
elig=[]
for g in allgenes:
    row={'gene':g,'hgnc_approved_protein_coding':g in pc,'in_broad':g in broadset,'in_historical733':g in histset,'in_inflamed142':g in inflset}
    for c in IND:
        row[c]=('not_uniquely_measured' if g not in full[c].index else 'nonfinite' if not np.isfinite(full[c].loc[g]).all() else 'constant' if full[c].loc[g].std(ddof=1)==0 else 'eligible_expression')
    elig.append(row)
pd.DataFrame(elig).to_csv(O/'gene_eligibility_audit.csv.gz',index=False,compression={'method':'gzip','mtime':0})
pd.DataFrame(audits).to_csv(O/'input_dimensions.csv',index=False)
pd.DataFrame(sourcechecks).to_csv(O/'input_validation.csv',index=False)
print('UNIVERSES', {k:len(v) for k,v in universes.items()},flush=True)

# Identical gene-set draws are scored in all four evaluation cohorts.
ys={c:clin[c].loc[clin[c].response.notna(),'response'].astype(int).values for c in IND}
panelauc={c:float(auc_rows(fz[c].loc[PANELS['Four-gene'],clin[c].response.notna()].mean(axis=0).values[None,:],ys[c])[0]) for c in IND}
nulls={}; summaries=[]; reproduce=[]
oldpath=register(OLD/'a07_random_gene_sets.npz','validation_target_only')
old=np.load(oldpath,allow_pickle=False)
for name,seed,oldname in [('broad_protein_coding',20261123,None),('historical733',20261102,'all_genes'),('inflamed142',20261103,'inflamed_genes')]:
    genes=universes[name];rng=np.random.default_rng(seed)
    draws=np.array([rng.choice(len(genes),4,replace=False) for _ in range(NS)])
    drawframe=pd.DataFrame(np.array(genes)[draws],columns=['gene1','gene2','gene3','gene4']);drawframe.insert(0,'draw',np.arange(1,NS+1))
    drawframe.to_csv(O/f'draws_{name}.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    vals={}
    for c in IND:
        base=fz[c] if oldname is None else cz[c]
        mat=base.loc[genes,clin[c].response.notna()].values
        assert np.isfinite(mat).all()
        vals[c]=auc_rows(mat[draws].mean(axis=1),ys[c])
        if oldname:
            delta=float(np.max(np.abs(vals[c]-old[f'{oldname}__{c}'])))
            assert delta<1e-12,(name,c,delta)
            reproduce.append({'background':name,'cohort':c,'draws_compared':NS,'max_abs_auc_difference':delta,'passed':True})
    vals['mean_independent4']=np.mean([vals[c] for c in IND],axis=0)
    if oldname:
        delta=float(np.max(np.abs(vals['mean_independent4']-old[f'{oldname}__mean_independent4'])))
        assert delta<1e-12
        reproduce.append({'background':name,'cohort':'mean_independent4','draws_compared':NS,'max_abs_auc_difference':delta,'passed':True})
    for c,v in vals.items():
        obs=panelauc[c] if c in IND else float(np.mean(list(panelauc.values())))
        summaries.append({'background':name,'cohort':c,'genes':len(genes),'draws':NS,'seed':seed,'panel_auc':obs,'random_mean':float(v.mean()),'random_median':float(np.median(v)),'random_q025':float(np.quantile(v,.025)),'random_q25':float(np.quantile(v,.25)),'random_q75':float(np.quantile(v,.75)),'random_q975':float(np.quantile(v,.975)),'panel_percentile':100*float(np.mean(v<obs)+.5*np.mean(v==obs)),'p_upper_nominal':float((1+(v>=obs).sum())/(NS+1))})
        nulls[f'{name}__{c}']=v
    vf=pd.DataFrame(vals);vf.insert(0,'draw',np.arange(1,NS+1));vf.to_csv(O/f'null_auc_{name}.csv.gz',index=False,compression={'method':'gzip','mtime':0})
pd.DataFrame(summaries).to_csv(O/'random_reference_summary.csv',index=False)
pd.DataFrame(reproduce).to_csv(O/'historical_reproduction.csv',index=False)
np.savez_compressed(O/'random_reference_auc.npz',**nulls)

# Frozen common membership: gene-content-matched z and within-sample rank scores.
coverage=[];common={};scoreframes={};rows=[];corr=[];delta_rows=[];bootstore={}
for p,gs in PANELS.items():
    common[p]=[g for g in gs if g in broad]
    assert common[p],p
    for c in IND:
        coverage.append({'panel':p,'cohort':c,'published_n':len(gs),'available_n':len([g for g in gs if g in full[c].index]),'common_n':len(common[p]),'available_genes':';'.join(g for g in gs if g in full[c].index),'absent_from_cohort':';'.join(g for g in gs if g not in full[c].index),'fixed_common_genes':';'.join(common[p]),'excluded_from_fixed_common':';'.join(g for g in gs if g not in common[p]),'negative_genes':'CD276' if p=='GEP-18' else ''})
pd.DataFrame(coverage).to_csv(O/'panel_coverage.csv',index=False)
for ci,c in enumerate(IND):
    a=clin[c];m=a.response.notna().values;y=ys[c]
    ranks=stats.rankdata(full[c].loc[broad].values,axis=0,method='average')
    ranks=pd.DataFrame((ranks-1)/(len(broad)-1)-.5,index=broad,columns=a.index)
    sf=pd.DataFrame({'patient_id':a.index,'response':a.response.values})
    scores={};boots={}
    rng=np.random.default_rng(20261200+ci);n1=int(y.sum());n0=len(y)-n1
    wp=rng.multinomial(n1,np.full(n1,1/n1),size=NB);wn=rng.multinomial(n0,np.full(n0,1/n0),size=NB)
    for p,gs in PANELS.items():
        for method,avail in [('z_available',[g for g in gs if g in full[c].index]),('z_common',common[p]),('rank_common',common[p])]:
            signs=np.array([-1. if p=='GEP-18' and g=='CD276' else 1. for g in avail])
            source=ranks if method=='rank_common' else fz[c]
            assert np.isfinite(source.loc[avail].values).all(),(c,p,method)
            s=(source.loc[avail].values*signs[:,None]).mean(axis=0)
            key=f'{p}__{method}';scores[key]=s;sf[key]=s
            auc=float(auc_rows(s[m][None,:],y)[0]);bs=bootstrap(s[m],y,wp,wn);boots[key]=bs
            rows.append({'cohort':c,'panel':p,'method':method,'genes_used':len(avail),'n':len(y),'responders':n1,'auc':auc,'ci95_low':float(np.quantile(bs,.025)),'ci95_high':float(np.quantile(bs,.975)),'bootstrap_draws':NB,'bootstrap_seed':20261200+ci})
    for p in PANELS:
        for method in ['z_available','z_common','rank_common']:
            v=scores[f'{p}__{method}']
            corr.append({'cohort':c,'panel':p,'method':method,'n_baseline':len(a),'rho_with_IFN4_same_method':rho(v,scores[f'IFN-4__{method}']),'rho_with_Four_gene_same_method':rho(v,scores[f'Four-gene__{method}'])})
        bs=boots[f'{p}__rank_common']-boots[f'{p}__z_common']
        da=float(auc_rows(scores[f'{p}__rank_common'][m][None,:],y)[0]-auc_rows(scores[f'{p}__z_common'][m][None,:],y)[0])
        delta_rows.append({'cohort':c,'panel':p,'n':len(y),'rank_minus_z_common_auc':da,'paired_ci95_low':float(np.quantile(bs,.025)),'paired_ci95_high':float(np.quantile(bs,.975)),'rho_rank_vs_z_common':rho(scores[f'{p}__rank_common'],scores[f'{p}__z_common'])})
    if args.private_audit: sf.to_csv(args.private_audit/f'patient_scores_{c}.csv',index=False)
    for k,b in boots.items():bootstore[f'{c}__{k}']=b
pd.DataFrame(rows).to_csv(O/'scoring_auc.csv',index=False)
pd.DataFrame(corr).to_csv(O/'scoring_correlations.csv',index=False)
pd.DataFrame(delta_rows).to_csv(O/'rank_vs_z_paired_auc.csv',index=False)
np.savez_compressed(O/'scoring_auc_bootstrap.npz',**bootstore)

metrics=pd.DataFrame(rows)
aggregate=metrics.groupby(['panel','method'],sort=False).agg(mean_cohort_auc=('auc','mean'),min_cohort_auc=('auc','min'),max_cohort_auc=('auc','max')).reset_index()
aggregate.to_csv(O/'scoring_mean_cohort_auc.csv',index=False)
metadata={'analysis':'v19 post hoc three-tier reference and single-sample-rank sensitivity','python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'command_template':'python run_reference_public.py --workspace WORKSPACE --v7-root PREPARED_V7 --v16-root V16_PUBLIC --output OUTPUT','hgnc_download':{'url':'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt','retrieved':'2026-09-07T03:32:49.738763+00:00','sha256':SOURCES['HGNC'][1]},'universe_sizes':{k:len(v) for k,v in universes.items()},'common_panel_membership':common,'source_validation_passed':True,'historical_random_AUC_reproduction_passed':True,'historical_background_is_conditional_on_five_cohort_expression':True,'new_broad_background_uses_four_cohort_expression_without_outcomes':True,'rank_method_note':'Mean signed within-sample centered percentiles over one frozen common protein-coding universe; source transformed values retained; not the original weighted/algorithmic signature implementations.'}
# Aggregate original AUC targets, not original patient-level predictions or scores.
old07=pd.read_csv(register(OLD/'a07_response_metrics.csv','aggregate_validation_target_only'))
old10=pd.read_csv(register(A10/'a10_response_metrics.csv','aggregate_validation_target_only'))
aggregate_checks=[]
for row in metrics[metrics.method.eq('z_available')].itertuples():
    target=old07 if row.panel in ['Four-gene','IFN-4'] else old10
    hit=target[(target.cohort==row.cohort)&target.measure.isin([row.panel,row.panel+' primary'])]
    assert len(hit)==1,(row.cohort,row.panel)
    d=abs(row.auc-float(hit.iloc[0].auc));assert d<1e-12,(row.cohort,row.panel,d)
    aggregate_checks.append({'cohort':row.cohort,'panel':row.panel,'max_abs_auc_difference':d,'passed':True})
pd.DataFrame(aggregate_checks).to_csv(O/'aggregate_AUC_reproduction.csv',index=False)
pd.DataFrame(derived_source_comparison).to_csv(O/'reconstructed_source_byte_comparison.csv',index=False)
metadata['derived_input_policy']='The archived HGNC snapshot is hash-locked. Derived expression files are recorded by SHA-256 but may differ bytewise after reconstruction; numerical concordance with compact inputs, full historical null AUC arrays and archived aggregate AUCs is required.'
metadata['patient_level_outputs_written']=bool(args.private_audit)
dump(metadata,O/'run_info.json');dump(manifest,O/'input_manifest.json')
print(pd.DataFrame(summaries).query("cohort=='mean_independent4'").to_string(index=False),flush=True)
print(aggregate.to_string(index=False),flush=True)
print('SUCCESS',flush=True)
