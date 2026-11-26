#!/usr/bin/env python3
"""Rebuild five bulk input populations from published processed source files.
No patient-level source files are distributed in this package. Files are read from
--raw-dir and written only to the caller-selected private --workspace.
"""
import argparse,csv,gzip,hashlib,io,json,re,shutil
from pathlib import Path
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--raw-dir',type=Path,required=True);p.add_argument('--workspace',type=Path,required=True);a=p.parse_args()
W=a.workspace;W.mkdir(parents=True,exist_ok=True);R=a.raw_dir;V=W/'prepared_v7'
for rec in json.loads((P/'resources/source_files.json').read_text()):
 f=R/rec['file'];assert f.is_file(),f'Missing source: {rec["file"]}'
 assert hashlib.sha256(f.read_bytes()).hexdigest()==rec['sha256'],f'Source hash mismatch: {rec["file"]}'
MAP={'CR':1,'PR':1,'SD':0,'PD':0};full={};clinical={};audits=[]
def write(frame,path,**kw):
 path.parent.mkdir(parents=True,exist_ok=True);frame.to_csv(path,**kw)
def cp(src,path):path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,path)
for name in ['hgnc_complete_set.txt','mapping_manifest.json','fixed_rank_universe.csv']:cp(P/'resources'/name,W/'v7_work/mapping'/name)
hg=pd.read_csv(P/'resources/hgnc_complete_set.txt',sep='\t',dtype=str);h=hg[hg.status.eq('Approved')]
u=pd.read_csv(P/'resources/fixed_rank_universe.csv')
aliases={g:{g} for g in h.symbol}
for row in h.itertuples():
 for field in ['prev_symbol','alias_symbol']:
  val=getattr(row,field)
  if isinstance(val,str):
   for alias in val.split('|'):aliases.setdefault(alias,set()).add(row.symbol)
required=set(u.symbol);relevant=required|{k for k,v in aliases.items() if v&required}
# Yale: read SOFT in deposited order and corroborate all expression values with matrix.
text=gzip.open(R/'GSE215868_family.soft.gz','rt',encoding='utf8').read();expressions={};rows=[]
for block in re.split(r'(?m)^\^SAMPLE = ',text)[1:]:
 identifier=block.splitlines()[0].strip();row={'gsm':identifier}
 for line in block.splitlines():
  if line.startswith('!Sample_characteristics_ch1 = '):
   key,value=line.split(' = ',1)[1].split(': ',1);row[key]=value
 table=block.split('!sample_table_begin\n')[1].split('!sample_table_end')[0]
 df=pd.read_csv(io.StringIO(table),sep='\t',dtype={'ID_REF':str});assert not df.ID_REF.duplicated().any()
 expressions[identifier]=df.set_index('ID_REF').VALUE;rows.append(row)
yale=pd.DataFrame(expressions);cy=pd.DataFrame(rows).set_index('gsm');cy['response']=cy.bor.map(MAP)
text=gzip.open(R/'GSE215868_series_matrix.txt.gz','rt',encoding='utf8').read();tab=text.split('!series_matrix_table_begin\n')[1].split('!series_matrix_table_end')[0]
assert np.array_equal(pd.read_csv(io.StringIO(tab),sep='\t',index_col=0).loc[yale.index,yale.columns].values,yale.values)
full['GSE215868']=yale;clinical['GSE215868']=cy[['response']]
# Riaz: retain author transformed scale and original Pre sample-table order.
cr=pd.read_csv(R/'riaz_sample_table.csv');cr=cr[cr.PreOn=='Pre'].set_index('Sample');assert cr.PatientID.is_unique
cr['response']=cr.BOR.map(MAP);er=pd.read_csv(R/'riaz_rld.csv',index_col=0)
clinical['GSE91061']=cr[['response']];rr=er.loc[u.entrez_id,cr.index].copy();rr.index=u.symbol;full['GSE91061']=rr
cp(R/'riaz_rld.csv',W/'melanoma_revision_work/external_audit/riaz_rld.csv')
# Hugo: source metadata determines pre-treatment eligibility and repeated-lesion aggregation.
meta={}
for line in (R/'GSE78220_series_matrix.txt').read_text().splitlines():
 if line.startswith('!Sample_'):
  key,*values=next(csv.reader([line],delimiter='\t'));meta.setdefault(key,[]).append(values)
rows=[]
for i,title in enumerate(meta['!Sample_title'][0]):
 d={'source_sample_id':title}
 for values in meta['!Sample_characteristics_ch1']:
  if ': ' in values[i]:key,val=values[i].split(': ',1);d[key]=val
 rows.append(d)
ch=pd.DataFrame(rows);eh=pd.read_excel(R/'GSE78220_PatientFPKM.xlsx',index_col=0);eh.index=eh.index.astype(str).str.strip();eh.columns=eh.columns.str.split('.').str[0]
assert set(eh.columns)==set(ch.source_sample_id) and eh.index.is_unique
log=np.log2(eh+1);out={};clin=[];br={'Complete Response':'CR','Partial Response':'PR','Stable Disease':'SD','Progressive Disease':'PD'}
for source_patient,g in ch[ch['biopsy time'].eq('pre-treatment')].groupby('patient id',sort=False):
 assert g['anti-pd-1 response'].nunique()==1
 pid='Hugo2016_'+source_patient;out[pid]=log[g.source_sample_id].mean(axis=1);clin.append({'sample_id':pid,'response':MAP.get(br.get(g.iloc[0]['anti-pd-1 response']))})
eh=pd.DataFrame(out);eh.index.name='gene_symbol';ch=pd.DataFrame(clin).set_index('sample_id')
write(eh,W/'v7_work/new_bulk/hugo2016/expression_gene_by_patient.csv');full['GSE78220']=eh;clinical['GSE78220']=ch
# Liu: expression-file patient order, pre-PD1 context and biopsy on/before treatment.
el=pd.read_csv(R/'RNA_TPM_matrix.txt',sep='\t',index_col=0);cl=pd.read_csv(R/'Supplemental_Table_1_wAge.tsv',sep='\t').set_index('Sample ID').reindex(el.index)
assert el.index.is_unique and el.columns.is_unique
keep=cl['biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)'].eq(3)&pd.to_numeric(cl.daysBiopsyToPD1).le(0)
cl=cl.loc[keep].copy();cl.index.name='patient_id';cl['response']=cl.BR.str.strip().str.upper().map(MAP)
el=np.log2(el.loc[cl.index].T+1);el.index.name='gene'
write(el,W/'v7_work/extra_bulk/liu2019/expression_baseline_log2TPM1.csv.gz',compression={'method':'gzip','mtime':0});full['Liu2019']=el;clinical['Liu2019']=cl[['response']]
# Gide: processed author reanalysis, exact RECIST categories (durable SD remains SD).
import rdata
g=rdata.read_rda(R/'Gide_data.RData')['Gide_data'];cg=g['Samples'].copy();cg.columns=cg.columns.astype(str);xx=g['TPM']
eg=pd.DataFrame(np.asarray(xx),index=[str(x).strip() for x in xx.coords['dim_0'].values],columns=[str(x).strip() for x in xx.coords['dim_1'].values])
assert eg.index.is_unique and eg.columns.is_unique
cg=cg[cg.PREEDT.eq('PRE')&cg.Treatment.isin(['PD1','ipiPD1'])].copy();cg['patient_id']='Gide2019_'+cg.Treatment.astype(str)+'_'+cg['Patient.no.'].astype(str);cg['response']=cg['Best.RECIST.response'].map(MAP)
eg=np.log2(eg[cg.Sample]+1);eg.columns=cg.patient_id;eg.index.name='gene_symbol'
write(eg,W/'v7_work/new_bulk/gide2019/expression_gene_by_patient.csv',float_format='%.10g')
# Re-read the historical precision used as the downstream source.
eg=pd.read_csv(W/'v7_work/new_bulk/gide2019/expression_gene_by_patient.csv',index_col=0)
clinical['Gide2019']=cg.set_index('patient_id')[['response']];full['Gide2019']=eg
expected={'GSE215868':(105,104,45),'GSE91061':(51,49,10),'GSE78220':(26,26,14),'Liu2019':(104,102,41),'Gide2019':(72,72,39)}
for c,e in full.items():
 cl=clinical[c];assert e.columns.tolist()==cl.index.tolist();assert e.index.is_unique and np.isfinite(e.to_numpy(float)).all()
 got=(len(cl),int(cl.response.notna().sum()),int(cl.response.sum()));assert got==expected[c],(c,got)
 write(cl,V/f'inputs/bulk/{c}/clinical.csv')
 compact=e.loc[e.index.isin(relevant)].copy()
 write(compact,V/f'inputs/rank/{c}/source_symbol_expression.csv.gz',compression={'method':'gzip','mtime':0})
 audits.append({'cohort':c,'baseline_n':got[0],'response_evaluable_n':got[1],'responders':got[2],'compact_gene_rows':len(compact)})
cp(P/'resources/fixed_rank_universe.csv',V/'mapping/fixed_rank_universe.csv')
write(pd.DataFrame(audits),W/'source_preparation_counts.csv',index=False)
(W/'source_preparation_receipt.json').write_text(json.dumps({'all_raw_hashes_matched':True,'sources':json.loads((P/'resources/source_files.json').read_text()),'counts':audits,'scope':'Processed source matrices and response annotations, not alignment of sequencing reads. Source files were supplied locally; no network download was performed by this script.'},indent=2))
print(pd.DataFrame(audits).to_string(index=False));print('PREPARED')
