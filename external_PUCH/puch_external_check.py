#!/usr/bin/env python3
"""Separate revision-stage PUCH check; numpy/pandas only.

Download original authors' files into a private input directory, then generate
aggregate results only. Never copy the input directory into a public archive.
The four original benchmark cohorts remain unchanged; PUCH is not pooled.
"""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
import numpy as np
import pandas as pd

COMMIT = '275be3c340e4cba698d085bc1a875b92dd45602c'
REPO = 'https://github.com/xmuyulab/ims_gene_signature'
EXP_BLOB = '52dce854ae5c5cc869dace96335884eccbf18acf'
EXP_BYTES = 16844340
CLI_BLOB = '82e7c0f3cfd5a66fc9957e02fbec3ed60ad3e2c2'
GENE_SETS = {'four_gene':['CD7','CXCL9','LYZ','STAT1'],
             'IFN4':['IFNG','IDO1','CXCL10','HLA-DRA']}
SEED = 20261123
N_BOOT = 5000

def download(name, directory):
    path = directory / name
    if path.exists():
        return path
    if name == 'mel_puch_exp_data.csv':
        url = f'https://api.github.com/repos/xmuyulab/ims_gene_signature/git/blobs/{EXP_BLOB}'
        req = urllib.request.Request(url, headers={'Accept':'application/vnd.github.raw+json','User-Agent':'PUCH-reproducibility-check'})
    else:
        url = f'https://raw.githubusercontent.com/xmuyulab/ims_gene_signature/{COMMIT}/data/{name}'
        req = urllib.request.Request(url, headers={'User-Agent':'PUCH-reproducibility-check'})
    temporary = path.with_suffix(path.suffix + '.partial')
    with urllib.request.urlopen(req, timeout=600) as response, temporary.open('wb') as out:
        while True:
            block = response.read(1024 * 1024)
            if not block: break
            out.write(block)
    temporary.replace(path)
    return path

def digest(path):
    b = path.read_bytes()
    return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),
            'git_blob_sha1':hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--download',action='store_true')
    a = parser.parse_args()
    a.input_dir.mkdir(parents=True,exist_ok=True); a.output_dir.mkdir(parents=True,exist_ok=True)
    names = ['mel_puch_exp_data.csv','mel_puch_cli_data.csv']
    if a.download:
        for name in names: download(name,a.input_dir)
    hashes = {name:digest(a.input_dir/name) for name in names}
    assert hashes[names[0]]['bytes']==EXP_BYTES, 'Expression file incomplete'
    assert hashes[names[0]]['git_blob_sha1']==EXP_BLOB, 'Expression file differs from pinned author blob'
    assert hashes[names[1]]['git_blob_sha1']==CLI_BLOB, 'Clinical file differs from pinned author blob'
    expression = pd.read_csv(a.input_dir/names[0],index_col=0)
    clinical = pd.read_csv(a.input_dir/names[1],index_col=0)
    assert expression.columns.is_unique and clinical.index.is_unique
    assert set(clinical['treatment'])=={'pre'}
    assert set(clinical['data_set'])=={'PUCH'}
    assert len(clinical)==60 and len(expression.columns)==55
    ids = expression.columns  # Preserve deposited matrix order; do not publish IDs.
    assert set(ids).issubset(set(clinical.index))
    clinical = clinical.loc[ids]
    labels = clinical['response'].astype(int)
    assert labels.value_counts().to_dict()=={0:35,1:14,-1:6}, 'Not original 55-patient objective-response population'
    genes = sum(GENE_SETS.values(),[])
    assert all((expression.index==g).sum()==1 for g in genes), 'Missing or duplicated target gene'
    assert np.isfinite(expression.loc[genes].to_numpy()).all()
    lo = float(np.nanmin(expression.to_numpy())); hi = float(np.nanmax(expression.to_numpy()))
    assert 0 <= lo and hi < 500, 'Deposited scale differs from the authors\' log-scale matrix'
    # Authors describe RSEM FPKM followed by log2; their loader preserves this
    # matrix on its deposited scale. No extra log or housekeeping subtraction.
    # Paper housekeeping normalization concerns its own signature comparison.
    selected = expression.loc[genes].T.astype(float)
    assert (selected.std(axis=0,ddof=1)>0).all()
    standardized = (selected-selected.mean(axis=0))/selected.std(axis=0,ddof=1)
    scores = pd.DataFrame({name:standardized[gs].mean(axis=1) for name,gs in GENE_SETS.items()})
    y = (labels.to_numpy()==1).astype(int)  # CR/PR vs SD/PD; no durable-benefit recoding.
    positives = np.flatnonzero(y==1); negatives = np.flatnonzero(y==0)
    rng = np.random.default_rng(SEED)
    ir = rng.integers(0,len(positives),size=(N_BOOT,len(positives)))
    inn = rng.integers(0,len(negatives),size=(N_BOOT,len(negatives)))
    aucs,boots = {},{}
    for name in GENE_SETS:
        s = scores[name].to_numpy()
        difference = s[positives,None] - s[negatives][None,:]
        comparison = (difference>0).astype(float) + 0.5*(difference==0)
        boots[name] = comparison[ir[:,:,None],inn[:,None,:]].mean(axis=(1,2))
        aucs[name] = {'auc':float(comparison.mean()),
                      'ci95':[float(x) for x in np.quantile(boots[name],[.025,.975])]}
    delta = boots['four_gene']-boots['IFN4']
    corr = float(np.corrcoef(scores['four_gene'].rank(method='average'),scores['IFN4'].rank(method='average'))[0,1])
    result = {'status':'eligible_separate_revision_stage_external_check',
      'cohort':'PUCH; Cui et al., npj Genomic Medicine 2021, DOI 10.1038/s41525-021-00169-w',
      'population':{'clinical_file_n':60,'expression_file_n':55,'intersection_n':55,
                    'clinical_without_expression_n':5,'CR_PR_n':14,'SD_n':6,'PD_n':35,
                    'responders_n':14,'nonresponders_n':41},
      'expression':{'gene_rows':len(expression),'scale':'Author-deposited processed log-scale expression; no additional transformation',
                    'global_min':lo,'global_max':hi,'four_gene_coverage':4,'IFN4_coverage':4},
      'scoring':'Positive mean of per-gene patient z-scores (ddof=1), estimated across all 55 patients without outcome labels',
      'spearman_four_gene_IFN4':corr,'standalone_AUC':aucs,
      'paired_delta_AUC_four_gene_minus_IFN4':{'estimate':aucs['four_gene']['auc']-aucs['IFN4']['auc'],
                                            'ci95':[float(x) for x in np.quantile(delta,[.025,.975])]},
      'bootstrap':{'replicates':N_BOOT,'seed':SEED,'method':'Outcome-stratified sampling with replacement; identical patients/resamples for both fixed scores; percentile CI'},
      'not_performed':['No gene reselection','No new cut-point','No fitting of an augmented prediction model','No pooling with the original four cohorts'],
      'software':{'numpy':np.__version__,'pandas':pd.__version__},
      'source_commit':COMMIT,'source_repository':REPO,'input_hashes':hashes}
    (a.output_dir/'PUCH_aggregate_results.json').write_text(json.dumps(result,indent=2)+'\n')
    rows = [{'comparison':name,**v} for name,v in aucs.items()]
    out = [{'quantity':r['comparison']+'_AUC','estimate':r['auc'],'CI95_low':r['ci95'][0],'CI95_high':r['ci95'][1]} for r in rows]
    d=result['paired_delta_AUC_four_gene_minus_IFN4']
    out += [{'quantity':'delta_AUC_four_gene_minus_IFN4','estimate':d['estimate'],'CI95_low':d['ci95'][0],'CI95_high':d['ci95'][1]},
            {'quantity':'spearman_four_gene_IFN4','estimate':corr,'CI95_low':None,'CI95_high':None}]
    pd.DataFrame(out).to_csv(a.output_dir/'PUCH_aggregate_results.csv',index=False)
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
