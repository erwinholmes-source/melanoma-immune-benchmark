"""Fixed paired ridge CV for Addendum 11; no feature/tuning search."""
import json
import numpy as np

GENES=['IFNG','IDO1','CXCL10','HLA-DRA','CD7','CXCL9','LYZ','STAT1']
CLINICAL=['LDH_elevated','M1c','priorCTLA4']

def expit(z): return np.exp(-np.logaddexp(0.,-np.asarray(z)))

def weighted_standardize(train,test,w):
    mu=np.average(train,axis=0,weights=w)
    var=np.sum(w[:,None]*(train-mu)**2,axis=0)/(w.sum()-1)
    const=var<1e-20; sd=np.sqrt(np.maximum(var,1e-20));sd[const]=1.
    a=(train-mu)/sd;b=(test-mu)/sd
    a[:,const]=0.;b[:,const]=0.
    return a,b,mu,sd,const

def preprocess(train,test,w,clinical=False):
    gz,gtest,gmu,gsd,gconst=weighted_standardize(train[:,:8],test[:,:8],w)
    scores=np.c_[gz[:,:4].mean(axis=1),gz[:,4:8].mean(axis=1)]
    tscores=np.c_[gtest[:,:4].mean(axis=1),gtest[:,4:8].mean(axis=1)]
    sz,st,smu,ssd,sconst=weighted_standardize(scores,tscores,w)
    if clinical:
        cz,ct,cmu,csd,cconst=weighted_standardize(train[:,8:11],test[:,8:11],w)
        a=np.c_[cz,sz];b=np.c_[ct,st]
    else:
        a,b=sz,st;cmu=csd=np.array([]);cconst=np.array([],bool)
    meta=dict(gene_mean=gmu.tolist(),gene_sd=gsd.tolist(),score_mean=smu.tolist(),score_sd=ssd.tolist(),
              clinical_mean=cmu.tolist(),clinical_sd=csd.tolist(),constant_genes=np.flatnonzero(gconst).tolist(),
              constant_scores=np.flatnonzero(sconst).tolist(),constant_clinical=np.flatnonzero(cconst).tolist())
    return a,b,meta

def ridge_fit(x,y,w,allow_fallback=True):
    X=np.c_[np.ones(len(x)),x];pen=np.ones(X.shape[1]);pen[0]=0.
    b=np.zeros(X.shape[1]);prev=np.average(y,weights=w);b[0]=np.log(prev/(1-prev))
    def obj(beta):
        z=X@beta
        return float(w@(np.logaddexp(0,z)-y*z)+.5*np.sum(pen*beta**2))
    def grad(beta):return X.T@(w*(expit(X@beta)-y))+pen*beta
    used_fallback=False
    for it in range(100):
        p=expit(X@b);g=X.T@(w*(p-y))+pen*b
        if np.max(np.abs(g))<1e-8:return b,dict(iterations=it,fallback=False,grad_max=float(np.max(np.abs(g))))
        H=X.T@((w*p*(1-p))[:,None]*X)+np.diag(pen)
        step=np.linalg.solve(H,g);old=obj(b);rate=1.
        for _ in range(35):
            new=b-rate*step
            if obj(new)<=old-1e-4*rate*(g@step)+1e-12:break
            rate*=.5
        b=new
        if np.max(np.abs(rate*step))<1e-11 and np.max(np.abs(grad(b)))<1e-6:
            return b,dict(iterations=it+1,fallback=False,grad_max=float(np.max(np.abs(grad(b)))))
    if allow_fallback:
        from scipy.optimize import minimize
        z=minimize(obj,b,jac=grad,method='BFGS',options={'gtol':1e-9,'maxiter':1000})
        if np.max(np.abs(grad(z.x)))<1e-6:
            return z.x,dict(iterations=int(z.nit),fallback=True,grad_max=float(np.max(np.abs(grad(z.x)))))
    raise RuntimeError('ridge convergence failed')

def weighted_metrics(y,p,w):
    pos=y==1;neg=y==0
    cmp=(p[pos,None]>p[None,neg]).astype(float)+.5*(p[pos,None]==p[None,neg])
    auc=float(w[pos]@cmp@w[neg]/(w[pos].sum()*w[neg].sum()))
    return auc,float(np.average((y-p)**2,weights=w))

def paired_cv(x,y,w,seed,clinical=False,repeats=5,details=False,preprocess_fn=None):
    if preprocess_fn is None:
        preprocess_fn=preprocess
    chosen=w>0;counts=[np.sum(chosen&(y==k)) for k in [0,1]];nfold=int(min(5,*counts))
    if nfold<2:raise ValueError('fewer than two unique IDs in a response class')
    rng=np.random.default_rng(seed);n=len(y)
    preds=np.full((repeats,n,2),np.nan);folds=np.full((repeats,n),-1,int)
    fitrows=[];preprows=[];fallbacks=0;constant_folds=0;maxgrad=0.
    for rep in range(repeats):
        f=np.full(n,-1,int)
        for label in [0,1]:
            ids=np.flatnonzero(chosen&(y==label));rng.shuffle(ids)
            for k,chunk in enumerate(np.array_split(ids,nfold)):f[chunk]=k
        folds[rep]=f
        assert np.all(f[chosen]>=0) and np.all(f[~chosen]==-1)
        for k in range(nfold):
            tr=np.flatnonzero(chosen&(f!=k));te=np.flatnonzero(chosen&(f==k))
            assert not np.intersect1d(tr,te).size and len(np.unique(y[tr]))==2
            a,b,prep=preprocess_fn(x[tr],x[te],w[tr],clinical)
            constant_folds+=bool(prep['constant_genes'] or prep['constant_scores'] or prep['constant_clinical'])
            if details:preprows.append(dict(repetition=rep,fold=k,training_ids=tr.tolist(),heldout_ids=te.tolist(),**prep))
            for model in [0,1]:
                aa=a[:,:-1] if model==0 else a;bb=b[:,:-1] if model==0 else b
                coef,diag=ridge_fit(aa,y[tr],w[tr]);fallbacks+=diag['fallback'];maxgrad=max(maxgrad,diag['grad_max'])
                preds[rep,te,model]=expit(coef[0]+bb@coef[1:])
                if details:fitrows.append(dict(repetition=rep,fold=k,model=f'M{model}',coefficients=coef.tolist(),**diag))
    assert np.isfinite(preds[:,chosen,:]).all()
    p=preds[:,chosen,:].mean(axis=0);yy=y[chosen];ww=w[chosen]
    a0,b0=weighted_metrics(yy,p[:,0],ww);a1,b1=weighted_metrics(yy,p[:,1],ww)
    metrics=dict(auc_M0=a0,auc_M1=a1,delta_auc=a1-a0,brier_M0=b0,brier_M1=b1,delta_brier=b1-b0,
                 nfold=nfold,unique_patients=int(chosen.sum()),unique_positive=int(counts[1]),unique_negative=int(counts[0]),
                 fallbacks=fallbacks,constant_folds=constant_folds,max_gradient=maxgrad)
    return metrics,preds,folds,fitrows,preprows
