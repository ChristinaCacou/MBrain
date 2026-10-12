from pathlib import Path
import json, zipfile, urllib.request, re
import numpy as np
from scipy.io import loadmat
from scipy.stats import rankdata, spearmanr, pearsonr

SEED=20261011
N_PERM=500000
N_BOOT=50000
URL='https://zenodo.org/records/6713616/files/analyses_scripts.zip?download=1'
ZIP=Path('chen_validation/analyses_scripts.zip')
ROOT=Path('chen_validation/analyses_scripts')
OUT=Path('analysis/chen_bridge'); OUT.mkdir(parents=True,exist_ok=True)
if not ZIP.exists(): urllib.request.urlretrieve(URL,ZIP)
ROOT.mkdir(parents=True,exist_ok=True)
if not any(ROOT.iterdir()):
    with zipfile.ZipFile(ZIP) as z: z.extractall(ROOT)

def native(x):
    if isinstance(x,np.ndarray):
        if x.size==1: return native(x.item())
        if np.issubdtype(x.dtype,np.number): return x.astype(float).ravel().tolist()
        if x.dtype==object: return [native(z) for z in x.ravel()]
    if isinstance(x,np.generic): return x.item()
    return x

def walk(obj,path='root',acc=None):
    if acc is None: acc=[]
    if hasattr(obj,'_fieldnames'):
        d={f:getattr(obj,f) for f in obj._fieldnames}; props=d.get('properties'); typ=d.get('type')
        rec={'path':path}
        if typ is not None:
            try: rec['type']=str(np.asarray(typ).squeeze().item())
            except: rec['type']=str(typ)
        if props is not None and hasattr(props,'_fieldnames'):
            for key in ['XData','YData','String']:
                if hasattr(props,key):
                    try: val=native(getattr(props,key))
                    except: continue
                    if val not in (None,[],""): rec[key]=val
        if len(rec)>1: acc.append(rec)
        for f,v in d.items():
            if f!='properties': walk(v,path+'.'+f,acc)
    elif isinstance(obj,np.ndarray) and (obj.dtype==object or obj.dtype.names):
        for i,v in enumerate(obj.ravel()): walk(v,f'{path}[{i}]',acc)
    return acc

def fig_objects(name):
    p=next(ROOT.rglob(name)); d=loadmat(p,squeeze_me=True,struct_as_record=False); out=[]
    for k,v in d.items():
        if not k.startswith('__'): walk(v,k,out)
    return out

def scatters_by_child(objs, child, n):
    token=f'children[{child}]'
    return [r for r in objs if r.get('type')=='specgraph.scattergroup' and token in r.get('path','') and len(r.get('XData',[]))==n and len(r.get('YData',[]))==n]

def one_sided_perm_spearman(x,y,n=N_PERM):
    x=np.asarray(x,float); y=np.asarray(y,float); obs=float(spearmanr(x,y).statistic)
    rx=rankdata(x); ry=rankdata(y); rx=rx-rx.mean(); denx=np.sqrt(np.sum(rx*rx)); rng=np.random.default_rng(SEED)
    ge=0; done=0; batch=10000
    while done<n:
        b=min(batch,n-done)
        # each row is an independent permutation of y ranks
        keys=rng.random((b,len(y))); idx=np.argsort(keys,axis=1); yp=ry[idx]; yp=yp-yp.mean(axis=1,keepdims=True)
        vals=(yp@rx)/(denx*np.sqrt(np.sum(yp*yp,axis=1)))
        ge += int(np.sum(vals>=obs-1e-12)); done += b
    return obs,float((ge+1)/(n+1))

def bootstrap_rho(x,y,n=N_BOOT):
    x=np.asarray(x,float); y=np.asarray(y,float); rng=np.random.default_rng(SEED); vals=[]
    for _ in range(n):
        idx=rng.integers(0,len(x),len(x))
        if len(np.unique(x[idx]))<2 or len(np.unique(y[idx]))<2: continue
        vals.append(float(spearmanr(x[idx],y[idx]).statistic))
    return [float(np.percentile(vals,2.5)),float(np.percentile(vals,97.5))]

def loo(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); pred=np.full(len(y),np.nan)
    for i in range(len(y)):
        tr=np.arange(len(y))!=i; X=np.column_stack([np.ones(tr.sum()),x[tr]])
        b=np.linalg.lstsq(X,y[tr],rcond=None)[0]; pred[i]=b[0]+b[1]*x[i]
    rho=float(spearmanr(pred,y).statistic); mse=float(np.mean((pred-y)**2)); var=float(np.mean((y-y.mean())**2)); r2=float(1-mse/var)
    return rho,r2,pred

# Mouse-level processed quantities are encoded in Figure S5A-E by the source MATLAB script.
# Source ordering: child 2 = cross-mouse neuronal recovery (SiL->right-ALM inhibition, SiR->left-ALM inhibition),
# child 3 = corresponding behavioral recovery, child 4 = asymmetry (SiL-SiR) vs control-trial bilateral rank correlation.
s5=fig_objects('FigureS5A_E.fig')
c2=scatters_by_child(s5,2,17); c3=scatters_by_child(s5,3,17); c4=scatters_by_child(s5,4,17)
if not (len(c2)==2 and len(c3)==2 and len(c4)==1): raise RuntimeError(('unexpected S5 structure',len(c2),len(c3),len(c4)))
SiL=np.asarray(c2[0]['XData'],float); neur_R=np.asarray(c2[0]['YData'],float)
SiR=np.asarray(c2[1]['XData'],float); neur_L=np.asarray(c2[1]['YData'],float)
beh_R=np.asarray(c3[0]['YData'],float); beh_L=np.asarray(c3[1]['YData'],float)
asym=np.asarray(c4[0]['XData'],float); rankCC=np.asarray(c4[0]['YData'],float)
if not np.allclose(asym,SiL-SiR,rtol=0,atol=1e-12): raise RuntimeError('S5 alignment check failed')

# Session-level cross-check from source figures: S4 control-trial rank correlation and Fig4 perturbation recovery share standard-task session order.
s4=fig_objects('FigureS4CF.fig'); f4=fig_objects('Figure4G_I.fig')
s4sc=[r for r in s4 if r.get('type')=='specgraph.scattergroup' and len(r.get('XData',[]))==49]
f4sc=[r for r in f4 if r.get('type')=='specgraph.scattergroup' and len(r.get('XData',[]))==49]
if len(s4sc)<1 or len(f4sc)<2: raise RuntimeError(('unexpected session structure',len(s4sc),len(f4sc)))
rankCC_session=np.asarray(s4sc[0]['YData'],float)
neur_L_session=np.asarray(f4sc[0]['XData'],float); beh_L_session=np.asarray(f4sc[0]['YData'],float)
neur_R_session=np.asarray(f4sc[1]['XData'],float); beh_R_session=np.asarray(f4sc[1]['YData'],float)

def summarize(x,y):
    sr,pp=one_sided_perm_spearman(x,y); pr=float(pearsonr(x,y).statistic); prp=float(pearsonr(x,y).pvalue); lr,lr2,pred=loo(x,y)
    return {'n':len(x),'spearman_rho':sr,'one_sided_permutation_p':pp,'pearson_r':pr,'pearson_two_sided_p':prp,'bootstrap_95_ci_spearman':bootstrap_rho(x,y),'loo_spearman_rho':lr,'loo_r2':lr2,'loo_predictions':pred.tolist()}

result={
 'dataset':'Chen et al., Cell 2021; Zenodo 6713616',
 'analysis_status':'source-constrained exploratory numerical reanalysis designed after inspection of the authors published scripts/figures; not preregistered and not outcome-blind to the publication',
 'observational_predictor':'control-trial bilateral ALM CD-projection rank correlation (rankCC) aggregated by mouse in source Figure S5',
 'mouse_level':{
   'right_ALM_neuronal_recovery':summarize(rankCC,neur_R),
   'right_ALM_behavior_recovery':summarize(rankCC,beh_R),
   'left_ALM_neuronal_recovery':summarize(rankCC,neur_L),
   'left_ALM_behavior_recovery':summarize(rankCC,beh_L),
 },
 'session_level_descriptive':{
   'right_ALM_neuronal_recovery':{'n':49,'spearman_rho':float(spearmanr(rankCC_session,neur_R_session).statistic),'two_sided_p':float(spearmanr(rankCC_session,neur_R_session).pvalue)},
   'right_ALM_behavior_recovery':{'n':49,'spearman_rho':float(spearmanr(rankCC_session,beh_R_session).statistic),'two_sided_p':float(spearmanr(rankCC_session,beh_R_session).pvalue)},
   'left_ALM_neuronal_recovery':{'n':49,'spearman_rho':float(spearmanr(rankCC_session,neur_L_session).statistic),'two_sided_p':float(spearmanr(rankCC_session,neur_L_session).pvalue)},
   'left_ALM_behavior_recovery':{'n':49,'spearman_rho':float(spearmanr(rankCC_session,beh_L_session).statistic),'two_sided_p':float(spearmanr(rankCC_session,beh_L_session).pvalue)},
 },
 'alignment_checks':{'n_mice':len(rankCC),'max_abs_asymmetry_minus_SiL_minus_SiR':float(np.max(np.abs(asym-(SiL-SiR))))},
 'scope_note':'The source archive exposes scripts and MATLAB figure objects but not a compact raw/model summary table. Values here are author-derived processed quantities recovered from the released figure objects. Mouse is the primary independent unit; the 49-session analysis is descriptive because sessions repeat mice.'
}
(OUT/'result.json').write_text(json.dumps(result,indent=2))
import csv
with open(OUT/'mouse_level.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['mouse_index','rankCC','SiL','SiR','asymmetry','right_neuronal_recovery','left_neuronal_recovery','right_behavior_recovery','left_behavior_recovery'])
    for i in range(len(rankCC)): w.writerow([i+1,rankCC[i],SiL[i],SiR[i],asym[i],neur_R[i],neur_L[i],beh_R[i],beh_L[i]])
print(json.dumps(result,indent=2))
