import json,csv,itertools
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,rankdata

SEED=20261012; NPERM=500000; NBOOT=50000
ROOT=Path('all-session-results'); OUT=Path('analysis/alm3_final'); OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for p in ROOT.rglob('*.json'):
    try:r=json.loads(p.read_text())
    except Exception:continue
    if 'error' not in r:rows.append(r)

def finite(v):
    try:return v is not None and np.isfinite(float(v))
    except Exception:return False
eligible=[r for r in rows if finite(r.get('predictor_cv_auc')) and finite(r.get('ipsi_recovery'))]
animals=sorted(set(r['animal'] for r in eligible)); mice=[]
for a in animals:
    q=[r for r in eligible if r['animal']==a]
    def mean(k):
        z=[float(r[k]) for r in q if finite(r.get(k))]
        return float(np.mean(z)) if z else None
    mice.append({'animal':a,'n_sessions':len(q),'predictor_cv_auc':mean('predictor_cv_auc'),'ipsi_recovery':mean('ipsi_recovery'),
                 'ipsi_behavior_recovery':mean('ipsi_behavior_recovery'),'contra_recovery':mean('contra_recovery'),
                 'contra_behavior_recovery':mean('contra_behavior_recovery'),'n_units':mean('n_units'),'control_balanced_accuracy':mean('control_balanced_accuracy')})

def perm(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);obs=float(spearmanr(x,y).statistic);n=len(x)
    rx=rankdata(x);ry=rankdata(y);rx-=rx.mean();den=np.sqrt(np.sum(rx*rx));ge=0
    if n<=9:
        tot=0
        for idx in itertools.permutations(range(n)):
            yp=ry[list(idx)];yp-=yp.mean();v=float(yp@rx/(den*np.sqrt(np.sum(yp*yp))));ge+=v>=obs-1e-12;tot+=1
        return obs,float(ge/tot),tot
    rng=np.random.default_rng(SEED);done=0;batch=10000
    while done<NPERM:
        b=min(batch,NPERM-done);idx=np.argsort(rng.random((b,n)),axis=1);yp=ry[idx];yp=yp-yp.mean(1,keepdims=True)
        v=(yp@rx)/(den*np.sqrt(np.sum(yp*yp,axis=1)));ge+=int(np.sum(v>=obs-1e-12));done+=b
    return obs,float((ge+1)/(NPERM+1)),NPERM

def boot(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);rng=np.random.default_rng(SEED);v=[]
    for _ in range(NBOOT):
        idx=rng.integers(0,len(x),len(x))
        if len(np.unique(x[idx]))<2 or len(np.unique(y[idx]))<2:continue
        z=spearmanr(x[idx],y[idx]).statistic
        if np.isfinite(z):v.append(float(z))
    return [float(np.percentile(v,2.5)),float(np.percentile(v,97.5))] if v else [None,None]
def loo(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);pred=np.full(len(y),np.nan)
    if len(y)<4:return None,None,[]
    for i in range(len(y)):
        tr=np.arange(len(y))!=i;A=np.column_stack([np.ones(tr.sum()),x[tr]]);b=np.linalg.lstsq(A,y[tr],rcond=None)[0];pred[i]=b[0]+b[1]*x[i]
    rho=float(spearmanr(pred,y).statistic);var=np.mean((y-y.mean())**2);r2=float(1-np.mean((pred-y)**2)/var) if var>0 else None
    return rho,r2,pred.tolist()
def assoc(ykey):
    q=[m for m in mice if finite(m.get('predictor_cv_auc')) and finite(m.get(ykey))];x=[m['predictor_cv_auc'] for m in q];y=[m[ykey] for m in q]
    if len(q)<3:return {'n_mice':len(q)}
    rho,p,nperm=perm(x,y);lr,lr2,pred=loo(x,y)
    return {'n_mice':len(q),'spearman_rho':rho,'one_sided_permutation_p':p,'n_permutations':nperm,'bootstrap_95_ci':boot(x,y),'loo_spearman_rho':lr,'loo_r2':lr2,'loo_predictions':pred}
def baseline(k):
    q=[m for m in mice if finite(m.get(k)) and finite(m.get('ipsi_recovery'))]
    return {'n_mice':len(q),'spearman_rho':float(spearmanr([m[k] for m in q],[m['ipsi_recovery'] for m in q]).statistic)} if len(q)>=3 else {'n_mice':len(q)}
primary=assoc('ipsi_recovery')
result={'analysis_status':'FROZEN_NUMERICAL_EXTERNAL_REPLICATION','dataset':'CRCNS alm-3 / Li et al. Nature 2016','lock':'raw_replication/analysis_lock.json',
        'n_session_files_parsed':len(rows),'n_primary_eligible_sessions':len(eligible),'n_primary_eligible_mice':len(mice),
        'primary_predictor':'control-trial cross-fitted population choice-code AUC, final 400 ms before cue',
        'primary_outcome':'ipsilateral early-delay photoinhibition (type 4) end-of-delay neuronal recovery ratio',
        'primary':primary,'secondary_behavior':assoc('ipsi_behavior_recovery'),'negative_control_contralateral_neuronal':assoc('contra_recovery'),
        'negative_control_contralateral_behavior':assoc('contra_behavior_recovery'),'baselines':{'n_units':baseline('n_units'),'control_behavior':baseline('control_balanced_accuracy')},
        'success_criterion_met':bool(primary.get('spearman_rho',-1)>0 and primary.get('one_sided_permutation_p',1)<0.05 and primary.get('loo_r2',-1)>0),
        'note':'Mouse is the inferential unit. Sessions are averaged equally within mouse; no imputation.'}
(OUT/'result.json').write_text(json.dumps(result,indent=2))
(OUT/'mouse_rows.json').write_text(json.dumps(mice,indent=2))
(OUT/'session_rows.json').write_text(json.dumps(rows,indent=2))
with open(OUT/'mouse_rows.csv','w',newline='') as f:
    if mice:
        w=csv.DictWriter(f,fieldnames=list(mice[0].keys()));w.writeheader();w.writerows(mice)
print(json.dumps(result,indent=2))
