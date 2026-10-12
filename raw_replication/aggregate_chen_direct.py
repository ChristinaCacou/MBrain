import json,csv
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,rankdata

SEED=20261012;NPERM=500000;NBOOT=50000
ROOT=Path('all-session-results');OUT=Path('analysis/chen_direct_final');OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for p in ROOT.rglob('*.json'):
    try:r=json.loads(p.read_text())
    except:continue
    if 'error' not in r:rows.append(r)
def finite(v):
    try:return v is not None and np.isfinite(float(v))
    except:return False
animals=sorted(set(r['animal'] for r in rows));mice=[]
for a in animals:
    q=[r for r in rows if r['animal']==a]
    def mn(k):
        z=[float(r[k]) for r in q if finite(r.get(k))];return float(np.mean(z)) if z else None
    mice.append({'animal':a,'n_sessions':len(q),'coupling':mn('control_bilateral_coupling'),'coupling_within_condition':mn('control_bilateral_within_condition_coupling'),
                 'transferable_local_cv_auc':mn('transferable_local_cv_auc'),'left_neuronal_recovery':mn('left_neuronal_recovery'),'right_neuronal_recovery':mn('right_neuronal_recovery'),
                 'left_behavior_recovery':mn('left_behavior_recovery'),'right_behavior_recovery':mn('right_behavior_recovery'),'n_units_left':mn('n_units_left'),'n_units_right':mn('n_units_right')})
def perm(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);obs=float(spearmanr(x,y).statistic);rx=rankdata(x);ry=rankdata(y);rx-=rx.mean();den=np.sqrt(np.sum(rx*rx));rng=np.random.default_rng(SEED);ge=0;done=0
    while done<NPERM:
        b=min(10000,NPERM-done);idx=np.argsort(rng.random((b,len(y))),axis=1);yp=ry[idx];yp=yp-yp.mean(1,keepdims=True);v=(yp@rx)/(den*np.sqrt(np.sum(yp*yp,axis=1)));ge+=int(np.sum(v>=obs-1e-12));done+=b
    return obs,float((ge+1)/(NPERM+1))
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
    if len(x)<4:return None,None,[]
    for i in range(len(y)):
        tr=np.arange(len(y))!=i;A=np.column_stack([np.ones(tr.sum()),x[tr]]);b=np.linalg.lstsq(A,y[tr],rcond=None)[0];pred[i]=b[0]+b[1]*x[i]
    rho=float(spearmanr(pred,y).statistic);var=np.mean((y-y.mean())**2);r2=float(1-np.mean((pred-y)**2)/var) if var>0 else None;return rho,r2,pred.tolist()
def assoc(ykey,xkey='coupling'):
    q=[m for m in mice if finite(m.get(xkey)) and finite(m.get(ykey))];x=[m[xkey] for m in q];y=[m[ykey] for m in q]
    if len(q)<3:return {'n_mice':len(q)}
    rho,p=perm(x,y);lr,lr2,pred=loo(x,y);return {'n_mice':len(q),'spearman_rho':rho,'one_sided_permutation_p':p,'bootstrap_95_ci':boot(x,y),'loo_spearman_rho':lr,'loo_r2':lr2,'loo_predictions':pred}
result={'analysis_status':'RAW_CHEN_REANALYSIS_FROM_SESSION_SPIKES','dataset':'Chen et al. Cell 2021 / Zenodo 6713616','n_raw_sessions_parsed':len(rows),'n_mice_any':len(mice),
        'direct_predictor':'control-trial bilateral CD-projection coupling from raw spikes, final 400 ms before cue',
        'right_neuronal_recovery':assoc('right_neuronal_recovery'),'right_behavior_recovery':assoc('right_behavior_recovery'),
        'left_neuronal_recovery':assoc('left_neuronal_recovery'),'left_behavior_recovery':assoc('left_behavior_recovery'),
        'within_condition_sensitivity_right_neuronal':assoc('right_neuronal_recovery','coupling_within_condition'),
        'transferable_metric_right_neuronal':assoc('right_neuronal_recovery','transferable_local_cv_auc'),
        'transferable_metric_left_neuronal':assoc('left_neuronal_recovery','transferable_local_cv_auc'),
        'note':'Chen is exploratory because the source publication and prior processed-figure results were already known. Mouse is the inferential unit; sessions are averaged within mouse.'}
(OUT/'result.json').write_text(json.dumps(result,indent=2));(OUT/'mouse_rows.json').write_text(json.dumps(mice,indent=2));(OUT/'session_rows.json').write_text(json.dumps(rows,indent=2))
with open(OUT/'mouse_rows.csv','w',newline='') as f:
    if mice:w=csv.DictWriter(f,fieldnames=list(mice[0].keys()));w.writeheader();w.writerows(mice)
print(json.dumps(result,indent=2))
