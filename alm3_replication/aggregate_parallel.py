import json,glob
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
SEED=20261012;NPER=500000;NBOOT=50000;OUT=Path('analysis/alm3_parallel');OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for f in glob.glob('all-session-results/**/*.json',recursive=True):
    try:
        d=json.load(open(f));
        if isinstance(d,dict) and 'path' in d:rows.append(d)
    except:pass
bykey={r['path']:r for r in rows};rows=list(bykey.values());valid=[r for r in rows if 'predictor' in r and 'impairment' in r]
by={}
for r in valid:by.setdefault(r['subject'],[]).append(r)
mice=[]
for s,rr in sorted(by.items()):mice.append({'subject':s,'n_sessions':len(rr),'predictor':float(np.mean([x['predictor'] for x in rr])),'impairment':float(np.mean([x['impairment'] for x in rr])),'median_laser_power':float(np.nanmedian([x['median_laser_power'] for x in rr if x.get('median_laser_power') is not None])) if any(x.get('median_laser_power') is not None for x in rr) else None})
def perm(x,y):
    obs=float(spearmanr(x,y).statistic);rng=np.random.default_rng(SEED);ge=0
    for _ in range(NPER):
        if spearmanr(x,rng.permutation(y)).statistic>=obs-1e-12:ge+=1
    return obs,float((ge+1)/(NPER+1))
def boot(x,y):
    rng=np.random.default_rng(SEED);v=[]
    for _ in range(NBOOT):
        ii=rng.integers(0,len(x),len(x));r=spearmanr(x[ii],y[ii]).statistic
        if np.isfinite(r):v.append(r)
    return [float(np.percentile(v,2.5)),float(np.percentile(v,97.5))] if v else [None,None]
def loo(x,y):
    pred=np.full(len(y),np.nan)
    for i in range(len(y)):
        tr=np.arange(len(y))!=i;b=np.linalg.lstsq(np.c_[np.ones(tr.sum()),x[tr]],y[tr],rcond=None)[0];pred[i]=b[0]+b[1]*x[i]
    return float(spearmanr(pred,y).statistic),float(1-np.mean((pred-y)**2)/np.mean((y-y.mean())**2)),pred.tolist()
if len(mice)>=4:
    x=np.array([m['predictor'] for m in mice]);y=np.array([m['impairment'] for m in mice]);rho,p=perm(x,y);lr,r2,pred=loo(x,y)
    result={'analysis_status':'LOCKED_EXTERNAL_REPLICATION','lock_commit':'89f6da46568eee087b67b9ac57abc36e82f71dfc','implementation_note':'parallel execution only; frozen eligibility, predictor, outcome and inference unchanged','n_candidate_assets':len(rows),'n_eligible_sessions':len(valid),'n_mice':len(mice),'spearman_rho':rho,'one_sided_permutation_p':p,'bootstrap_95_ci':boot(x,y),'loo_spearman_rho':lr,'loo_r2':r2,'success_criterion_met':bool(rho>0 and p<0.05 and r2>0),'subjects':[m['subject'] for m in mice]}
else:result={'analysis_status':'LOCKED_EXTERNAL_REPLICATION','lock_commit':'89f6da46568eee087b67b9ac57abc36e82f71dfc','n_candidate_assets':len(rows),'n_eligible_sessions':len(valid),'n_mice':len(mice),'success_criterion_met':False,'note':'fewer than four eligible mice; primary correlation not run'}
(OUT/'session_results.json').write_text(json.dumps(rows,indent=2));(OUT/'mouse_results.json').write_text(json.dumps(mice,indent=2));(OUT/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
