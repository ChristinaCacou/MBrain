import json, glob, math, csv
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

SEED=20261012
N_PERM=500000
N_BOOT=50000

def finite(x):
    try:return np.isfinite(float(x))
    except:return False

def rho(x,y):
    if len(x)<3:return np.nan
    return float(spearmanr(np.asarray(x,float),np.asarray(y,float)).statistic)

def loo_linear(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(x); pred=np.full(n,np.nan)
    if n<4:return {'rho':np.nan,'r2':np.nan,'predictions':pred.tolist()}
    for i in range(n):
        tr=np.arange(n)!=i
        X=np.c_[np.ones(tr.sum()),x[tr]]
        beta=np.linalg.lstsq(X,y[tr],rcond=None)[0]
        pred[i]=beta[0]+beta[1]*x[i]
    rr=rho(pred,y)
    sse=float(np.sum((y-pred)**2)); sst=float(np.sum((y-y.mean())**2))
    r2=float(1-sse/sst) if sst>0 else np.nan
    return {'rho':rr,'r2':r2,'predictions':pred.tolist()}

files=sorted(glob.glob('downloaded-results/**/*.json',recursive=True))
rows=[]
for f in files:
    try:
        r=json.load(open(f))
        if 'asset_path' in r: rows.append(r)
    except Exception: pass
# Remove exact duplicates if artifact download creates nested repeats.
uniq={r['asset_path']:r for r in rows}; rows=list(uniq.values())
# Lock rule: callosotomy-only sessions excluded only if explicitly identified by metadata.
def explicit_callosotomy(r):
    txt=' '.join(str(r.get(k,'')) for k in ['session_description','subject_description','subject_genotype']).lower()
    return 'callos' in txt or 'commissuro' in txt
eligible=[r for r in rows if r.get('eligible') and finite(r.get('predictor_local_choice_selectivity')) and finite(r.get('behavioral_impairment')) and not explicit_callosotomy(r)]
# equal-session-within-mouse aggregation, as frozen
mice={}
for r in eligible:
    mice.setdefault(str(r['mouse']),[]).append(r)
mouse_rows=[]
for m,rr in sorted(mice.items()):
    mouse_rows.append({'mouse':m,'n_sessions':len(rr),
                       'predictor':float(np.mean([x['predictor_local_choice_selectivity'] for x in rr])),
                       'outcome':float(np.mean([x['behavioral_impairment'] for x in rr])),
                       'hemispheres':','.join(sorted(set(str(x.get('recorded_hemisphere')) for x in rr)))})
x=np.array([r['predictor'] for r in mouse_rows],float); y=np.array([r['outcome'] for r in mouse_rows],float)
obs=rho(x,y)
rng=np.random.default_rng(SEED)
if len(x)>=3 and np.isfinite(obs):
    ge=0
    for _ in range(N_PERM):
        if rho(x,rng.permutation(y))>=obs-1e-15: ge+=1
    p=(ge+1)/(N_PERM+1)
    boots=[]
    for _ in range(N_BOOT):
        idx=rng.integers(0,len(x),len(x))
        rr=rho(x[idx],y[idx])
        if np.isfinite(rr): boots.append(rr)
    ci=[float(np.quantile(boots,.025)),float(np.quantile(boots,.975))] if boots else [np.nan,np.nan]
else:
    p=np.nan;ci=[np.nan,np.nan]
loo=loo_linear(x,y)
# descriptive session-level association
sx=np.array([r['predictor_local_choice_selectivity'] for r in eligible],float)
sy=np.array([r['behavioral_impairment'] for r in eligible],float)
session_rho=rho(sx,sy)
criterion=bool(np.isfinite(obs) and p<0.05 and loo['r2']>0)
out={
 'analysis_status':'FROZEN_EXTERNAL_REPLICATION',
 'dataset':'DANDI 000011 version 0.220126.1907',
 'n_session_files_processed':len(rows),
 'n_eligible_sessions':len(eligible),
 'n_mice':len(mouse_rows),
 'n_explicit_callosotomy_sessions_excluded':sum(explicit_callosotomy(r) for r in rows),
 'primary_mouse_spearman_rho':obs,
 'primary_one_sided_permutation_p':p,
 'mouse_bootstrap_95_ci':ci,
 'leave_one_mouse_out_rho':loo['rho'],
 'leave_one_mouse_out_r2':loo['r2'],
 'success_criterion_met':criterion,
 'descriptive_session_spearman_rho':session_rho,
 'mouse_rows':mouse_rows,
 'session_exclusion_summary':{},
 'lock_commit':'743e81be81064fe7e3bfa500b2a15e0472eae79d',
 'seed':SEED,'n_permutations':N_PERM,'n_bootstrap':N_BOOT
}
for r in rows:
    key='eligible' if r in eligible else ('explicit_callosotomy' if explicit_callosotomy(r) else str(r.get('reason','other')))
    out['session_exclusion_summary'][key]=out['session_exclusion_summary'].get(key,0)+1
Path('analysis/alm3_replication').mkdir(parents=True,exist_ok=True)
Path('analysis/alm3_replication/result.json').write_text(json.dumps(out,indent=2,allow_nan=True))
with open('analysis/alm3_replication/mouse_table.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['mouse','n_sessions','predictor','outcome','hemispheres']);w.writeheader();w.writerows(mouse_rows)
with open('analysis/alm3_replication/session_table.csv','w',newline='') as f:
    cols=['asset_path','mouse','recorded_hemisphere','eligible','reason','predictor_local_choice_selectivity','behavioral_impairment','n_eligible_units','n_valid_trials','n_control_trials','n_same_hemi_stim_trials','n_contra_stim_trials','n_bilateral_stim_trials','n_photostim_events']
    w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows([{k:r.get(k) for k in cols} for r in rows])
print(json.dumps(out,indent=2,allow_nan=True))
