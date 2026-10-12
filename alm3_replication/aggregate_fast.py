import json, glob, csv
from pathlib import Path
import numpy as np
from scipy.stats import rankdata, spearmanr

SEED=20261012
N_PERM=500000
N_BOOT=50000
LOCK='743e81be81064fe7e3bfa500b2a15e0472eae79d'

def finite(x):
    try:return np.isfinite(float(x))
    except:return False

def rho(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    return float(spearmanr(x,y).statistic) if len(x)>=3 else np.nan

def corr_rows(A,b):
    A=np.asarray(A,float);b=np.asarray(b,float)
    A=A-A.mean(axis=1,keepdims=True); b=b-b.mean()
    den=np.sqrt((A*A).sum(axis=1)*(b*b).sum())
    return (A@b)/den

def perm_p(x,y,obs,rng,n=N_PERM,chunk=25000):
    rx=rankdata(x).astype(float); ry=rankdata(y).astype(float)
    ge=0; done=0
    while done<n:
        m=min(chunk,n-done)
        P=np.vstack([rng.permutation(ry) for _ in range(m)])
        vals=corr_rows(P,rx)
        ge += int(np.sum(vals>=obs-1e-15)); done += m
    return (ge+1)/(n+1)

def boot_ci(x,y,rng,n=N_BOOT,chunk=10000):
    x=np.asarray(x,float);y=np.asarray(y,float);N=len(x);vals=[]
    # Bootstrap Spearman with ties handled per resample; n is small, loop is acceptable.
    for _ in range(n):
        idx=rng.integers(0,N,N); rr=rho(x[idx],y[idx])
        if np.isfinite(rr): vals.append(rr)
    return [float(np.quantile(vals,.025)),float(np.quantile(vals,.975))] if vals else [np.nan,np.nan]

def loo_linear(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);n=len(x);pred=np.full(n,np.nan)
    if n<4:return {'rho':np.nan,'r2':np.nan,'predictions':pred.tolist()}
    for i in range(n):
        tr=np.arange(n)!=i; X=np.c_[np.ones(tr.sum()),x[tr]]
        beta=np.linalg.lstsq(X,y[tr],rcond=None)[0]; pred[i]=beta[0]+beta[1]*x[i]
    sse=float(np.sum((y-pred)**2));sst=float(np.sum((y-y.mean())**2))
    return {'rho':rho(pred,y),'r2':float(1-sse/sst) if sst>0 else np.nan,'predictions':pred.tolist()}

files=sorted(glob.glob('downloaded-results/**/*.json',recursive=True)); rows=[]
for f in files:
    try:
        r=json.load(open(f));
        if 'asset_path' in r: rows.append(r)
    except Exception: pass
rows=list({r['asset_path']:r for r in rows}.values())
def explicit_callosotomy(r):
    txt=' '.join(str(r.get(k,'')) for k in ['session_description','subject_description','subject_genotype']).lower()
    return 'callos' in txt or 'commissuro' in txt
eligible=[r for r in rows if r.get('eligible') and finite(r.get('predictor_local_choice_selectivity')) and finite(r.get('behavioral_impairment')) and not explicit_callosotomy(r)]
mice={}
for r in eligible:mice.setdefault(str(r['mouse']),[]).append(r)
mouse_rows=[]
for m,rr in sorted(mice.items()):
    mouse_rows.append({'mouse':m,'n_sessions':len(rr),'predictor':float(np.mean([x['predictor_local_choice_selectivity'] for x in rr])),'outcome':float(np.mean([x['behavioral_impairment'] for x in rr])),'hemispheres':','.join(sorted(set(str(x.get('recorded_hemisphere')) for x in rr)))})
x=np.array([r['predictor'] for r in mouse_rows],float);y=np.array([r['outcome'] for r in mouse_rows],float)
obs=rho(x,y);rng=np.random.default_rng(SEED)
p=perm_p(x,y,obs,rng) if len(x)>=3 and np.isfinite(obs) else np.nan
ci=boot_ci(x,y,rng) if len(x)>=3 else [np.nan,np.nan]
loo=loo_linear(x,y)
sx=np.array([r['predictor_local_choice_selectivity'] for r in eligible],float);sy=np.array([r['behavioral_impairment'] for r in eligible],float)
out={'analysis_status':'FROZEN_EXTERNAL_REPLICATION','dataset':'DANDI 000011 version 0.220126.1907','n_session_files_processed':len(rows),'n_eligible_sessions':len(eligible),'n_mice':len(mouse_rows),'n_explicit_callosotomy_sessions_excluded':sum(explicit_callosotomy(r) for r in rows),'primary_mouse_spearman_rho':obs,'primary_one_sided_permutation_p':p,'mouse_bootstrap_95_ci':ci,'leave_one_mouse_out_rho':loo['rho'],'leave_one_mouse_out_r2':loo['r2'],'success_criterion_met':bool(np.isfinite(obs) and p<.05 and loo['r2']>0),'descriptive_session_spearman_rho':rho(sx,sy),'mouse_rows':mouse_rows,'session_exclusion_summary':{},'lock_commit':LOCK,'seed':SEED,'n_permutations':N_PERM,'n_bootstrap':N_BOOT,'implementation_note':'Vectorized permutation implementation; frozen metric, thresholds, inferential unit and success criterion unchanged.'}
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
