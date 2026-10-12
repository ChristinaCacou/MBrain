import json, glob, os
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
SEED=20261012; NPERM=500000; NBOOT=50000
OUT=Path('analysis/chen_raw_final'); OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for f in glob.glob('all-session-results/**/*.json',recursive=True):
    try:
        d=json.load(open(f))
        if isinstance(d,dict) and 'animal' in d: rows.append(d)
    except: pass
# dedupe
bykey={}
for r in rows: bykey[(r.get('animal'),r.get('session'))]=r
rows=list(bykey.values())
valid=[]
for r in rows:
    if 'error' in r: continue
    for hemi in ['left','right']:
        p=r.get('predictor',{}).get(hemi); y=r.get('outcome',{}).get(hemi,{}).get('impairment')
        if p is not None and y is not None and np.isfinite(p) and np.isfinite(y):
            valid.append({'animal':r['animal'],'session':r['session'],'hemisphere':hemi,'predictor':float(p),'impairment':float(y),'n_units':r.get('n_eligible_units',{}).get(hemi)})
# average sessions per mouse x hemisphere
agg=[]
for key in sorted(set((x['animal'],x['hemisphere']) for x in valid)):
    q=[x for x in valid if (x['animal'],x['hemisphere'])==key]
    agg.append({'animal':key[0],'hemisphere':key[1],'n_sessions':len(q),'predictor':float(np.mean([z['predictor'] for z in q])),'impairment':float(np.mean([z['impairment'] for z in q]))})
animals=sorted(set(x['animal'] for x in agg))
# pooled mouse-hemisphere association; block permutations shuffle outcome mouse blocks, keeping L/R values paired within mouse
pairs={a:{x['hemisphere']:x for x in agg if x['animal']==a} for a in animals}
animals2=[a for a in animals if 'left' in pairs[a] and 'right' in pairs[a]]
X=np.array([[pairs[a][h]['predictor'] for h in ['left','right']] for a in animals2],float)
Y=np.array([[pairs[a][h]['impairment'] for h in ['left','right']] for a in animals2],float)
def rho_flat(A,B): return float(spearmanr(A.ravel(),B.ravel()).statistic)
obs=rho_flat(X,Y) if len(animals2)>=3 else np.nan
rng=np.random.default_rng(SEED); ge=0
if np.isfinite(obs):
    for _ in range(NPERM):
        yp=Y[rng.permutation(len(Y)),:]
        if rho_flat(X,yp)>=obs-1e-12: ge+=1
    p=float((ge+1)/(NPERM+1))
else:p=np.nan
boot=[]
if len(animals2)>=3:
    rng=np.random.default_rng(SEED)
    for _ in range(NBOOT):
        ii=rng.integers(0,len(X),len(X)); rr=rho_flat(X[ii],Y[ii])
        if np.isfinite(rr):boot.append(rr)
ci=[float(np.percentile(boot,2.5)),float(np.percentile(boot,97.5))] if boot else [None,None]
# leave-one-mouse-out pooled linear prediction
pred=np.full_like(Y,np.nan,float)
for i in range(len(X)):
    tr=np.arange(len(X))!=i
    xx=X[tr].ravel(); yy=Y[tr].ravel(); A=np.c_[np.ones(len(xx)),xx]
    b=np.linalg.lstsq(A,yy,rcond=None)[0]; pred[i]=b[0]+b[1]*X[i]
loo_rho=rho_flat(pred,Y) if len(X)>=4 else np.nan
loo_r2=float(1-np.mean((pred-Y)**2)/np.mean((Y-Y.mean())**2)) if len(X)>=4 else np.nan
# side-specific descriptive/inferential mouse associations
side={}
for h in ['left','right']:
    q=[x for x in agg if x['hemisphere']==h]; x=np.array([z['predictor'] for z in q]); y=np.array([z['impairment'] for z in q]); r=float(spearmanr(x,y).statistic) if len(q)>=3 else np.nan
    rng=np.random.default_rng(SEED); gg=0
    if np.isfinite(r):
        for _ in range(NPERM):
            if spearmanr(x,rng.permutation(y)).statistic>=r-1e-12:gg+=1
        pp=float((gg+1)/(NPERM+1))
    else:pp=np.nan
    side[h]={'n_mice':len(q),'rho':r,'one_sided_permutation_p':pp}
result={'analysis_status':'EXPLORATORY_RAW_CHEN_REANALYSIS','n_raw_session_files':len(rows),'n_eligible_session_hemisphere_rows':len(valid),'n_mouse_hemisphere_rows':len(agg),'n_complete_bilateral_mice':len(animals2),'pooled_mouse_block_spearman_rho':obs,'pooled_mouse_block_one_sided_p':p,'mouse_block_bootstrap_95_ci':ci,'leave_one_mouse_out_rho':loo_rho,'leave_one_mouse_out_r2':loo_r2,'side_specific':side,'animals_complete':animals2,'note':'Predictor was computed from raw spikes on control trials only in the source-native early-delay interval (-1.7 to -0.9 s relative to cue); causal outcome was same-session, same-hemisphere unilateral photoinhibition impairment. Inference treats mouse as the resampling/permutation block.'}
(OUT/'session_hemisphere_rows.json').write_text(json.dumps(valid,indent=2)); (OUT/'mouse_hemisphere_rows.json').write_text(json.dumps(agg,indent=2)); (OUT/'result.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
