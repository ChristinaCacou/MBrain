"""Execute GUO2017_CROSSLOOP_LOCK.md after the numerical lock was committed.

Primary: pooled within-direction rank association between control-trial late-delay
population decoding and contralateral behavioral impairment caused by reciprocal-node
delay photoinhibition. Mouse is the inferential unit.
"""
from pathlib import Path
from itertools import permutations
import json, requests, h5py, numpy as np, pandas as pd
from remfile import File as RemFile
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.pipeline import make_pipeline

DS='000009'; VER='0.220126.1903'; SEED=20261012; NPERM=500000; NBOOT=50000
OUT=Path('analysis/guo2017_crossloop'); OUT.mkdir(parents=True,exist_ok=True)


def decode(x):
    a=np.asarray(x).reshape(-1); out=[]
    for v in a:
        if isinstance(v,bytes): out.append(v.decode(errors='replace'))
        else: out.append(str(v))
    return np.asarray(out,dtype=object)

def asset_url(aid):
    m=requests.get(f'https://api.dandiarchive.org/api/assets/{aid}/',timeout=60).json(); cu=m.get('contentUrl') or []
    if isinstance(cu,str): cu=[cu]
    return next((u for u in cu if 'amazonaws.com' in u),cu[0] if cu else None)

def list_assets():
    u=f'https://api.dandiarchive.org/api/dandisets/{DS}/versions/{VER}/assets/?page_size=200'; out=[]
    while u:
        j=requests.get(u,timeout=60).json(); out.extend(j['results']); u=j.get('next')
    return out

def rec_region(f):
    p='general/extracellular_ephys/electrodes/location'
    if p not in f:return None
    vals=decode(f[p])
    s=' '.join(vals.tolist())
    if 'brain_region: ALM' in s:return 'ALM'
    if 'brain_region: Thalamus' in s:return 'Thalamus'
    return None

def target_site(f):
    p='general/optogenetics'
    if p not in f:return None
    names=list(f[p].keys())
    return names[0] if len(names)==1 else '|'.join(sorted(names))

def unit_spikes(f):
    st=np.asarray(f['units/spike_times']).astype(float).reshape(-1)
    ix=np.asarray(f['units/spike_times_index']).astype(int).reshape(-1)
    starts=np.r_[0,ix[:-1]]
    return [st[a:b] for a,b in zip(starts,ix)]

def count_matrix(spikes, centers, before, after):
    X=np.empty((len(centers),len(spikes)),float)
    lo=centers+before; hi=centers+after
    for j,s in enumerate(spikes):
        X[:,j]=np.searchsorted(s,hi,side='left')-np.searchsorted(s,lo,side='left')
    return X

def cv_auc(X,y):
    y=np.asarray(y,int); nmin=int(min((y==0).sum(),(y==1).sum()))
    if nmin<3 or X.shape[1]<1:return np.nan
    k=min(5,nmin)
    cv=StratifiedKFold(n_splits=k,shuffle=True,random_state=SEED)
    model=make_pipeline(StandardScaler(),LogisticRegression(C=1,solver='liblinear',max_iter=2000,random_state=SEED))
    try:
        p=cross_val_predict(model,X,y,cv=cv,method='predict_proba')[:,1]
        return float(roc_auc_score(y,p))
    except Exception:
        return np.nan

def process_asset(a):
    aid=a['asset_id']; url=asset_url(aid); r={'asset_id':aid,'path':a.get('path')}
    if not url:return r
    try:
      with h5py.File(RemFile(url),'r') as f:
        rr=rec_region(f); ts=target_site(f)
        if not ((rr=='ALM' and ts=='left-Thalamus') or (rr=='Thalamus' and ts=='left-ALM')): return r
        r['recording_region']=rr; r['target_site']=ts
        sid=f['general/subject/subject_id'][()]; r['subject']=sid.decode() if isinstance(sid,bytes) else str(sid)
        t=f['intervals/trials']
        typ=decode(t['type']); resp=decode(t['response']); stim=np.asarray(t['stim_present']).reshape(-1).astype(bool)
        per=decode(t['photo_stim_period']); cue=np.asarray(t['cue_start_time']).reshape(-1).astype(float); pole=np.asarray(t['pole_in_time']).reshape(-1).astype(float)
        ctrl=(~stim)&np.isin(typ,['lick left','lick right'])&np.isfinite(cue)
        y=(typ[ctrl]=='lick right').astype(int)
        sp=unit_spikes(f); X=count_matrix(sp,cue[ctrl],-0.5,0.0)
        r['late_auc']=cv_auc(X,y); r['n_control_decode']=int(ctrl.sum()); r['n_units']=int(len(sp))
        ctrl_early=(~stim)&np.isin(typ,['lick left','lick right'])&np.isfinite(pole)
        ye=(typ[ctrl_early]=='lick right').astype(int); Xe=count_matrix(sp,pole[ctrl_early],0.0,0.5)
        r['early_auc']=cv_auc(Xe,ye)
        contra=(typ=='lick right')
        c=contra&(~stim)&(resp!='')
        s=contra&stim&(np.char.lower(np.char.strip(per.astype(str)))=='delay')&(resp!='')
        r['n_behavior_control']=int(c.sum());r['n_behavior_stim']=int(s.sum())
        if c.sum()>=10 and s.sum()>=10:
            cor=np.char.lower(np.char.strip(resp.astype(str)))=='correct'
            r['control_accuracy']=float(cor[c].mean());r['stim_accuracy']=float(cor[s].mean());r['impairment']=float(r['control_accuracy']-r['stim_accuracy'])
        else:r['impairment']=np.nan
    except Exception as e:r['error']=repr(e)
    return r

def centered_ranks(v):
    q=rankdata(np.asarray(v,float),method='average');return q-q.mean()

def pooled_stat(df,pred='late_auc'):
    xs=[];ys=[]
    for g in ['ALM','Thalamus']:
        q=df[df.recording_region==g]
        if len(q)<2:continue
        xs.extend(centered_ranks(q[pred]).tolist());ys.extend(centered_ranks(q.impairment).tolist())
    if len(xs)<3:return np.nan
    return float(np.corrcoef(xs,ys)[0,1])

def perm_primary(df,obs):
    rng=np.random.default_rng(SEED); groups=[]
    for g in ['ALM','Thalamus']:
        q=df[df.recording_region==g];groups.append((centered_ranks(q.late_auc),centered_ranks(q.impairment)))
    x=np.concatenate([z[0] for z in groups]);denx=np.sqrt(np.sum(x*x));ge=0;done=0
    while done<NPERM:
        b=min(10000,NPERM-done); yparts=[]
        for _,y in groups:
            keys=rng.random((b,len(y))); idx=np.argsort(keys,axis=1);yparts.append(y[idx])
        Y=np.concatenate(yparts,axis=1); vals=(Y@x)/(np.sqrt(np.sum(Y*Y,axis=1))*denx);ge+=int(np.sum(vals>=obs-1e-15));done+=b
    return float((ge+1)/(NPERM+1))

def boot_ci(df):
    rng=np.random.default_rng(SEED+1); vals=[]
    groups=[df[df.recording_region==g].reset_index(drop=True) for g in ['ALM','Thalamus']]
    for _ in range(NBOOT):
        parts=[q.iloc[rng.integers(0,len(q),len(q))] for q in groups]
        vals.append(pooled_stat(pd.concat(parts,ignore_index=True)))
    a=np.asarray(vals,float);return [float(x) for x in np.nanpercentile(a,[2.5,97.5])]

def exact_dir(q):
    x=q.late_auc.to_numpy(float);y=q.impairment.to_numpy(float);obs=float(spearmanr(x,y).statistic);ge=0;tot=0
    for yp in permutations(y):
        ge+=int(float(spearmanr(x,np.asarray(yp,float)).statistic)>=obs-1e-15);tot+=1
    return {'n_mice':int(len(q)),'rho':obs,'exact_one_sided_p':float(ge/tot)}

def loo(df):
    n=len(df);pred=np.full(n,np.nan);base=np.full(n,np.nan);y=df.impairment.to_numpy(float);x=df.late_auc.to_numpy(float);g=(df.recording_region=='Thalamus').astype(float).to_numpy()
    for i in range(n):
        tr=np.arange(n)!=i; X=np.column_stack([x[tr],g[tr]]);m=LinearRegression().fit(X,y[tr]);pred[i]=m.predict([[x[i],g[i]]])[0]
        same=tr&(g==g[i]);base[i]=y[same].mean()
    den=np.sum((y-base)**2);r2=float(1-np.sum((y-pred)**2)/den) if den>0 else np.nan
    return {'spearman_rho':float(spearmanr(pred,y).statistic),'r2_vs_fold_stratum_mean':r2,'predictions':pred.tolist(),'baseline_predictions':base.tolist()}

def main():
    rows=[]
    for a in list_assets():
        if not a.get('path','').endswith('.nwb'):continue
        r=process_asset(a)
        if r.get('recording_region') in ('ALM','Thalamus'):rows.append(r)
        if r.get('recording_region'):print(r.get('subject'),r.get('recording_region'),r.get('target_site'),r.get('late_auc'),r.get('impairment'))
    sess=pd.DataFrame(rows);sess.to_csv(OUT/'session_results.csv',index=False)
    eligible=sess[np.isfinite(pd.to_numeric(sess.late_auc,errors='coerce'))&np.isfinite(pd.to_numeric(sess.impairment,errors='coerce'))].copy()
    mouse=eligible.groupby(['subject','recording_region','target_site'],as_index=False).agg(late_auc=('late_auc','mean'),early_auc=('early_auc','mean'),impairment=('impairment','mean'),n_sessions=('asset_id','size'),mean_control_accuracy=('control_accuracy','mean'),mean_stim_accuracy=('stim_accuracy','mean'))
    mouse.to_csv(OUT/'mouse_results.csv',index=False)
    obs=pooled_stat(mouse);p=perm_primary(mouse,obs) if len(mouse)>=3 else np.nan;ci=boot_ci(mouse) if len(mouse)>=3 else [np.nan,np.nan];pred=loo(mouse) if len(mouse)>=3 else {}
    sec={}
    for g in ['ALM','Thalamus']:
        q=mouse[mouse.recording_region==g]
        sec[g]=exact_dir(q) if 2<len(q)<=8 else {'n_mice':int(len(q))}
    early=pooled_stat(mouse,'early_auc') if mouse.early_auc.notna().all() else np.nan
    counts=mouse.groupby('recording_region').size().to_dict();success=bool(len(mouse)>=12 and all(counts.get(g,0)>=5 for g in ['ALM','Thalamus']) and np.isfinite(obs) and obs>0 and np.isfinite(p) and p<.05 and pred.get('r2_vs_fold_stratum_mean',-np.inf)>0)
    result={'analysis_status':'FROZEN_NUMERICAL_GUO2017_CROSSLOOP_TEST','lock_commit':'561215ed975cff693af2e93795d81dbbeae7aafc','dataset':f'DANDI {DS} version {VER}','n_eligible_sessions':int(len(eligible)),'n_eligible_mice':int(len(mouse)),'mice_by_recording_region':counts,'primary':{'pooled_within_stratum_rank_association':obs,'one_sided_permutation_p':p,'stratified_bootstrap95':ci,'success':success,**pred},'secondary_by_recording_region':sec,'temporal_control':{'early_window_pooled_association':early,'late_minus_early':float(obs-early) if np.isfinite(early) else None},'interpretation':'Frozen numerical source-data test. The publication qualitative conclusions were public, so this is not a blinded preregistration.'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
