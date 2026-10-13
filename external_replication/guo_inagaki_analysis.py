import json, requests, h5py, tempfile, os, math, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import balanced_accuracy_score

SEED=20261012; NPERM=500000; NBOOT=50000
OUT=Path('analysis/guo_inagaki_locked'); OUT.mkdir(parents=True,exist_ok=True)
API='https://api.dandiarchive.org/api/dandisets/000009/versions/0.220126.1903/assets/'

def dec(v):
    if isinstance(v,bytes): return v.decode('utf-8','ignore')
    return str(v)
def arrstr(ds): return np.array([dec(x) for x in ds[()]])
def alias(s):
    s=str(s).lower()
    if 'alm' in s: return 'ALM'
    if 'thalam' in s or 'vm' in s or 'val' in s: return 'THAL'
    return None

def get_assets():
    A=[]; u=API
    while u:
        j=requests.get(u,timeout=60).json(); A+=j.get('results',[]); u=j.get('next')
    return [dict(asset_id=a.get('asset_id') or a.get('identifier'),path=a.get('path'),size=a.get('size')) for a in A if 'ecephys+ogen.nwb' in a.get('path','')]

def asset_url(asset_id):
    u=f'https://api.dandiarchive.org/api/dandisets/000009/versions/0.220126.1903/assets/{asset_id}/'
    j=requests.get(u,timeout=60).json(); urls=j.get('contentUrl') or j.get('content_url') or []
    if isinstance(urls,str): urls=[urls]
    if not urls:
        j=requests.get(f'https://api.dandiarchive.org/api/assets/{asset_id}/',timeout=60).json(); urls=j.get('contentUrl') or []
        if isinstance(urls,str): urls=[urls]
    return next((x for x in urls if 'amazonaws.com' in x or 's3' in x), urls[0] if urls else None)

def spikes_by_unit(f):
    st=np.asarray(f['units/spike_times'][()],float); idx=np.asarray(f['units/spike_times_index'][()],int)
    starts=np.r_[0,idx[:-1]]
    return [st[a:b] for a,b in zip(starts,idx)]

def session_extract(fp,path,asset_id):
    subj=path.split('/')[0].replace('sub-','')
    with h5py.File(fp,'r') as f:
        # recorded region from electrodes used by units
        locs=arrstr(f['general/extracellular_ephys/electrodes/location'])
        ue=np.asarray(f['units/electrodes'][()],int).reshape(-1)
        used=[]
        for e in ue:
            if 0 <= e < len(locs): used.append(locs[e])
        rregs=[alias(x) for x in used if alias(x)]
        rec=max(set(rregs),key=rregs.count) if rregs else None
        # optogenetic target from subgroup names and locations
        sregs=[]
        if 'general/optogenetics' in f:
            g=f['general/optogenetics']
            for nm in g.keys():
                sg=g[nm]; vals=[nm]
                if 'location' in sg:
                    try: vals.append(dec(sg['location'][()]))
                    except: pass
                for v in vals:
                    a=alias(v)
                    if a:sregs.append(a)
        stim=max(set(sregs),key=sregs.count) if sregs else None
        if rec not in ('ALM','THAL') or stim!=rec:
            return {'asset_id':asset_id,'path':path,'subject':subj,'record_region':rec,'stim_region':stim,'eligible':False,'reason':'region_mismatch'}
        T=f['intervals/trials']
        n=len(T['id'])
        good=np.asarray(T['is_good'][()]).astype(bool)
        stim_present=np.asarray(T['stim_present'][()]).astype(bool)
        typ=np.char.lower(arrstr(T['type']).astype(str))
        resp=np.char.lower(arrstr(T['response']).astype(str))
        per=np.char.lower(arrstr(T['photo_stim_period']).astype(str))
        cue=np.asarray(T['cue_start_time'][()],float)
        left=np.array(['left' in x for x in typ]); right=np.array(['right' in x for x in typ]); choice=left|right
        correct=np.array(['correct' in x and 'incorrect' not in x for x in resp])
        delay=np.array(['delay' in x for x in per])
        ctrl=good & (~stim_present) & correct & choice & np.isfinite(cue)
        if ctrl.sum()<16 or (ctrl&left).sum()<8 or (ctrl&right).sum()<8:
            return {'asset_id':asset_id,'path':path,'subject':subj,'record_region':rec,'stim_region':stim,'eligible':False,'reason':'control_trials','n_ctrl':int(ctrl.sum()),'n_left_ctrl':int((ctrl&left).sum()),'n_right_ctrl':int((ctrl&right).sum())}
        units=spikes_by_unit(f); nu=len(units)
        if nu<5:
            return {'asset_id':asset_id,'path':path,'subject':subj,'record_region':rec,'stim_region':stim,'eligible':False,'reason':'unit_count','n_units':nu}
        ids=np.where(ctrl)[0]; X=np.zeros((len(ids),nu),float); y=right[ids].astype(int)
        for ri,t in enumerate(ids):
            a=cue[t]-0.6; b=cue[t]
            for ui,sp in enumerate(units): X[ri,ui]=np.searchsorted(sp,b)-np.searchsorted(sp,a)
        if len(np.unique(y))<2 or min(np.bincount(y))<8:
            return {'asset_id':asset_id,'path':path,'subject':subj,'record_region':rec,'stim_region':stim,'eligible':False,'reason':'choice_balance'}
        cv=StratifiedKFold(n_splits=5,shuffle=True,random_state=SEED)
        model=make_pipeline(StandardScaler(),LogisticRegression(penalty='l2',C=1,max_iter=2000,random_state=SEED))
        scores=cross_val_score(model,X,y,cv=cv,scoring='balanced_accuracy')
        pred=float(np.mean(scores)-0.5)
        sd=X.std(0,ddof=1); sd[sd==0]=1
        z=(X-X.mean(0))/sd
        diag=float(np.nanmean(np.abs(z[y==1].mean(0)-z[y==0].mean(0))))
        # outcome: matched-site delay photostim impairment, side first then average
        base=good & (~stim_present) & choice
        st=good & stim_present & delay & choice
        if base.sum()<8 or st.sum()<8:
            return {'asset_id':asset_id,'path':path,'subject':subj,'record_region':rec,'stim_region':stim,'eligible':False,'reason':'stim_trials','n_base':int(base.sum()),'n_stim':int(st.sum())}
        imps=[]
        side_counts={}
        for name,mask in [('left',left),('right',right)]:
            b=base&mask; s=st&mask
            side_counts[f'n_{name}_base']=int(b.sum()); side_counts[f'n_{name}_stim']=int(s.sum())
            if b.sum()>0 and s.sum()>0:
                imps.append(float(correct[b].mean()-correct[s].mean()))
        if not imps:
            return {'asset_id':asset_id,'path':path,'subject':subj,'record_region':rec,'stim_region':stim,'eligible':False,'reason':'no_outcome'}
        impairment=float(np.mean(imps))
        return {'asset_id':asset_id,'path':path,'subject':subj,'record_region':rec,'stim_region':stim,'eligible':True,'predictor':pred,'diagnostic_single_unit':diag,'impairment':impairment,'n_units':nu,'n_ctrl':int(ctrl.sum()),'n_stim':int(st.sum()),**side_counts}

def rho(x,y):
    r=spearmanr(x,y).statistic
    return float(r) if np.isfinite(r) else np.nan

def perm_p(x,y,n=NPERM,seed=SEED):
    obs=rho(x,y); rx=rankdata(x); ry=rankdata(y); rx=rx-rx.mean(); ry=ry-ry.mean(); den=np.sqrt((rx*rx).sum()*(ry*ry).sum()); rng=np.random.default_rng(seed); ge=0
    batch=5000
    for start in range(0,n,batch):
        k=min(batch,n-start); P=np.empty((k,len(y)))
        for z in range(k): P[z]=rng.permutation(ry)
        rr=(P@rx)/den; ge+=int(np.sum(rr>=obs-1e-15))
    return obs,float((ge+1)/(n+1))

def bootstrap_ci(x,y,n=NBOOT,seed=SEED+1):
    rng=np.random.default_rng(seed); vals=[]; N=len(x)
    for _ in range(n):
        ix=rng.integers(0,N,N); r=rho(x[ix],y[ix])
        if np.isfinite(r):vals.append(r)
    return [float(np.percentile(vals,2.5)),float(np.percentile(vals,97.5))]

def loo(x,y):
    n=len(x); pred=np.full(n,np.nan)
    for i in range(n):
        tr=np.arange(n)!=i
        m=LinearRegression().fit(x[tr,None],y[tr]); pred[i]=m.predict([[x[i]]])[0]
    mse=np.mean((pred-y)**2); var=np.mean((y-y.mean())**2)
    return {'rho':rho(pred,y),'r2':float(1-mse/var) if var>0 else None,'mae':float(np.mean(np.abs(pred-y)))}

def infer(df,seed=SEED):
    x=df.predictor.to_numpy(float); y=df.impairment.to_numpy(float)
    r,p=perm_p(x,y,seed=seed); ci=bootstrap_ci(x,y,seed=seed+1); lp=loo(x,y) if len(df)>=4 else {'rho':None,'r2':None,'mae':None}
    return {'n':int(len(df)),'n_mice':int(df.subject.nunique()),'rho':r,'p_one_sided':p,'bootstrap95':ci,'loo':lp}

def main():
    assets=get_assets(); rows=[]
    for k,a in enumerate(assets):
        try:
            u=asset_url(a['asset_id'])
            if not u: raise RuntimeError('no content URL')
            with tempfile.NamedTemporaryFile(suffix='.nwb') as tf:
                r=requests.get(u,timeout=180); r.raise_for_status(); tf.write(r.content); tf.flush()
                z=session_extract(tf.name,a['path'],a['asset_id'])
            rows.append(z); print(k+1,'/',len(assets),z)
        except Exception as e:
            rows.append({'asset_id':a['asset_id'],'path':a['path'],'eligible':False,'reason':'exception','error':repr(e)}); print('ERROR',a['path'],repr(e))
    tab=pd.DataFrame(rows); tab.to_csv(OUT/'session_inventory.csv',index=False)
    d=tab[tab.eligible==True].copy() if 'eligible' in tab else pd.DataFrame()
    result={'status':'FROZEN_GUO_INAGAKI_MATCHED_SYSTEM_TEST','dataset':'DANDI 000009 v0.220126.1903','lock':'external_replication/GUO_INAGAKI_2017_LOCK.md','n_assets_ecephys_ogen':len(assets),'n_eligible_sessions':int(len(d)),'n_eligible_mice':int(d.subject.nunique()) if len(d) else 0}
    if len(d)>=3:
        result['primary']=infer(d,SEED)
        success=len(d)>=10 and d.subject.nunique()>=6 and result['primary']['rho']>0 and result['primary']['p_one_sided']<.05 and result['primary']['loo']['r2'] is not None and result['primary']['loo']['r2']>0
        result['primary']['success']=bool(success)
        # mouse-mean cluster sensitivity
        md=d.groupby('subject',as_index=False).agg(predictor=('predictor','mean'),impairment=('impairment','mean'))
        if len(md)>=3: result['mouse_mean_sensitivity']=infer(md,SEED+200)
        # anatomical strata
        result['strata']={}
        for j,reg in enumerate(['ALM','THAL']):
            q=d[d.record_region==reg]
            if len(q)>=5: result['strata'][reg]=infer(q,SEED+500+j*100)
        # diagnostic predictor, descriptive
        if d.diagnostic_single_unit.notna().sum()>=3:
            q=d.dropna(subset=['diagnostic_single_unit']).copy(); q['predictor']=q['diagnostic_single_unit']; result['diagnostic_single_unit']=infer(q,SEED+900)
    else:
        result['primary']={'status':'STRUCTURALLY_UNDERPOWERED_OR_INELIGIBLE'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True))
    print(json.dumps(result,indent=2,sort_keys=True)); print(d.to_string(index=False) if len(d) else 'NO ELIGIBLE')
if __name__=='__main__': main()
