import hashlib, itertools, json, math, shutil
from pathlib import Path

import mat73
import numpy as np
import pandas as pd
from remotezip import RemoteZip
from scipy.io import loadmat
from scipy.stats import rankdata, spearmanr

SEED=20261011
URL='https://zenodo.org/records/17973874/files/2026_Ghaderi_NCOMMS_Zenodo_data_code.zip?download=1'
TARGETS={
 'decoder':'processed_data/Decoding_NbNeurons5_50.mat',
 'opto':'processed_data/Optoinhibition_mat.mat',
}
REGIONS=['A1','wS1','wS2','wM2','ALM']
EPOCHS={'audio':(-0.05,0.30),'delay':(0.40,0.85),'whisker':(0.95,1.20)}
OUT=Path('analysis/ghaderi_matched_validation'); OUT.mkdir(parents=True,exist_ok=True)
TMP=Path('ghaderi_validation/_tmp'); TMP.mkdir(parents=True,exist_ok=True)


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def extract_targets():
    out={}
    with RemoteZip(URL) as rz:
        names=rz.namelist()
        for key,suffix in TARGETS.items():
            matches=[n for n in names if n.endswith(suffix)]
            if len(matches)!=1: raise RuntimeError((key,matches))
            name=matches[0]; local=TMP/Path(suffix).name
            with rz.open(name) as src, open(local,'wb') as dst: shutil.copyfileobj(src,dst,1<<20)
            out[key]=local
            print(key,name,local.stat().st_size,sha256(local),flush=True)
    return out


def nonempty_cell(x):
    if x is None: return False
    if isinstance(x,str): return len(x.strip())>0
    try: return np.asarray(x).size>0
    except Exception: return True


def decoder_metrics(path):
    d=loadmat(path,squeeze_me=True,struct_as_record=False)
    A=d['Accuracy']; nb=getattr(A,'Nb_neurons50')
    times=np.round(np.arange(-0.99,2.0001,0.01),2)
    rows=[]
    for region in REGIONS:
        arr=np.asarray(getattr(nb,region),dtype=float)
        if arr.ndim==2:
            if arr.shape[1]==len(times): sr=arr
            elif arr.shape[0]==len(times): sr=arr.T
            else: raise RuntimeError((region,'decoder shape',arr.shape))
        elif arr.ndim==3:
            dims=list(arr.shape)
            ta=[i for i,n in enumerate(dims) if n==len(times)]
            if len(ta)!=1: raise RuntimeError((region,'time axis',arr.shape,ta))
            ta=ta[0]
            ra=[i for i,n in enumerate(dims) if n==100 and i!=ta]
            if len(ra)==1: ra=ra[0]
            else: ra=[i for i in range(3) if i!=ta][-1]
            sa=[i for i in range(3) if i not in (ta,ra)][0]
            x=np.moveaxis(arr,(sa,ta,ra),(0,1,2))
            sr=np.nanmean(x,axis=2)
        else: raise RuntimeError((region,'decoder ndim',arr.shape))
        try:
            sadd=np.atleast_1d(getattr(getattr(nb,'sessionaddress'),region)).ravel()
            if len(sadd)==sr.shape[0]: valid=np.array([nonempty_cell(z) for z in sadd],bool)
            else: valid=np.any(np.isfinite(sr),axis=1)
        except Exception:
            valid=np.any(np.isfinite(sr),axis=1)
        sr=sr[valid]
        for si,row in enumerate(sr):
            rec={'region':region,'decoder_session_index':si}
            for ep,(lo,hi) in EPOCHS.items():
                mask=(times>=lo-1e-9)&(times<=hi+1e-9)
                rec[f'decoder_{ep}']=float(np.nanmean(row[mask]))
            rows.append(rec)
    sess=pd.DataFrame(rows)
    reg=sess.groupby('region',as_index=False).agg(**{f'decoder_{ep}':(f'decoder_{ep}','mean') for ep in EPOCHS},n_decoder_sessions=('decoder_session_index','count'))
    return sess,reg


def unwrap(x):
    while isinstance(x,(list,tuple)) and len(x)==1: x=x[0]
    if isinstance(x,np.ndarray) and x.size==1: return unwrap(x.item())
    return x


def to_strings(x):
    x=unwrap(x)
    if isinstance(x,str): return np.asarray([x],dtype=object)
    a=np.atleast_1d(x).ravel(); out=[]
    for z in a:
        z=unwrap(z)
        if isinstance(z,str): out.append(z)
        elif isinstance(z,bytes): out.append(z.decode())
        elif isinstance(z,np.ndarray) and z.dtype.kind in 'US': out.append(''.join(z.astype(str).ravel().tolist()))
        else: out.append(str(z))
    return np.asarray(out,dtype=object)


def scalar_string(x):
    x=unwrap(x)
    if isinstance(x,str): return x
    if isinstance(x,bytes): return x.decode()
    if isinstance(x,np.ndarray) and x.dtype.kind in 'US': return ''.join(x.astype(str).ravel().tolist())
    a=np.asarray(x)
    if a.size==1: return str(a.item())
    return str(x)


def fld(o,name):
    if isinstance(o,dict): return o[name]
    return getattr(o,name)


def opto_metrics(path):
    d=mat73.loadmat(str(path))
    objs=d['optomat']
    if isinstance(objs,dict):
        # mat73 can transpose a MATLAB struct array into a dict of field lists
        n=max(len(v) if isinstance(v,list) else 1 for v in objs.values())
        oo=[]
        for i in range(n):
            z={}
            for k,v in objs.items():
                if isinstance(v,list) and len(v)==n: z[k]=v[i]
                else: z[k]=v
            oo.append(z)
        objs=oo
    objs=np.atleast_1d(objs).ravel().tolist()
    rows=[]
    print('opto objects',len(objs),type(objs[0]) if objs else None,flush=True)
    for si,o in enumerate(objs):
        region=scalar_string(fld(o,'fiber_location')).strip()
        if region not in REGIONS: continue
        trial=np.asarray(fld(o,'trial_type')).astype(float).ravel()
        lick=np.asarray(fld(o,'lick_flag')).astype(bool).ravel()
        early=np.asarray(fld(o,'early_lick')).astype(bool).ravel()
        wins=to_strings(fld(o,'opto_window'))
        n=min(len(trial),len(lick),len(early),len(wins)); trial=trial[:n]; lick=lick[:n]; early=early[:n]; wins=wins[:n]
        completed=~early
        base=(completed)&(trial==1)&(wins=='nolight')
        if base.sum()==0: continue
        base_p=lick[base].mean()*100.0
        for ep in EPOCHS:
            cur=(completed)&(trial==1)&(wins==ep)
            if cur.sum()==0: continue
            delta=lick[cur].mean()*100.0-base_p
            rows.append({'region':region,'epoch':ep,'opto_session_index':si,'n_light_trials':int(cur.sum()),'n_nolight_trials':int(base.sum()),'delta_lick_pp':float(delta),'causal_impairment_pp':float(-delta)})
    sess=pd.DataFrame(rows)
    if sess.empty: raise RuntimeError('no optogenetic rows parsed')
    reg=sess.groupby(['region','epoch'],as_index=False).agg(causal_impairment_pp=('causal_impairment_pp','mean'),delta_lick_pp=('delta_lick_pp','mean'),n_opto_sessions=('opto_session_index','nunique'),n_light_trials=('n_light_trials','sum'),n_nolight_trials=('n_nolight_trials','sum'))
    return sess,reg


def within_epoch_rank_stat(tab):
    xs=[]; ys=[]
    for ep in EPOCHS:
        q=tab[tab.epoch==ep].set_index('region').loc[REGIONS]
        rx=rankdata(q.decoder_accuracy.to_numpy(float)); ry=rankdata(q.causal_impairment_pp.to_numpy(float))
        rx=rx-rx.mean(); ry=ry-ry.mean(); xs.append(rx); ys.append(ry)
    x=np.concatenate(xs); y=np.concatenate(ys)
    den=np.sqrt((x@x)*(y@y))
    return float((x@y)/den) if den>0 else np.nan


def exact_common_region_permutation(tab):
    obs=within_epoch_rank_stat(tab); vals=[]
    base=tab.copy()
    for perm in itertools.permutations(REGIONS):
        pmap=dict(zip(REGIONS,perm))
        causal=base[['region','epoch','causal_impairment_pp']].copy(); causal['region']=causal['region'].map(pmap)
        t=base.drop(columns=['causal_impairment_pp']).merge(causal,on=['region','epoch'],how='left')
        vals.append(within_epoch_rank_stat(t))
    vals=np.asarray(vals,float); p=float(np.mean(vals>=obs-1e-12))
    return obs,p,vals


def exact_spearman(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); obs=float(spearmanr(x,y).statistic); vals=[]
    for perm in itertools.permutations(range(len(y))): vals.append(float(spearmanr(x,y[list(perm)]).statistic))
    vals=np.asarray(vals); return obs,float(np.mean(vals>=obs-1e-12))


def loo_epoch_fixed(tab):
    df=tab.copy().reset_index(drop=True); y=df.causal_impairment_pp.to_numpy(float); pred=np.full(len(df),np.nan)
    ep=pd.get_dummies(df.epoch,drop_first=True,dtype=float)
    dec=(df.decoder_accuracy.to_numpy(float)-df.decoder_accuracy.mean())/(df.decoder_accuracy.std(ddof=0)+1e-12)
    X=np.column_stack([np.ones(len(df)),dec,ep.to_numpy(float)])
    for r in REGIONS:
        te=(df.region==r).to_numpy(); tr=~te
        b=np.linalg.lstsq(X[tr],y[tr],rcond=None)[0]; pred[te]=X[te]@b
    rho=float(spearmanr(pred,y).statistic); p=float(spearmanr(pred,y).pvalue)
    mse=float(np.mean((pred-y)**2)); var=float(np.mean((y-y.mean())**2)); r2=float(1-mse/var) if var>0 else np.nan
    return {'rho':rho,'p':p,'r2':r2,'mse':mse,'predictions':pred.tolist()}


def cluster_bootstrap(tab,n=20000):
    rng=np.random.default_rng(SEED); vals=[]
    for _ in range(n):
        idx=rng.integers(0,len(REGIONS),len(REGIONS)); sampled=[]
        for j,k in enumerate(idx):
            r=REGIONS[k]; q=tab[tab.region==r].copy(); q['region']=f'{r}_{j}'; sampled.append(q)
        s=pd.concat(sampled,ignore_index=True); xs=[];ys=[]
        for ep in EPOCHS:
            q=s[s.epoch==ep]
            rx=rankdata(q.decoder_accuracy); ry=rankdata(q.causal_impairment_pp); rx-=rx.mean(); ry-=ry.mean(); xs.append(rx);ys.append(ry)
        x=np.concatenate(xs); y=np.concatenate(ys); den=np.sqrt((x@x)*(y@y))
        if den>0: vals.append(float((x@y)/den))
    return [float(np.percentile(vals,2.5)),float(np.percentile(vals,97.5))],len(vals)


def main():
    paths=extract_targets()
    dec_sess,dec_reg=decoder_metrics(paths['decoder']); opt_sess,opt_reg=opto_metrics(paths['opto'])
    rows=[]
    for _,r in dec_reg.iterrows():
        region=r.region
        for ep in EPOCHS:
            q=opt_reg[(opt_reg.region==region)&(opt_reg.epoch==ep)]
            if q.empty: continue
            z=q.iloc[0]
            rows.append({'region':region,'epoch':ep,'decoder_accuracy':float(r[f'decoder_{ep}']),'causal_impairment_pp':float(z.causal_impairment_pp),'delta_lick_pp':float(z.delta_lick_pp),'n_decoder_sessions':int(r.n_decoder_sessions),'n_opto_sessions':int(z.n_opto_sessions),'n_light_trials':int(z.n_light_trials),'n_nolight_trials':int(z.n_nolight_trials)})
    tab=pd.DataFrame(rows)
    if len(tab)!=15: raise RuntimeError(('expected 15 region-epoch cells',len(tab),tab))
    primary,pp,perms=exact_common_region_permutation(tab)
    ci,nboot=cluster_bootstrap(tab)
    epoch_results={}
    for ep in EPOCHS:
        q=tab[tab.epoch==ep].set_index('region').loc[REGIONS]
        rho,p=exact_spearman(q.decoder_accuracy,q.causal_impairment_pp); epoch_results[ep]={'spearman_rho':rho,'exact_one_sided_p':p}
    rg=tab.groupby('region',as_index=False).agg(decoder_accuracy=('decoder_accuracy','mean'),causal_impairment_pp=('causal_impairment_pp','mean')).set_index('region').loc[REGIONS]
    rr,rp=exact_spearman(rg.decoder_accuracy,rg.causal_impairment_pp)
    loo=loo_epoch_fixed(tab)
    result={
      'dataset':'Ghaderi et al. Nature Communications 2026','doi':'10.1038/s41467-026-73622-y','zenodo_record':'17973874',
      'common_regions':REGIONS,'n_region_epoch_cells':len(tab),'primary_statistic':'within-epoch rank correlation with common region-label permutation across epochs',
      'primary_rho':primary,'primary_exact_one_sided_p':pp,'primary_cluster_bootstrap_95_ci':ci,'primary_cluster_bootstrap_valid_draws':nboot,
      'success_criterion':'positive primary statistic and exact one-sided p<0.05','success':bool(primary>0 and pp<0.05),
      'epoch_specific':epoch_results,'region_average_spearman_rho':rr,'region_average_exact_one_sided_p':rp,
      'leave_one_region_out_epoch_fixed':{k:v for k,v in loo.items() if k!='predictions'},
      'file_sha256':{k:sha256(v) for k,v in paths.items()},'source_archive_url':URL,
      'important_scope_note':'Recording and inactivation experiments are from the same published behavioral system but are not the same trials; numeric analysis was prospectively specified before extracting processed values, but the publication qualitative conclusions were known.'
    }
    tab['loo_prediction']=loo['predictions']
    dec_sess.to_csv(OUT/'decoder_session_metrics.csv',index=False); opt_sess.to_csv(OUT/'opto_session_metrics.csv',index=False); tab.to_csv(OUT/'region_epoch_table.csv',index=False)
    (OUT/'primary_result.json').write_text(json.dumps(result,indent=2,sort_keys=True))
    (OUT/'RESULT.md').write_text('# Ghaderi et al. matched-system validation\n\n'+('**PASS**' if result['success'] else '**FAIL**')+f"\n\nPrimary within-epoch rank association = {primary:.6f}; exact one-sided P = {pp:.6g}; cluster-bootstrap 95% CI = [{ci[0]:.4f}, {ci[1]:.4f}].\n\nRegion-average rho = {rr:.6f}; exact one-sided P = {rp:.6g}.\n\nLOO epoch-fixed rho = {loo['rho']:.6f}; R2 = {loo['r2']:.6f}.\n")
    print(json.dumps(result,indent=2),flush=True); print(tab.to_string(index=False),flush=True)

if __name__=='__main__': main()
