from pathlib import Path
import json, os, re, time, math
import numpy as np
import h5py, requests
from scipy.stats import rankdata, spearmanr

SEED=20261012
N_PERM=500000
N_BOOT=50000
OUT=Path('analysis/alm3_locked'); OUT.mkdir(parents=True,exist_ok=True)
CACHE=Path('alm3_replication/cache'); CACHE.mkdir(parents=True,exist_ok=True)
MANIFEST_URL='https://raw.githubusercontent.com/dandisets/000011/draft/.dandi/assets.json'


def dec(x):
    if isinstance(x, bytes): return x.decode(errors='replace')
    if isinstance(x, np.bytes_): return bytes(x).decode(errors='replace')
    return str(x)

def find_recorded_hemi(f):
    vals=[]
    base='general/extracellular_ephys'
    if base in f:
        for k,g in f[base].items():
            if isinstance(g,h5py.Group):
                loc=g.attrs.get('location',None)
                if loc is not None: vals.append(dec(loc).lower())
    if 'general/extracellular_ephys/electrodes/location' in f:
        vals += [dec(x).lower() for x in f['general/extracellular_ephys/electrodes/location'][:]]
    txt=' '.join(vals)
    left=('left alm' in txt) or ('"brain_regions": "left alm"' in txt)
    right=('right alm' in txt) or ('"brain_regions": "right alm"' in txt)
    if left and not right: return 'left'
    if right and not left: return 'right'
    return None

def read_event_times(f,name):
    p=f'acquisition/BehavioralEvents/{name}/timestamps'
    if p in f: return np.asarray(f[p][:],float)
    return np.array([],float)

def map_events_to_trials(events, starts, stops):
    out=np.full(len(starts),np.nan)
    for x in events:
        i=np.searchsorted(starts,x,side='right')-1
        if i>=0 and i<len(starts) and x<=stops[i]+1e-6 and np.isnan(out[i]): out[i]=x
    return out

def laser_epochs(f):
    eps=[]
    base='stimulus/presentation'
    if base in f:
        for key in f[base].keys():
            if not key.endswith('_laser_power'): continue
            hemi='left' if key.lower().startswith('left') else ('right' if key.lower().startswith('right') else None)
            g=f[f'{base}/{key}']
            if 'data' not in g or 'timestamps' not in g: continue
            p=np.asarray(g['data'][:],float); ts=np.asarray(g['timestamps'][:],float)
            if p.size<2 or ts.size!=p.size: continue
            p=np.nan_to_num(p,nan=0.0,posinf=0.0,neginf=0.0)
            hi=float(np.max(p)) if p.size else 0.0
            if hi<=0: continue
            on=p>0.15*hi
            d=np.diff(on.astype(np.int8))
            ss=list(np.where(d==1)[0]+1); ee=list(np.where(d==-1)[0]+1)
            if on[0]: ss=[0]+ss
            if on[-1]: ee=ee+[len(on)-1]
            raw=[]
            for s,e in zip(ss,ee):
                if e<=s or not np.isfinite(ts[s]) or not np.isfinite(ts[e]) or ts[e]<=ts[s]: continue
                raw.append([float(ts[s]),float(ts[e]),[(s,e)]])
            merged=[]
            for r in raw:
                if merged and r[0]-merged[-1][1] < 0.030:
                    merged[-1][1]=r[1]; merged[-1][2].extend(r[2])
                else: merged.append(r)
            for s_t,e_t,segs in merged:
                vals=np.concatenate([p[a:b+1] for a,b in segs])
                eps.append({'start':s_t,'stop':e_t,'hemi':hemi,'power':float(np.max(vals))})
    # fallback to behavioral photostim events if no continuous laser trace
    if not eps:
        ss=read_event_times(f,'photostim_start_time'); ee=read_event_times(f,'photostim_stop_time')
        sites=[]
        if 'general/optogenetics' in f: sites=[str(k).lower() for k in f['general/optogenetics'].keys()]
        hemi='left' if any('left alm' in x for x in sites) else ('right' if any('right alm' in x for x in sites) else None)
        for s,e in zip(ss,ee):
            if np.isfinite(s) and np.isfinite(e) and e>s: eps.append({'start':float(s),'stop':float(e),'hemi':hemi,'power':float('nan')})
    return eps

def auc_binary(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    z=np.r_[a,b]; r=rankdata(z,method='average'); n1=len(a); n2=len(b)
    if n1==0 or n2==0: return np.nan
    u=float(np.sum(r[:n1])-n1*(n1+1)/2)
    return u/(n1*n2)

def spike_lists(f, good_idx):
    st=np.asarray(f['units/spike_times'][:],float); ix=np.asarray(f['units/spike_times_index'][:],int)
    starts=np.r_[0,ix[:-1]]
    return [st[starts[i]:ix[i]] for i in good_idx]

def count_windows(spikes, lo, hi):
    return np.array([np.searchsorted(spikes,hi,side='left')-np.searchsorted(spikes,lo,side='left') for lo,hi in zip(lo,hi)],float)

def process_session(path, asset_path):
    with h5py.File(path,'r') as f:
        if 'intervals/trials' not in f or 'units' not in f: return None
        t=f['intervals/trials']; u=f['units']
        starts=np.asarray(t['start_time'][:],float); stops=np.asarray(t['stop_time'][:],float)
        instr=np.array([dec(x).strip().lower() for x in t['trial_instruction'][:]],object)
        outcome=np.array([dec(x).strip().lower() for x in t['outcome'][:]],object)
        early=np.array([dec(x).strip().lower() for x in t['early_lick'][:]],object)
        go=map_events_to_trials(read_event_times(f,'go'),starts,stops)
        rec_hemi=find_recorded_hemi(f)
        if rec_hemi is None: return {'path':asset_path,'excluded':'recorded hemisphere not identified'}
        quality=np.array([dec(x).strip().lower() for x in u['quality'][:]],object) if 'quality' in u else np.array(['good']*len(u['id']),object)
        good=np.flatnonzero(quality=='good')
        if len(good)<10: return {'path':asset_path,'subject':asset_path.split('/')[0],'excluded':'<10 good units','n_good':int(len(good))}
        eps=laser_epochs(f)
        match=[e for e in eps if e['hemi'] in (None,rec_hemi)]
        if not match: return {'path':asset_path,'subject':asset_path.split('/')[0],'excluded':'no same-hemisphere photostim epochs','recorded_hemi':rec_hemi}
        any_stim=np.zeros(len(starts),bool); late_stim=np.zeros(len(starts),bool); powers=[]
        for i in range(len(starts)):
            if not np.isfinite(go[i]): continue
            w0,w1=go[i]-0.6,go[i]-0.1
            for e in eps:
                ov=max(0.0,min(e['stop'],stops[i])-max(e['start'],starts[i]))
                if ov>0: any_stim[i]=True
            for e in match:
                ov=max(0.0,min(e['stop'],w1)-max(e['start'],w0))
                if ov>=0.25:
                    late_stim[i]=True
                    if np.isfinite(e['power']): powers.append(e['power'])
                    break
        performing=np.isin(instr,['left','right']) & np.isfinite(go)
        noearly=np.array([('early' not in x) or ('no early' in x) for x in early])
        hit=(outcome=='hit')
        control=performing & noearly & (~any_stim)
        # predictor uses only successful control trials as locked
        predsel=control & hit
        L=np.flatnonzero(predsel & (instr=='left')); R=np.flatnonzero(predsel & (instr=='right'))
        if len(L)<10 or len(R)<10: return {'path':asset_path,'subject':asset_path.split('/')[0],'excluded':'insufficient successful control trials','nL':int(len(L)),'nR':int(len(R))}
        sp=spike_lists(f,good)
        dvals=[]
        for s in sp:
            cL=count_windows(s,go[L]-0.6,go[L]-0.1); cR=count_windows(s,go[R]-0.6,go[R]-0.1)
            a=auc_binary(cL,cR)
            if np.isfinite(a): dvals.append(2*abs(a-0.5))
        if not dvals: return {'path':asset_path,'subject':asset_path.split('/')[0],'excluded':'no valid unit discriminability'}
        predictor=float(np.mean(dvals))
        stim=performing & noearly & late_stim
        if int(stim.sum())<10: return {'path':asset_path,'subject':asset_path.split('/')[0],'excluded':'<10 matched late-delay photostim trials','n_stim':int(stim.sum()),'predictor':predictor}
        impair=[]; detail={}
        for lab in ['left','right']:
            ci=np.flatnonzero(control & (instr==lab)); si=np.flatnonzero(stim & (instr==lab))
            if len(ci)==0 or len(si)<3: continue
            ca=float(np.mean(hit[ci])); sa=float(np.mean(hit[si])); imp=ca-sa
            impair.append(imp); detail[lab]={'n_control':len(ci),'n_stim':len(si),'control_hit_rate':ca,'stim_hit_rate':sa,'impairment':imp}
        if not impair: return {'path':asset_path,'subject':asset_path.split('/')[0],'excluded':'no instruction has >=3 matched stim trials'}
        return {'path':asset_path,'subject':asset_path.split('/')[0],'recorded_hemi':rec_hemi,'n_good_units':int(len(good)),
                'n_control_hit_left':int(len(L)),'n_control_hit_right':int(len(R)),'n_late_stim':int(stim.sum()),
                'predictor':predictor,'impairment':float(np.mean(impair)),'detail':detail,
                'median_laser_power':float(np.median(powers)) if powers else None}

def download(url,dest):
    for attempt in range(4):
        try:
            with requests.get(url,stream=True,timeout=(30,180)) as r:
                r.raise_for_status()
                with open(dest,'wb') as f:
                    for ch in r.iter_content(chunk_size=8*1024*1024):
                        if ch: f.write(ch)
            return
        except Exception:
            if attempt==3: raise
            time.sleep(5*(attempt+1))

def perm_p(x,y,n=N_PERM):
    x=np.asarray(x,float); y=np.asarray(y,float); obs=float(spearmanr(x,y).statistic); rng=np.random.default_rng(SEED); ge=0
    for _ in range(n):
        if spearmanr(x,rng.permutation(y)).statistic >= obs-1e-12: ge+=1
    return obs,float((ge+1)/(n+1))
def boot_ci(x,y,n=N_BOOT):
    x=np.asarray(x,float); y=np.asarray(y,float); rng=np.random.default_rng(SEED); v=[]
    for _ in range(n):
        ix=rng.integers(0,len(x),len(x)); r=spearmanr(x[ix],y[ix]).statistic
        if np.isfinite(r): v.append(r)
    return [float(np.percentile(v,2.5)),float(np.percentile(v,97.5))] if v else [None,None]
def loo(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); pred=np.full(len(y),np.nan)
    for i in range(len(y)):
        tr=np.arange(len(y))!=i; X=np.c_[np.ones(tr.sum()),x[tr]]; b=np.linalg.lstsq(X,y[tr],rcond=None)[0]; pred[i]=b[0]+b[1]*x[i]
    rho=float(spearmanr(pred,y).statistic); r2=float(1-np.mean((pred-y)**2)/np.mean((y-y.mean())**2))
    return rho,r2,pred

manifest=requests.get(MANIFEST_URL,timeout=120).json()
assets=[]
for a in manifest:
    p=a.get('path','')
    if not p.endswith('_behavior+ecephys+ogen.nwb'): continue
    urls=a.get('metadata',{}).get('contentUrl',[])
    s3=next((u for u in urls if 'dandiarchive.s3.amazonaws.com' in u),None)
    if s3: assets.append((p,s3,int(a.get('size',0))))
print('candidate assets',len(assets),flush=True)
recs=[]
for j,(p,url,size) in enumerate(assets,1):
    tmp=CACHE/'session.nwb'
    try:
        print(f'[{j}/{len(assets)}] {p} {size/1e6:.1f} MB',flush=True)
        download(url,tmp)
        rec=process_session(tmp,p); recs.append(rec or {'path':p,'excluded':'no result'})
        print(json.dumps(recs[-1],default=str),flush=True)
    except Exception as e:
        recs.append({'path':p,'error':repr(e)}); print('ERROR',repr(e),flush=True)
    finally:
        if tmp.exists(): tmp.unlink()
(OUT/'session_results.json').write_text(json.dumps(recs,indent=2))
valid=[r for r in recs if isinstance(r,dict) and 'predictor' in r and 'impairment' in r]
by={}
for r in valid: by.setdefault(r['subject'],[]).append(r)
mice=[]
for sub,rr in sorted(by.items()):
    mice.append({'subject':sub,'n_sessions':len(rr),'predictor':float(np.mean([x['predictor'] for x in rr])),'impairment':float(np.mean([x['impairment'] for x in rr])),'median_laser_power':float(np.nanmedian([x['median_laser_power'] for x in rr if x.get('median_laser_power') is not None])) if any(x.get('median_laser_power') is not None for x in rr) else None})
(OUT/'mouse_results.json').write_text(json.dumps(mice,indent=2))
if len(mice)>=4:
    x=np.array([m['predictor'] for m in mice]); y=np.array([m['impairment'] for m in mice])
    rho,p=perm_p(x,y); lr,lr2,pred=loo(x,y)
    result={'analysis_status':'LOCKED_EXTERNAL_REPLICATION','lock_commit':'89f6da46568eee087b67b9ac57abc36e82f71dfc','n_candidate_assets':len(assets),'n_eligible_sessions':len(valid),'n_mice':len(mice),'spearman_rho':rho,'one_sided_permutation_p':p,'bootstrap_95_ci':boot_ci(x,y),'loo_spearman_rho':lr,'loo_r2':lr2,'success_criterion_met':bool(rho>0 and p<0.05 and lr2>0),'subjects':[m['subject'] for m in mice]}
else:
    result={'analysis_status':'LOCKED_EXTERNAL_REPLICATION','lock_commit':'89f6da46568eee087b67b9ac57abc36e82f71dfc','n_candidate_assets':len(assets),'n_eligible_sessions':len(valid),'n_mice':len(mice),'success_criterion_met':False,'note':'fewer than four eligible mice; primary correlation not run'}
(OUT/'result.json').write_text(json.dumps(result,indent=2)); print('FINAL',json.dumps(result,indent=2),flush=True)
