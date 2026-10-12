import argparse, json, math, re, os
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO
import requests

VERSION='0.220126.1907'
WIN=(-1.7,-0.9)
MIN_CTRL_PER_CHOICE=10
MIN_STIM_PER_CHOICE=5

def _auc(scores, labels):
    scores=np.asarray(scores,float); labels=np.asarray(labels,int)
    n1=int(labels.sum()); n0=int((labels==0).sum())
    if n1==0 or n0==0: return np.nan
    r=rankdata(scores,method='average')
    u=r[labels==1].sum()-n1*(n1+1)/2
    return float(u/(n1*n0))

def _inside_intervals(a,b,intervals):
    try:
        arr=np.asarray(intervals,float).reshape(-1,2)
    except Exception:
        return True
    if arr.size==0: return True
    return bool(np.any((a>=arr[:,0]-1e-8)&(b<=arr[:,1]+1e-8)))

def _parse_hemi(nwb):
    locs=[]
    for eg in nwb.electrode_groups.values():
        locs.append(str(getattr(eg,'location','')).lower())
    s=' '.join(locs)
    if 'left alm' in s: return 'left'
    if 'right alm' in s: return 'right'
    # fallback from unit electrode_group text
    if nwb.units is not None:
        try:
            u=nwb.units.to_dataframe()
            txt=' '.join(str(x).lower() for x in u.get('electrode_group',[]))
            if 'left alm' in txt: return 'left'
            if 'right alm' in txt: return 'right'
        except Exception: pass
    return None

def _event_map_to_trials(event_ts, starts, stops):
    out={}
    for t in np.asarray(event_ts,float):
        if not np.isfinite(t): continue
        idx=int(np.searchsorted(starts,t,side='right')-1)
        if idx>=0 and idx<len(starts) and starts[idx]-1e-8<=t<=stops[idx]+1e-8:
            out.setdefault(idx,[]).append(float(t))
    return out

def _series_max(ts, t0, t1):
    try:
        tt=np.asarray(ts.timestamps[:],float)
        if tt.size==0: return -np.inf
        i0=int(np.searchsorted(tt,t0,side='left')); i1=int(np.searchsorted(tt,t1,side='right'))
        i0=max(0,i0); i1=min(len(tt),max(i0+1,i1))
        x=np.asarray(ts.data[i0:i1],float).ravel()
        x=x[np.isfinite(x)]
        return float(np.max(x)) if x.size else -np.inf
    except Exception:
        return -np.inf

def _photostim_sites(nwb, starts, stops):
    if 'BehavioralEvents' not in nwb.acquisition: return {}, {}
    be=nwb.acquisition['BehavioralEvents']
    if not hasattr(be,'time_series') or 'photostim_start_time' not in be.time_series: return {}, {}
    st=np.asarray(be.time_series['photostim_start_time'].timestamps[:],float)
    if 'photostim_stop_time' in be.time_series:
        en=np.asarray(be.time_series['photostim_stop_time'].timestamps[:],float)
    else:
        en=np.full_like(st,np.nan)
    sites=list(nwb.ogen_sites.keys()) if getattr(nwb,'ogen_sites',None) is not None else []
    laser={}
    for name,ts in nwb.stimulus.items():
        if name.endswith('_laser_power'):
            laser[name[:-len('_laser_power')]]=ts
    trial_site={}; detail={}
    for j,t0 in enumerate(st):
        if not np.isfinite(t0): continue
        idx=int(np.searchsorted(starts,t0,side='right')-1)
        if idx<0 or idx>=len(starts) or not (starts[idx]-1e-8<=t0<=stops[idx]+1e-8): continue
        t1=float(en[j]) if j<len(en) and np.isfinite(en[j]) and en[j]>t0 else min(float(stops[idx]),float(t0)+0.9)
        if len(sites)==1:
            site=sites[0]
            scores={site:None}
        elif laser:
            scores={k:_series_max(ts,float(t0),float(t1)) for k,ts in laser.items()}
            site=max(scores,key=scores.get)
        else:
            site='unknown'
            scores={}
        trial_site[idx]=site
        detail[idx]={'start':float(t0),'stop':float(t1),'site':site,'trace_max':scores}
    return trial_site, detail

def analyze(nwb_path, asset_path):
    with NWBHDF5IO(str(nwb_path),'r',load_namespaces=True) as io:
        nwb=io.read()
        mouse=str(getattr(getattr(nwb,'subject',None),'subject_id','unknown'))
        trials=nwb.trials.to_dataframe().reset_index(drop=True)
        units=nwb.units.to_dataframe().reset_index(drop=True)
        starts=trials['start_time'].to_numpy(float); stops=trials['stop_time'].to_numpy(float)
        hemi=_parse_hemi(nwb)
        if 'BehavioralEvents' not in nwb.acquisition or 'go' not in nwb.acquisition['BehavioralEvents'].time_series:
            return {'asset_path':asset_path,'mouse':mouse,'eligible':False,'reason':'no_go_events'}
        go_map=_event_map_to_trials(nwb.acquisition['BehavioralEvents'].time_series['go'].timestamps[:],starts,stops)
        stim_map, stim_detail=_photostim_sites(nwb,starts,stops)
        instr=trials['trial_instruction'].astype(str).str.lower().to_numpy()
        early=trials['early_lick'].astype(str).str.lower().to_numpy()
        outcome=trials['outcome'].astype(str).str.lower().to_numpy()
        valid=np.array([(i in go_map) and instr[i] in ('left','right') and early[i]=='no early' and outcome[i] in ('hit','miss') for i in range(len(trials))])
        go=np.full(len(trials),np.nan)
        for i,v in go_map.items(): go[i]=v[0]
        control=valid & np.array([i not in stim_map for i in range(len(trials))])
        same=np.zeros(len(trials),bool); contra=np.zeros(len(trials),bool); bilateral=np.zeros(len(trials),bool)
        for i,s in stim_map.items():
            sl=str(s).lower()
            if 'bilateral' in sl: bilateral[i]=True
            elif hemi and hemi in sl: same[i]=True
            elif ('left' in sl or 'right' in sl): contra[i]=True
        same &= valid; contra &= valid; bilateral &= valid
        # Frozen outcome requires minimum control and same-hemisphere perturbation counts for each choice.
        counts={}
        diffs=[]
        for choice in ('left','right'):
            c=(instr==choice)
            ci=np.where(control & c)[0]; si=np.where(same & c)[0]
            counts[choice]={'control':int(len(ci)),'same_hemi_stim':int(len(si))}
            if len(ci)>=MIN_CTRL_PER_CHOICE and len(si)>=MIN_STIM_PER_CHOICE:
                ac=float(np.mean(outcome[ci]=='hit')); ap=float(np.mean(outcome[si]=='hit'))
                counts[choice].update({'control_accuracy':ac,'stim_accuracy':ap,'impairment':ac-ap})
                diffs.append(ac-ap)
        behavioral_impairment=float(np.mean(diffs)) if len(diffs)==2 else np.nan
        # Frozen predictor, control trials only. Unit-specific coverage is respected via obs_intervals.
        unit_scores=[]; unit_meta=[]
        for ui,row in units.iterrows():
            quality=str(row.get('quality','')).lower(); cell=str(row.get('cell_type','')).lower()
            if quality and quality!='good': continue
            if cell and cell not in ('pyr','pyramidal'): continue
            obs=row.get('obs_intervals',None)
            idx=[]; lab=[]; sc=[]
            sp=np.asarray(row['spike_times'],float)
            for ti in np.where(control)[0]:
                a=go[ti]+WIN[0]; b=go[ti]+WIN[1]
                if not np.isfinite(a) or not _inside_intervals(a,b,obs): continue
                idx.append(int(ti)); lab.append(1 if instr[ti]=='right' else 0)
                sc.append(int(np.searchsorted(sp,b,'left')-np.searchsorted(sp,a,'left')))
            lab=np.asarray(lab,int)
            if int(np.sum(lab==0))<MIN_CTRL_PER_CHOICE or int(np.sum(lab==1))<MIN_CTRL_PER_CHOICE: continue
            auc=_auc(sc,lab); sel=2*abs(auc-0.5)
            unit_scores.append(sel); unit_meta.append({'unit_index':int(ui),'auc':auc,'selectivity':sel,'n_left':int(np.sum(lab==0)),'n_right':int(np.sum(lab==1))})
        predictor=float(np.mean(unit_scores)) if unit_scores else np.nan
        eligible=bool(np.isfinite(predictor) and np.isfinite(behavioral_impairment))
        return {
          'asset_path':asset_path,'mouse':mouse,'recorded_hemisphere':hemi,'eligible':eligible,
          'reason':None if eligible else ('insufficient_predictor_or_outcome'),
          'predictor_local_choice_selectivity':predictor,'behavioral_impairment':behavioral_impairment,
          'n_eligible_units':int(len(unit_scores)),'n_valid_trials':int(valid.sum()),'n_control_trials':int(control.sum()),
          'n_same_hemi_stim_trials':int(same.sum()),'n_contra_stim_trials':int(contra.sum()),'n_bilateral_stim_trials':int(bilateral.sum()),
          'choice_counts':counts,'n_photostim_events':int(len(stim_map)),
          'stim_site_counts':{s:int(sum(1 for v in stim_map.values() if v==s)) for s in sorted(set(stim_map.values()))},
          'unit_scores':unit_meta,
          'session_identifier':str(nwb.identifier),
          'session_description':str(nwb.session_description),
          'subject_description':str(getattr(getattr(nwb,'subject',None),'description','')),
          'subject_genotype':str(getattr(getattr(nwb,'subject',None),'genotype','')),
          'analysis_window_relative_go':list(WIN)
        }

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--asset',required=True); ap.add_argument('--out',required=True); args=ap.parse_args()
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    tmp=Path('alm3_replication')/('tmp_'+re.sub(r'[^A-Za-z0-9_.-]','_',args.asset.split('/')[-1]))
    with DandiAPIClient() as client:
        ds=client.get_dandiset('000011',VERSION)
        asset=next((a for a in ds.get_assets() if a.path==args.asset),None)
        if asset is None: raise RuntimeError('asset not found: '+args.asset)
        url=asset.get_content_url(follow_redirects=1, strip_query=False)
        with requests.get(url,stream=True,timeout=300) as r:
            r.raise_for_status()
            with open(tmp,'wb') as f:
                for chunk in r.iter_content(8*1024*1024):
                    if chunk: f.write(chunk)
    try:
        res=analyze(tmp,args.asset)
        Path(args.out).write_text(json.dumps(res,indent=2,allow_nan=True))
        print(json.dumps({k:v for k,v in res.items() if k not in ('unit_scores',)},indent=2,allow_nan=True))
    finally:
        tmp.unlink(missing_ok=True)
if __name__=='__main__': main()
