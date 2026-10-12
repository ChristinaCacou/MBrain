from pathlib import Path
import os,json,time
import numpy as np
import h5py, requests
from scipy.stats import rankdata
OUT=Path('analysis/alm3_session'); OUT.mkdir(parents=True,exist_ok=True)
PATH=os.environ['ASSET_PATH']; URL=os.environ['ASSET_URL']; TMP=Path('alm3_replication/session.nwb')
def dec(x):
    if isinstance(x,bytes):return x.decode(errors='replace')
    if isinstance(x,np.bytes_):return bytes(x).decode(errors='replace')
    return str(x)
def hemi(f):
    vals=[]
    b='general/extracellular_ephys'
    if b in f:
        for _,g in f[b].items():
            if isinstance(g,h5py.Group) and g.attrs.get('location',None) is not None:vals.append(dec(g.attrs['location']).lower())
    p='general/extracellular_ephys/electrodes/location'
    if p in f:vals += [dec(x).lower() for x in f[p][:]]
    s=' '.join(vals); l='left alm' in s; r='right alm' in s
    return 'left' if l and not r else ('right' if r and not l else None)
def events(f,n):
    p=f'acquisition/BehavioralEvents/{n}/timestamps'; return np.asarray(f[p][:],float) if p in f else np.array([])
def to_trials(ev,st,sp):
    o=np.full(len(st),np.nan)
    for x in ev:
        i=np.searchsorted(st,x,side='right')-1
        if 0<=i<len(st) and x<=sp[i]+1e-6 and np.isnan(o[i]):o[i]=x
    return o
def lasers(f):
    eps=[]; base='stimulus/presentation'
    if base in f:
        for k in f[base].keys():
            if not k.endswith('_laser_power'):continue
            hh='left' if k.lower().startswith('left') else ('right' if k.lower().startswith('right') else None)
            g=f[f'{base}/{k}']
            if 'data' not in g or 'timestamps' not in g:continue
            p=np.nan_to_num(np.asarray(g['data'][:],float)); ts=np.asarray(g['timestamps'][:],float)
            if len(p)<2 or len(p)!=len(ts) or np.max(p)<=0:continue
            on=p>0.15*np.max(p); d=np.diff(on.astype(np.int8)); ss=list(np.where(d==1)[0]+1); ee=list(np.where(d==-1)[0]+1)
            if on[0]:ss=[0]+ss
            if on[-1]:ee=ee+[len(on)-1]
            rr=[]
            for a,b in zip(ss,ee):
                if b>a and np.isfinite(ts[a]) and np.isfinite(ts[b]) and ts[b]>ts[a]:rr.append([float(ts[a]),float(ts[b]),[(a,b)]])
            mm=[]
            for r in rr:
                if mm and r[0]-mm[-1][1]<0.03:mm[-1][1]=r[1];mm[-1][2]+=r[2]
                else:mm.append(r)
            for a,b,segs in mm:
                vals=np.concatenate([p[s:e+1] for s,e in segs]);eps.append({'start':a,'stop':b,'hemi':hh,'power':float(np.max(vals))})
    if not eps:
        ss=events(f,'photostim_start_time');ee=events(f,'photostim_stop_time');sites=[str(k).lower() for k in f['general/optogenetics'].keys()] if 'general/optogenetics' in f else []
        hh='left' if any('left alm' in x for x in sites) else ('right' if any('right alm' in x for x in sites) else None)
        for a,b in zip(ss,ee):
            if np.isfinite(a) and np.isfinite(b) and b>a:eps.append({'start':float(a),'stop':float(b),'hemi':hh,'power':np.nan})
    return eps
def auc(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float);n1=len(a);n2=len(b)
    if not n1 or not n2:return np.nan
    r=rankdata(np.r_[a,b],method='average');u=np.sum(r[:n1])-n1*(n1+1)/2;return float(u/(n1*n2))
def counts(sp,lo,hi):return np.array([np.searchsorted(sp,b,side='left')-np.searchsorted(sp,a,side='left') for a,b in zip(lo,hi)],float)
def process():
    with h5py.File(TMP,'r') as f:
        if 'intervals/trials' not in f or 'units' not in f:return {'path':PATH,'excluded':'missing trials/units'}
        t=f['intervals/trials'];u=f['units'];st=np.asarray(t['start_time'][:],float);sp=np.asarray(t['stop_time'][:],float)
        instr=np.array([dec(x).strip().lower() for x in t['trial_instruction'][:]],object);out=np.array([dec(x).strip().lower() for x in t['outcome'][:]],object);early=np.array([dec(x).strip().lower() for x in t['early_lick'][:]],object)
        go=to_trials(events(f,'go'),st,sp);hh=hemi(f)
        if hh is None:return {'path':PATH,'subject':PATH.split('/')[0],'excluded':'recorded hemisphere not identified'}
        q=np.array([dec(x).strip().lower() for x in u['quality'][:]],object) if 'quality' in u else np.array(['good']*len(u['id']),object);good=np.flatnonzero(q=='good')
        if len(good)<10:return {'path':PATH,'subject':PATH.split('/')[0],'excluded':'<10 good units','n_good':int(len(good))}
        eps=lasers(f);match=[e for e in eps if e['hemi'] in (None,hh)]
        if not match:return {'path':PATH,'subject':PATH.split('/')[0],'excluded':'no same-hemisphere photostim epochs','recorded_hemi':hh}
        any_stim=np.zeros(len(st),bool);late=np.zeros(len(st),bool);powers=[]
        for i in range(len(st)):
            if not np.isfinite(go[i]):continue
            w0,w1=go[i]-0.6,go[i]-0.1
            for e in eps:
                if max(0,min(e['stop'],sp[i])-max(e['start'],st[i]))>0:any_stim[i]=True
            for e in match:
                if max(0,min(e['stop'],w1)-max(e['start'],w0))>=0.25:
                    late[i]=True
                    if np.isfinite(e['power']):powers.append(e['power'])
                    break
        perf=np.isin(instr,['left','right'])&np.isfinite(go);noearly=np.array([('early' not in x) or ('no early' in x) for x in early]);hit=out=='hit';control=perf&noearly&(~any_stim);predsel=control&hit
        L=np.flatnonzero(predsel&(instr=='left'));R=np.flatnonzero(predsel&(instr=='right'))
        if len(L)<10 or len(R)<10:return {'path':PATH,'subject':PATH.split('/')[0],'excluded':'insufficient successful control trials','nL':int(len(L)),'nR':int(len(R))}
        allsp=np.asarray(u['spike_times'][:],float);ix=np.asarray(u['spike_times_index'][:],int);starts=np.r_[0,ix[:-1]];d=[]
        for j in good:
            s=allsp[starts[j]:ix[j]];a=auc(counts(s,go[L]-0.6,go[L]-0.1),counts(s,go[R]-0.6,go[R]-0.1))
            if np.isfinite(a):d.append(2*abs(a-0.5))
        predictor=float(np.mean(d))
        stim=perf&noearly&late
        if int(stim.sum())<10:return {'path':PATH,'subject':PATH.split('/')[0],'excluded':'<10 matched late-delay photostim trials','n_stim':int(stim.sum()),'predictor':predictor}
        imp=[];detail={}
        for lab in ['left','right']:
            ci=np.flatnonzero(control&(instr==lab));si=np.flatnonzero(stim&(instr==lab))
            if len(ci)==0 or len(si)<3:continue
            ca=float(np.mean(hit[ci]));sa=float(np.mean(hit[si]));imp.append(ca-sa);detail[lab]={'n_control':len(ci),'n_stim':len(si),'control_hit_rate':ca,'stim_hit_rate':sa,'impairment':ca-sa}
        if not imp:return {'path':PATH,'subject':PATH.split('/')[0],'excluded':'no instruction has >=3 matched stim trials'}
        return {'path':PATH,'subject':PATH.split('/')[0],'recorded_hemi':hh,'n_good_units':int(len(good)),'n_control_hit_left':int(len(L)),'n_control_hit_right':int(len(R)),'n_late_stim':int(stim.sum()),'predictor':predictor,'impairment':float(np.mean(imp)),'detail':detail,'median_laser_power':float(np.median(powers)) if powers else None}
def download():
    for a in range(4):
        try:
            with requests.get(URL,stream=True,timeout=(30,180)) as r:
                r.raise_for_status()
                with open(TMP,'wb') as f:
                    for c in r.iter_content(8*1024*1024):
                        if c:f.write(c)
            return
        except Exception:
            if a==3:raise
            time.sleep(5*(a+1))
try:download();res=process()
except Exception as e:res={'path':PATH,'subject':PATH.split('/')[0],'error':repr(e)}
finally:
    if TMP.exists():TMP.unlink()
name=PATH.replace('/','__').replace('+','_');(OUT/f'{name}.json').write_text(json.dumps(res,indent=2,default=str));print(json.dumps(res,indent=2,default=str))
