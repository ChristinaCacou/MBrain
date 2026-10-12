import json, urllib.request
from pathlib import Path
import numpy as np, h5py

ASSET='0eab806c-c5c3-4d01-bd7c-15e328a7e923'
URL=f'https://api.dandiarchive.org/api/assets/{ASSET}/download/'
OUT=Path('analysis/chen2024_deep_schema'); OUT.mkdir(parents=True,exist_ok=True)
FP=OUT/'sample.nwb'; urllib.request.urlretrieve(URL,FP)

def dec(v):
    if isinstance(v,bytes): return v.decode(errors='replace')
    if isinstance(v,np.generic): return v.item()
    return v

def small_values(ds, n=20):
    try:
        x=ds[()]
        a=np.asarray(x)
        if a.size>5000:return {'shape':list(a.shape),'dtype':str(a.dtype),'sample':None}
        flat=[dec(z) for z in a.ravel()[:n]]
        u=[]
        for z in [dec(q) for q in a.ravel()]:
            s=str(z)
            if s not in u:u.append(s)
            if len(u)>=30:break
        return {'shape':list(a.shape),'dtype':str(a.dtype),'sample':flat,'unique_first30':u}
    except Exception as e:return {'error':repr(e)}

with h5py.File(FP,'r') as f:
    out={
      'identifier': dec(f['identifier'][()]) if 'identifier' in f else None,
      'session_description': dec(f['session_description'][()]) if 'session_description' in f else None,
      'session_start_time': dec(f['session_start_time'][()]) if 'session_start_time' in f else None,
      'paths':[], 'selected':{}
    }
    def visit(name,obj):
        if name.count('/')<=4: out['paths'].append({'path':name,'type':'dataset' if isinstance(obj,h5py.Dataset) else 'group','shape':list(obj.shape) if isinstance(obj,h5py.Dataset) else None})
        low=name.lower()
        if isinstance(obj,h5py.Dataset) and any(k in low for k in ['trial','lick','cue','go','stim','laser','instruction','outcome','unit_quality','classification','spike_times','electrode_group']):
            if obj.size<=5000 or any(k in low for k in ['trial_instruction','outcome','photostim','unit_quality','classification']):
                out['selected'][name]=small_values(obj)
    f.visititems(visit)
    if 'general/subject' in f:
        out['subject']={}
        for k,v in f['general/subject'].items():
            if isinstance(v,h5py.Dataset) and v.size<50:out['subject'][k]=small_values(v)
    if 'general/optogenetics' in f:
        out['optogenetics']={}
        for k,g in f['general/optogenetics'].items():
            out['optogenetics'][k]={}
            for kk,v in g.items():
                if isinstance(v,h5py.Dataset) and v.size<1000:out['optogenetics'][k][kk]=small_values(v)
    if 'intervals/trials' in f:
        out['trials']={}
        for k,v in f['intervals/trials'].items():
            if isinstance(v,h5py.Dataset):out['trials'][k]=small_values(v,30)
    if 'acquisition' in f:
        out['acquisition_groups']=list(f['acquisition'].keys())
    if 'processing' in f:
        out['processing_groups']=list(f['processing'].keys())
(OUT/'deep_schema.json').write_text(json.dumps(out,indent=2,default=str))
print(json.dumps(out,indent=2,default=str))
FP.unlink(missing_ok=True)
