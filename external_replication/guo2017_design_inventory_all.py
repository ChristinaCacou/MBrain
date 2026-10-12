"""Outcome-blind structural inventory over all DANDI 000009 NWB assets.

Reads no response column and no spike times. Used only to determine whether the
preplanned recording/perturbation design is structurally estimable before freezing
numerical outcomes.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import json,re,requests,h5py,numpy as np,pandas as pd
from remfile import File as RemFile
DS='000009'; VER='0.220126.1903'; OUT=Path('analysis/guo2017_design_inventory_all'); OUT.mkdir(parents=True,exist_ok=True)

def assets():
 u=f'https://api.dandiarchive.org/api/dandisets/{DS}/versions/{VER}/assets/?page_size=200'; out=[]
 while u:
  j=requests.get(u,timeout=60).json();out+=j['results'];u=j.get('next')
 return out

def dec(x):
 a=np.asarray(x).reshape(-1);return [v.decode(errors='replace') if isinstance(v,bytes) else str(v) for v in a]

def one(a):
 meta=requests.get(f"https://api.dandiarchive.org/api/assets/{a['asset_id']}/",timeout=60).json();cu=meta.get('contentUrl') or []
 if isinstance(cu,str):cu=[cu]
 url=next((u for u in cu if 'amazonaws.com' in u),cu[0] if cu else None)
 r={'asset_id':a['asset_id'],'path':a['path']}
 if not url:return r
 try:
  with h5py.File(RemFile(url),'r') as f:
   sid=f['general/subject/subject_id'][()];r['subject']=sid.decode() if isinstance(sid,bytes) else str(sid)
   loc=dec(f['general/extracellular_ephys/electrodes/location']) if 'general/extracellular_ephys/electrodes/location' in f else []
   r['recording_location']='|'.join(sorted(set(loc)))
   t=f['intervals/trials']
   sp=np.asarray(t['stim_present']).reshape(-1).astype(bool) if 'stim_present' in t else np.zeros(len(t['id']),bool)
   r['n_trials']=int(len(sp));r['n_stim']=int(sp.sum());r['has_delay_stim']=False
   if 'photo_stim_period' in t:
    periods=np.asarray(dec(t['photo_stim_period']));r['has_delay_stim']=bool(np.any((periods=='delay') & sp))
   for k in ['photo_loc_galvo_x','photo_loc_galvo_y','photo_loc_galvo_z']:
    if k in t:
     v=np.asarray(t[k]).reshape(-1).astype(float);v=v[sp & np.isfinite(v)] if len(v)==len(sp) else v[np.isfinite(v)]
     r[k]=','.join(map(str,sorted(set(np.round(v,4).tolist()))))
   if 'type' in t:r['trial_types']='|'.join(sorted(set(dec(t['type']))))
   r['has_units']='units' in f;r['n_units']=int(len(f['units/id'])) if 'units/id' in f else 0
 except Exception as e:r['error']=repr(e)
 return r

A=assets();rows=[]
with ThreadPoolExecutor(max_workers=12) as ex:
 fut=[ex.submit(one,a) for a in A]
 for f in as_completed(fut):rows.append(f.result())
df=pd.DataFrame(rows);df.to_csv(OUT/'session_design.csv',index=False)
stim=df[(df.n_stim.fillna(0)>0)&(df.has_delay_stim.fillna(False))].copy()
# Parse coarse recording region from source location string.
def reg(s):
 s=str(s)
 m=re.search(r'brain_region:\s*([^;|]+)',s)
 return m.group(1).strip() if m else 'unknown'
stim['recording_region']=stim.recording_location.map(reg)
summary={
 'dandiset':DS,'version':VER,'n_assets':len(df),'n_subjects':int(df.subject.nunique()),
 'n_delay_stim_sessions':int(len(stim)),'n_delay_stim_subjects':int(stim.subject.nunique()),
 'delay_stim_by_recording_region':stim.groupby('recording_region').size().to_dict(),
 'delay_stim_subjects_by_recording_region':stim.groupby('recording_region').subject.nunique().to_dict(),
 'stim_coordinate_combinations':stim.groupby(['recording_region','photo_loc_galvo_x','photo_loc_galvo_y','photo_loc_galvo_z']).size().reset_index(name='n').to_dict('records'),
 'note':'No behavioral response values or spike times were read.'
}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
