"""Outcome-blind inventory of source-exported optogenetic target labels for DANDI 000009.
Reads no response values and no spikes.
"""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import json,requests,h5py,numpy as np,pandas as pd
from remfile import File as RemFile
DS='000009';VER='0.220126.1903';OUT=Path('analysis/guo2017_target_inventory');OUT.mkdir(parents=True,exist_ok=True)

def assets():
 u=f'https://api.dandiarchive.org/api/dandisets/{DS}/versions/{VER}/assets/?page_size=200';o=[]
 while u:
  j=requests.get(u,timeout=60).json();o+=j['results'];u=j.get('next')
 return o

def dec(v):
 if isinstance(v,bytes):return v.decode(errors='replace')
 a=np.asarray(v)
 if a.size==1:return dec(a.reshape(-1)[0])
 return str(v)

def one(a):
 m=requests.get(f"https://api.dandiarchive.org/api/assets/{a['asset_id']}/",timeout=60).json();cu=m.get('contentUrl') or []
 if isinstance(cu,str):cu=[cu]
 u=next((q for q in cu if 'amazonaws.com' in q),cu[0] if cu else None);r={'asset_id':a['asset_id'],'path':a['path']}
 if not u:return r
 try:
  with h5py.File(RemFile(u),'r') as f:
   r['subject']=dec(f['general/subject/subject_id'][()]) if 'general/subject/subject_id' in f else None
   if 'general/optogenetics' in f:
    sites=[]
    for name in f['general/optogenetics'].keys():
     g=f['general/optogenetics'][name];sites.append({'name':name,'location':dec(g['location'][()]) if 'location' in g else None,'description':dec(g['description'][()]) if 'description' in g else None})
    r['sites']=sites
   if 'general/extracellular_ephys/electrodes/location' in f:
    vals=np.asarray(f['general/extracellular_ephys/electrodes/location']).reshape(-1);r['recording_location']=dec(vals[0]) if len(vals) else None
 except Exception as e:r['error']=repr(e)
 return r
A=assets();rows=[]
with ThreadPoolExecutor(max_workers=12) as ex:
 for f in as_completed([ex.submit(one,a) for a in A]):rows.append(f.result())
(OUT/'targets.json').write_text(json.dumps(rows,indent=2))
flat=[]
for r in rows:
 for s in r.get('sites',[]) or []:flat.append({'asset_id':r['asset_id'],'path':r['path'],'subject':r.get('subject'),'recording_location':r.get('recording_location'),**s})
pd.DataFrame(flat).to_csv(OUT/'target_sites.csv',index=False)
summary={'n_assets':len(rows),'n_assets_with_optogenetic_site':sum(bool(r.get('sites')) for r in rows),'site_name_counts':pd.Series([s['name'] for r in rows for s in (r.get('sites') or [])]).value_counts().to_dict(),'note':'No behavioral responses or spikes were read.'}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
