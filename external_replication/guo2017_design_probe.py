"""Design-only probe for DANDI 000009.

Reads only session/subject/region metadata, trial timing columns, photostimulation
metadata and counts. It never reads the `response` trial column or spike times.
"""
from pathlib import Path
import json, requests, h5py, numpy as np
from remfile import File as RemFile

DS='000009'; VER='0.220126.1903'; OUT=Path('analysis/guo2017_design_probe'); OUT.mkdir(parents=True,exist_ok=True)
api=f'https://api.dandiarchive.org/api/dandisets/{DS}/versions/{VER}/assets/?page_size=200'
j=requests.get(api,timeout=60).json(); assets=j['results']
while j.get('next'):
    j=requests.get(j['next'],timeout=60).json(); assets+=j['results']

def decode_arr(ds):
    x=np.asarray(ds)
    out=[]
    for v in x.reshape(-1):
        if isinstance(v,bytes):out.append(v.decode(errors='replace'))
        else:out.append(str(v))
    return out

rows=[]
for a in assets[:20]:
    meta=requests.get(f"https://api.dandiarchive.org/api/assets/{a['asset_id']}/",timeout=60).json()
    cu=meta.get('contentUrl') or []
    if isinstance(cu,str):cu=[cu]
    url=next((u for u in cu if 'amazonaws.com' in u),cu[0] if cu else None)
    if not url:continue
    r={'path':a['path'],'asset_id':a['asset_id']}
    try:
      with h5py.File(RemFile(url),'r') as f:
        def scalar(path):
            if path not in f:return None
            x=f[path][()]
            if isinstance(x,bytes):return x.decode(errors='replace')
            if np.asarray(x).size==1:return str(np.asarray(x).reshape(-1)[0])
            return str(x)
        r['subject_id']=scalar('general/subject/subject_id')
        r['session_description']=scalar('session_description')
        if 'general/extracellular_ephys/electrodes/location' in f:
            r['recording_locations']=sorted(set(decode_arr(f['general/extracellular_ephys/electrodes/location'])))
        t=f['intervals/trials']
        for k in ['photo_stim_period','photo_stim_type']:
            if k in t:r[k]=sorted(set(decode_arr(t[k])))
        for k in ['photo_loc_galvo_x','photo_loc_galvo_y','photo_loc_galvo_z']:
            if k in t:
                vals=np.asarray(t[k]).astype(float).reshape(-1); vals=vals[np.isfinite(vals)]
                r[k+'_unique']=sorted(set(np.round(vals,6).tolist()))[:20]
        if 'stim_present' in t:
            vals=np.asarray(t['stim_present']).reshape(-1); r['n_trials']=int(len(vals)); r['n_stim']=int(np.sum(vals.astype(bool))); r['n_control']=int(len(vals)-np.sum(vals.astype(bool)))
        for k in ['pole_in_time','pole_out_time','cue_start_time','cue_end_time']:
            if k in t:
                vals=np.asarray(t[k]).astype(float).reshape(-1); r[k+'_finite_n']=int(np.isfinite(vals).sum())
    except Exception as e:r['error']=repr(e)
    rows.append(r)
(OUT/'design_probe.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
