import json, urllib.request
from pathlib import Path
import numpy as np, h5py
from remfile import File as RemoteFile

BASE='https://api.dandiarchive.org/api'; DS='000363'; VER='0.230822.0128'
OUT=Path('analysis/chen2024_remote_inventory'); OUT.mkdir(parents=True,exist_ok=True)

def getj(url):
    req=urllib.request.Request(url,headers={'User-Agent':'MBrain-remote-inventory/1.0','Accept':'application/json'})
    with urllib.request.urlopen(req,timeout=120) as r:return json.load(r)
def dec(x):
    if isinstance(x,bytes):return x.decode(errors='replace')
    if isinstance(x,np.generic):return x.item()
    return x
def vals(ds,limit=2000):
    try:
        a=np.asarray(ds[()]).ravel()
        if len(a)>limit:a=a[:limit]
        return [dec(x) for x in a]
    except:return []
def uniq(ds,limit=2000,nmax=20):
    out=[]
    for x in vals(ds,limit):
        s=str(x)
        if s not in out:out.append(s)
        if len(out)>=nmax:break
    return out

url=f'{BASE}/dandisets/{DS}/versions/{VER}/assets/?page_size=100'; A=[]
while url:
    d=getj(url);A+=d.get('results',[]);url=d.get('next')
A=[a for a in A if a.get('path','').endswith('.nwb')]
rows=[]
for ii,a in enumerate(A):
    aid=a.get('asset_id') or a.get('id');
    if aid and aid.startswith('dandiasset:'):aid=aid.split(':',1)[1]
    p=a.get('path',''); row={'asset_id':aid,'path':p,'size':a.get('size')}
    try:
        u=f'{BASE}/assets/{aid}/download/'
        rf=RemoteFile(u)
        with h5py.File(rf,'r') as f:
            if 'general/subject/subject_id' in f:row['subject_id']=dec(f['general/subject/subject_id'][()])
            if 'session_start_time' in f:row['session_start_time']=dec(f['session_start_time'][()])
            loc=[]
            if 'general/extracellular_ephys/electrodes/location' in f:loc=uniq(f['general/extracellular_ephys/electrodes/location'],limit=100000,nmax=50)
            row['locations']=loc
            row['has_left_ALM']=any('left ALM' in x for x in loc)
            row['has_right_ALM']=any('right ALM' in x for x in loc)
            if 'units/id' in f:row['n_units']=int(f['units/id'].shape[0])
            if 'units/unit_quality' in f:row['unit_quality_values']=uniq(f['units/unit_quality'],limit=100000,nmax=20)
            if 'units/classification' in f:row['classification_values']=uniq(f['units/classification'],limit=100000,nmax=20)
            if 'intervals/trials/id' in f:row['n_trials']=int(f['intervals/trials/id'].shape[0])
            for k in ['task','task_protocol','trial_instruction','outcome','photostim_onset','photostim_duration','photostim_power']:
                q=f'intervals/trials/{k}'
                if q in f:row[k+'_values']=uniq(f[q],limit=100000,nmax=30)
            # record any path name that could encode photostim side
            names=[]
            def visit(name,obj):
                low=name.lower()
                if any(k in low for k in ['opto','photostim','laser','stim_site','stimsite']):names.append(name)
            f.visititems(visit);row['opto_paths']=names[:100]
        rf.close()
    except Exception as e:row['error']=repr(e)
    rows.append(row);print(ii+1,len(A),json.dumps(row),flush=True)
(OUT/'session_inventory.json').write_text(json.dumps(rows,indent=2,default=str))
summary={'n_assets':len(rows),'n_ok':sum('error' not in r for r in rows),'n_bilateral_ALM':sum(r.get('has_left_ALM') and r.get('has_right_ALM') for r in rows),'bilateral_paths':[r['path'] for r in rows if r.get('has_left_ALM') and r.get('has_right_ALM')], 'bilateral_subjects':sorted(set(str(r.get('subject_id')) for r in rows if r.get('has_left_ALM') and r.get('has_right_ALM')))}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2));print('SUMMARY',json.dumps(summary,indent=2),flush=True)
