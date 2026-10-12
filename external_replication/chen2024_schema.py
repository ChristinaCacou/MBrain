import json, os, tempfile, urllib.request
from pathlib import Path
import h5py

BASE='https://api.dandiarchive.org/api'
DANDISET='000363'; VERSION='0.230822.0128'
OUT=Path('analysis/chen2024_schema'); OUT.mkdir(parents=True, exist_ok=True)

def get_json(url):
    req=urllib.request.Request(url,headers={'User-Agent':'MBrain-schema/1.0','Accept':'application/json'})
    with urllib.request.urlopen(req,timeout=120) as r:return json.load(r)

url=f'{BASE}/dandisets/{DANDISET}/versions/{VERSION}/assets/?page_size=100'; A=[]
while url:
    d=get_json(url); A+=d.get('results',[]); url=d.get('next')
A=[a for a in A if a.get('path','').endswith('.nwb') and '+ogen' in a.get('path','')]
A=sorted(A,key=lambda a:a.get('size',10**18))[:3]
res=[]
for a in A:
    aid=a.get('asset_id') or a.get('id'); path=a['path'];
    if aid and aid.startswith('dandiasset:'): aid=aid.split(':',1)[1]
    fp=OUT/(path.split('/')[-1])
    urllib.request.urlretrieve(f'{BASE}/assets/{aid}/download/',fp)
    R={'asset_id':aid,'path':path,'size':a.get('size'),'groups':[]}
    with h5py.File(fp,'r') as f:
        R['groups']=list(f.keys())
        if 'units' in f:
            R['units_columns']=list(f['units'].keys())
            for k in ['location','quality','unit_name','electrodes','waveform_mean','spike_times_index']:
                if k in f['units']:
                    R[f'units_{k}_shape']=list(f['units'][k].shape)
        if 'general' in f and 'extracellular_ephys' in f['general']:
            ge=f['general/extracellular_ephys']
            R['electrode_groups']=list(ge.keys())
            eg={}
            for name in ge.keys():
                g=ge[name]; eg[name]={}
                for k in ['location','description','device']:
                    if k in g:
                        try:
                            v=g[k][()]
                            if isinstance(v,bytes):v=v.decode(errors='replace')
                            eg[name][k]=str(v)
                        except Exception: pass
            R['electrode_group_meta']=eg
        if 'general' in f and 'extracellular_ephys' in f['general'] and 'electrodes' in f['general/extracellular_ephys']:
            e=f['general/extracellular_ephys/electrodes']
            R['electrodes_columns']=list(e.keys())
            if 'location' in e:
                vals=e['location'][()]
                out=[]
                for v in vals[:]:
                    if isinstance(v,bytes):v=v.decode(errors='replace')
                    out.append(str(v))
                R['electrode_locations_unique']=sorted(set(out))
        if 'intervals' in f:
            R['interval_groups']=list(f['intervals'].keys())
            iv={}
            for name in f['intervals'].keys():
                g=f['intervals'][name]; iv[name]=list(g.keys())
            R['interval_columns']=iv
        if 'stimulus' in f:
            R['stimulus_groups']=list(f['stimulus'].keys())
            if 'presentation' in f['stimulus']:
                R['stimulus_presentation']=list(f['stimulus/presentation'].keys())
    res.append(R)
    fp.unlink(missing_ok=True)
(OUT/'schema.json').write_text(json.dumps(res,indent=2))
print(json.dumps(res,indent=2))
