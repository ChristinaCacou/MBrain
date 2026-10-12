"""Outcome-blind structural inventory for DANDI 000009.

This script deliberately does not compute or print any perturbation effect, response
probability, decoder statistic, or other numerical outcome. It only inventories asset
paths, subject identifiers encoded in metadata, and schema/column names needed to
write a frozen numerical analysis plan before outcome extraction.
"""
from pathlib import Path
import json, re
import requests
import h5py
from remfile import File as RemFile

DANDISET='000009'
VERSION='0.220126.1903'
API=f'https://api.dandiarchive.org/api/dandisets/{DANDISET}/versions/{VERSION}/assets/'
OUT=Path('analysis/guo2017_schema_inventory'); OUT.mkdir(parents=True, exist_ok=True)


def list_assets():
    url=API+'?page_size=200'
    out=[]
    while url:
        j=requests.get(url,timeout=60).json()
        out.extend(j['results']); url=j.get('next')
    return out


def walk_names(g, prefix='', depth=0, max_depth=3):
    out=[]
    if depth>max_depth:return out
    for k in list(g.keys()):
        p=f'{prefix}/{k}' if prefix else k
        out.append(p)
        try:
            if isinstance(g[k], h5py.Group): out.extend(walk_names(g[k],p,depth+1,max_depth))
        except Exception: pass
    return out

assets=list_assets()
rows=[]
for a in assets:
    path=a.get('path','')
    rows.append({'asset_id':a.get('asset_id'),'path':path,'size':a.get('size')})
(OUT/'assets.json').write_text(json.dumps(rows,indent=2))

nwb=[a for a in assets if a.get('path','').endswith('.nwb')]
sample=[]
for a in nwb[:3]:
    aid=a['asset_id']
    meta=requests.get(f'https://api.dandiarchive.org/api/assets/{aid}/',timeout=60).json()
    content=meta.get('contentUrl') or []
    if isinstance(content,str):content=[content]
    url=next((u for u in content if 'amazonaws.com' in u or 'dandiarchive' in u), content[0] if content else None)
    entry={'path':a['path'],'asset_id':aid,'content_url_present':bool(url)}
    if url:
        try:
            rf=RemFile(url)
            with h5py.File(rf,'r') as f:
                entry['top_level_keys']=list(f.keys())
                entry['schema_paths']=walk_names(f,max_depth=3)
                if 'intervals' in f and 'trials' in f['intervals']:
                    entry['trial_columns']=list(f['intervals/trials'].keys())
                if 'units' in f:
                    entry['unit_columns']=list(f['units'].keys())
                if 'stimulus' in f:
                    entry['stimulus_keys']=list(f['stimulus'].keys())
                if 'processing' in f:
                    entry['processing_keys']=list(f['processing'].keys())
        except Exception as e:
            entry['schema_error']=repr(e)
    sample.append(entry)

summary={
 'dandiset':DANDISET,'version':VERSION,'n_assets':len(assets),'n_nwb_assets':len(nwb),
 'note':'No perturbation-effect or behavioral-outcome values were extracted in this inventory.',
 'sample':sample
}
(OUT/'inventory.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
