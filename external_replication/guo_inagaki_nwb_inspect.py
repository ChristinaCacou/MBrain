import json, requests, h5py, os
from pathlib import Path
OUT=Path('analysis/guo_inagaki_nwb_inspect'); OUT.mkdir(parents=True,exist_ok=True)
assets=json.load(open('analysis/guo_inagaki_schema/assets.json')) if Path('analysis/guo_inagaki_schema/assets.json').exists() else None
if assets is None:
    url='https://api.dandiarchive.org/api/dandisets/000009/versions/0.220126.1903/assets/'; assets=[]
    while url:
        j=requests.get(url,timeout=60).json(); assets+=j.get('results',[]); url=j.get('next')
        assets=[{'asset_id':a.get('asset_id') or a.get('identifier'),'path':a.get('path'),'size':a.get('size')} for a in assets]
sel=[a for a in assets if 'ecephys+ogen.nwb' in a['path']][:4]
report=[]
for i,a in enumerate(sel):
    meta=requests.get(f"https://api.dandiarchive.org/api/assets/{a['asset_id']}/",timeout=60).json()
    urls=meta.get('contentUrl') or meta.get('content_url') or []
    if isinstance(urls,str): urls=[urls]
    url=next((u for u in urls if 'amazonaws.com' in u or 's3' in u), urls[0] if urls else None)
    if not url:
        # version asset endpoint
        meta2=requests.get(f"https://api.dandiarchive.org/api/dandisets/000009/versions/0.220126.1903/assets/{a['asset_id']}/",timeout=60).json()
        urls=meta2.get('contentUrl') or []
        if isinstance(urls,str): urls=[urls]
        url=next((u for u in urls if 'amazonaws.com' in u or 's3' in u), urls[0] if urls else None)
    fp=OUT/f'sample_{i}.nwb'
    r=requests.get(url,timeout=180); r.raise_for_status(); fp.write_bytes(r.content)
    item={'path':a['path'],'asset_id':a['asset_id'],'size':len(r.content),'datasets':[],'groups':[],'attrs':{}}
    with h5py.File(fp,'r') as f:
        def visit(name,obj):
            if isinstance(obj,h5py.Dataset):
                item['datasets'].append({'name':name,'shape':list(obj.shape),'dtype':str(obj.dtype)})
            else:item['groups'].append(name)
        f.visititems(visit)
        for key in ['identifier','session_description']:
            if key in f:
                try:item['attrs'][key]=f[key][()].decode() if hasattr(f[key][()],'decode') else str(f[key][()])
                except:pass
    report.append(item)
(OUT/'report.json').write_text(json.dumps(report,indent=2))
for it in report:
    print('\nFILE',it['path']);
    for d in it['datasets']:
        n=d['name'].lower()
        if any(k in n for k in ['trial','unit','spike','stim','brain','region','event','response','lick']): print(d)
