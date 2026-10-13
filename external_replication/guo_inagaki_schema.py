import json, requests, os, re
from pathlib import Path

OUT=Path('analysis/guo_inagaki_schema'); OUT.mkdir(parents=True, exist_ok=True)
BASE='https://api.dandiarchive.org/api/dandisets/000009/versions/0.220126.1903/assets/'
assets=[]; url=BASE
while url:
    j=requests.get(url,timeout=60).json(); assets += j.get('results',[]); url=j.get('next')
rows=[]
for a in assets:
    p=a.get('path','');
    rows.append({'asset_id':a.get('asset_id') or a.get('identifier'),'path':p,'size':a.get('size')})
(OUT/'assets.json').write_text(json.dumps(rows,indent=2))
print('n assets',len(rows))
print('\n'.join(f"{r['size']}\t{r['path']}\t{r['asset_id']}" for r in rows[:200]))
