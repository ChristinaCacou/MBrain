import csv, json, re, sys, urllib.parse, urllib.request
from pathlib import Path

BASE='https://api.dandiarchive.org/api'
DANDISET='000363'
VERSION='0.230822.0128'
OUT=Path('analysis/chen2024_external')
OUT.mkdir(parents=True, exist_ok=True)

def get_json(url):
    req=urllib.request.Request(url, headers={'User-Agent':'MBrain-external-replication/1.0','Accept':'application/json'})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

url=f'{BASE}/dandisets/{DANDISET}/versions/{VERSION}/assets/?page_size=100'
assets=[]
while url:
    d=get_json(url)
    assets.extend(d.get('results',[]))
    url=d.get('next')

rows=[]
for a in assets:
    path=a.get('path','')
    aid=a.get('asset_id') or a.get('id')
    size=a.get('size') or a.get('blob',{}).get('size')
    msub=re.search(r'(sub-[^/]+)', path)
    mses=re.search(r'(ses-[^/]+)', path)
    rows.append({
        'asset_id':aid,
        'path':path,
        'size':size,
        'subject':msub.group(1) if msub else '',
        'session':mses.group(1) if mses else '',
        'modified':a.get('date_modified','')
    })

(OUT/'assets.json').write_text(json.dumps(rows, indent=2))
with open(OUT/'assets.csv','w',newline='') as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ['asset_id','path','size','subject','session','modified'])
    w.writeheader(); w.writerows(rows)

print('N_ASSETS', len(rows))
print('TOTAL_GB', sum((r['size'] or 0) for r in rows)/1e9)
for r in rows[:25]: print(json.dumps(r))

# Pull metadata for a small deterministic sample, including NWB assets.
sample=[r for r in rows if r['path'].lower().endswith('.nwb')][:8]
md=[]
for r in sample:
    aid=r['asset_id']
    try:
        d=get_json(f'{BASE}/dandisets/{DANDISET}/versions/{VERSION}/assets/{aid}/')
        md.append(d)
    except Exception as e:
        md.append({'asset_id':aid,'path':r['path'],'error':repr(e)})
(OUT/'sample_asset_metadata.json').write_text(json.dumps(md, indent=2))
print('N_SAMPLE_METADATA', len(md))
if md:
    print('SAMPLE_KEYS', sorted(md[0].keys()))
