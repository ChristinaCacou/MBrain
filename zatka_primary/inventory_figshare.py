import json, requests
from pathlib import Path
OUT=Path('analysis/zatka_primary'); OUT.mkdir(parents=True,exist_ok=True)
url='https://api.figshare.com/v2/articles/13008038'
r=requests.get(url,timeout=60); r.raise_for_status(); data=r.json()
files=data.get('files',[])
summary=[]
for f in files:
    summary.append({k:f.get(k) for k in ['id','name','size','download_url','supplied_md5']})
(OUT/'figshare_inventory.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
