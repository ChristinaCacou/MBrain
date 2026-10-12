import json, requests
url='https://raw.githubusercontent.com/dandisets/000011/draft/.dandi/assets.json'
manifest=requests.get(url,timeout=120).json(); out=[]
for a in manifest:
    p=a.get('path','')
    if not p.endswith('_behavior+ecephys+ogen.nwb'): continue
    urls=a.get('metadata',{}).get('contentUrl',[])
    s3=next((u for u in urls if 'dandiarchive.s3.amazonaws.com' in u),None)
    if s3: out.append({'path':p,'url':s3})
print(json.dumps(out,separators=(',',':')))
