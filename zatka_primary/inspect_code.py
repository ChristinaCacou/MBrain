import zipfile, re, json, os, requests
from pathlib import Path
OUT=Path('analysis/zatka_primary/code_inspect'); OUT.mkdir(parents=True,exist_ok=True)
url='https://ndownloader.figshare.com/files/24786170'
zpath=OUT/'code.zip'
with requests.get(url,stream=True,timeout=120) as r:
    r.raise_for_status()
    with open(zpath,'wb') as f:
        for chunk in r.iter_content(1<<20):
            if chunk:f.write(chunk)
with zipfile.ZipFile(zpath) as z:
    names=z.namelist()
    (OUT/'file_list.txt').write_text('\n'.join(names))
    hits=[]
    pats=re.compile(r'(exp47|opto|laser|inactiv|choice|atlas|allen|coordinate)',re.I)
    for n in names:
        if not n.lower().endswith(('.m','.py','.txt','.md')): continue
        try: txt=z.read(n).decode('utf-8','ignore')
        except Exception: continue
        if pats.search(n) or pats.search(txt):
            ms=[]
            for i,line in enumerate(txt.splitlines(),1):
                if pats.search(line): ms.append(f'{i}: {line}')
            hits.append({'file':n,'matches':ms[:200]})
            outname=n.replace('/','__')+'.txt'
            (OUT/outname).write_text(txt)
    (OUT/'hits.json').write_text(json.dumps(hits,indent=2))
print('files',len(names),'hits',len(hits))
for h in hits:
    print('\n###',h['file'])
    print('\n'.join(h['matches'][:30]))
