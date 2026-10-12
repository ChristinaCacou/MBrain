import subprocess,zipfile,json,re
from pathlib import Path
import numpy as np
from scipy.io import loadmat
subprocess.run(['aws','s3','cp','--no-sign-request','s3://crcnsarchive/alm-3/live/data/meta_data_files.zip','meta.zip'],check=True)
with zipfile.ZipFile('meta.zip') as z:z.extractall('meta')
out=Path('analysis/alm3_meta_inspection');out.mkdir(parents=True,exist_ok=True)
def strings(x,depth=0):
    if depth>6:return []
    if x is None:return []
    if isinstance(x,(str,bytes)):return [x.decode(errors='ignore') if isinstance(x,bytes) else x]
    if isinstance(x,dict):
        r=[]
        for k,v in x.items():r+=[str(k)];r+=strings(v,depth+1)
        return r
    if hasattr(x,'_fieldnames'):
        r=[]
        for k in x._fieldnames:r+=[k];r+=strings(getattr(x,k,None),depth+1)
        return r
    if isinstance(x,np.ndarray):
        if x.dtype.kind in 'US':return [str(x)]
        r=[]
        for z in x.ravel()[:500]:r+=strings(z,depth+1)
        return r
    if isinstance(x,(list,tuple)):
        r=[]
        for z in x[:500]:r+=strings(z,depth+1)
        return r
    return []
rows=[]
for p in sorted(Path('meta').rglob('*.mat')):
    try:d=loadmat(p,squeeze_me=True,struct_as_record=False)
    except Exception as e:
        rows.append({'file':p.name,'error':repr(e)});continue
    ss=strings({k:v for k,v in d.items() if not k.startswith('__')});txt=' | '.join(s for s in ss if s and len(s)<500)
    m=re.search(r'(ANM\d+)_?(\d{8})?',p.name)
    rows.append({'file':p.name,'animal':m.group(1) if m else None,'date':m.group(2) if m else None,'callosotomy_keyword':bool(re.search(r'callos|corpus call|cc.?cut|bisection',txt,re.I)),'excerpt':txt[:5000]})
(out/'metadata_scan.json').write_text(json.dumps(rows,indent=2))
kw=[r for r in rows if r.get('callosotomy_keyword')]
(out/'callosotomy_keyword_hits.json').write_text(json.dumps(kw,indent=2))
# compact per-animal excerpts, useful even if explicit keyword is absent
by={}
for r in rows:
    a=r.get('animal');
    if a:by.setdefault(a,[]).append({'file':r['file'],'hit':r.get('callosotomy_keyword'),'excerpt':r.get('excerpt','')[:1200]})
(out/'by_animal.json').write_text(json.dumps(by,indent=2))
print(json.dumps({'n_metadata':len(rows),'n_keyword_hits':len(kw),'hit_animals':sorted(set(r.get('animal') for r in kw if r.get('animal')))},indent=2))
