import json, urllib.request, pathlib, sys, os
import scipy.io as sio
import numpy as np

OUT=pathlib.Path('analysis/pinto2019_source');OUT.mkdir(parents=True,exist_ok=True)
base='https://api.data.mendeley.com'
urls=[
 f'{base}/datasets/d2dkk9647b/files?version=2',
 f'{base}/datasets/publics/d2dkk9647b/files?version=2'
]
files=None;used=None
for u in urls:
    try:
        req=urllib.request.Request(u,headers={'User-Agent':'MBrain-pinto/1.0','Accept':'application/json'})
        with urllib.request.urlopen(req,timeout=60) as r: files=json.load(r);used=u;break
    except Exception as e: print('FAIL',u,repr(e))
if files is None: raise RuntimeError('No public file listing')
print('USED',used);print(json.dumps(files,indent=2)[:20000])
(OUT/'files.json').write_text(json.dumps(files,indent=2))
if isinstance(files,dict): vals=files.get('files') or files.get('results') or files.get('items') or []
else: vals=files
cand=[x for x in vals if 'Pinto2019_figSourceData.mat' in str(x.get('filename',''))]
if not cand: raise RuntimeError('source MAT not found in listing')
f=cand[0]; dl=(f.get('content_details') or {}).get('download_url') or f.get('download_url')
if not dl:
    fid=f.get('id');u=f'{base}/datasets/d2dkk9647b/files/{fid}/file_downloaded?version=2'
    req=urllib.request.Request(u,headers={'User-Agent':'MBrain-pinto/1.0'});resp=urllib.request.urlopen(req,timeout=60);dl=resp.geturl();resp.close()
print('DOWNLOAD',dl)
fp=OUT/'Pinto2019_figSourceData.mat';urllib.request.urlretrieve(dl,fp)
print('SIZE',fp.stat().st_size)
# inventory without loading all numeric arrays into JSON
M=sio.loadmat(fp,squeeze_me=True,struct_as_record=False)

def describe(x,depth=0,maxdepth=4):
    if depth>maxdepth:return {'type':type(x).__name__}
    if hasattr(x,'_fieldnames'):
        return {'type':'mat_struct','fields':{k:describe(getattr(x,k),depth+1,maxdepth) for k in x._fieldnames}}
    if isinstance(x,np.ndarray):
        d={'type':'ndarray','shape':list(x.shape),'dtype':str(x.dtype)}
        if x.dtype==object and x.size<=30:d['items']=[describe(z,depth+1,maxdepth) for z in x.ravel()[:30]]
        elif x.size<=30:
            try:d['values']=np.asarray(x).tolist()
            except:pass
        return d
    if isinstance(x,(str,int,float,np.number)):
        try:return {'type':type(x).__name__,'value':x.item() if hasattr(x,'item') else x}
        except:return {'type':type(x).__name__,'value':str(x)}
    return {'type':type(x).__name__}
inv={k:describe(v) for k,v in M.items() if not k.startswith('__')}
(OUT/'structure.json').write_text(json.dumps(inv,indent=2,default=str))
print(json.dumps(inv,indent=2,default=str)[:50000])
