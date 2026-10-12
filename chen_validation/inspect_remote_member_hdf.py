import os, json, shutil
from pathlib import Path
import numpy as np
from remotezip import RemoteZip
import mat73
url='https://zenodo.org/records/6713616/files/datafiles_1_of_8.zip?download=1'
member='datafiles_1_of_8/data_structure_BAYLORGC12_2018_04_05.mat'
tmp=Path('chen_validation/test_raw.mat')
with RemoteZip(url) as rz, rz.open(member) as src, open(tmp,'wb') as dst:
    shutil.copyfileobj(src,dst,length=8*1024*1024)
d=mat73.loadmat(tmp)
obj=d.get('obj')
def summary(v,depth=0):
    if depth>3:return {'type':type(v).__name__}
    if isinstance(v,dict):
        return {'type':'dict','keys':list(v.keys()),'children':{k:summary(v[k],depth+1) for k in list(v.keys())[:30]}}
    if isinstance(v,list):
        return {'type':'list','len':len(v),'preview_types':[type(x).__name__ for x in v[:10]],'first':summary(v[0],depth+1) if v else None}
    if isinstance(v,np.ndarray):
        out={'type':'ndarray','shape':list(v.shape),'dtype':str(v.dtype)}
        if v.size<=30: out['preview']=v.ravel().tolist()
        return out
    if isinstance(v,(str,int,float,bool,np.integer,np.floating)) or v is None:
        return {'type':type(v).__name__,'value':v}
    return {'type':type(v).__name__,'repr':repr(v)[:500]}
rep=summary(obj)
os.makedirs('analysis/chen_hdf_schema',exist_ok=True)
open('analysis/chen_hdf_schema/schema.json','w').write(json.dumps(rep,indent=2,default=str))
print(json.dumps(rep,indent=2,default=str))
tmp.unlink()
