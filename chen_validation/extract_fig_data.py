from pathlib import Path
import json, zipfile, urllib.request
import numpy as np
from scipy.io import loadmat

URL='https://zenodo.org/records/6713616/files/analyses_scripts.zip?download=1'
zip_path=Path('chen_validation/analyses_scripts.zip')
if not zip_path.exists(): urllib.request.urlretrieve(URL,zip_path)
root=Path('chen_validation/analyses_scripts')
if not root.exists() or not any(root.iterdir()):
    root.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(zip_path) as z: z.extractall(root)
out=Path('analysis/chen_fig_data'); out.mkdir(parents=True,exist_ok=True)

# MATLAB .fig files are MAT files containing nested handle-graphics structs.
def native(x):
    if isinstance(x,np.ndarray):
        if x.dtype.names:
            return [native(z) for z in x.ravel()]
        if x.dtype==object:
            return [native(z) for z in x.ravel()]
        if x.size==1: return native(x.item())
        if np.issubdtype(x.dtype,np.number): return x.astype(float).ravel().tolist()
        return x.ravel().tolist()
    if hasattr(x,'_fieldnames'):
        return {f:native(getattr(x,f)) for f in x._fieldnames}
    if isinstance(x,(np.generic,)): return x.item()
    return x

def walk(obj,path='root',acc=None):
    if acc is None: acc=[]
    if hasattr(obj,'_fieldnames'):
        d={f:getattr(obj,f) for f in obj._fieldnames}
        props=d.get('properties',None)
        ptype=d.get('type',None)
        rec={'path':path}
        if ptype is not None:
            try: rec['type']=str(np.asarray(ptype).squeeze().item())
            except: rec['type']=str(ptype)
        if props is not None and hasattr(props,'_fieldnames'):
            for key in ['XData','YData','ZData','String','Tag','DisplayName','Marker','LineStyle']:
                if hasattr(props,key):
                    val=getattr(props,key)
                    try: nv=native(val)
                    except: continue
                    if nv not in (None,[],""):
                        rec[key]=nv
        if len(rec)>1: acc.append(rec)
        for f,v in d.items():
            if f=='properties': continue
            walk(v,path+'.'+f,acc)
    elif isinstance(obj,np.ndarray):
        if obj.dtype==object or obj.dtype.names:
            for i,v in enumerate(obj.ravel()): walk(v,f'{path}[{i}]',acc)
    return acc

figs=['Figure4G_I.fig','FigureS4CF.fig','FigureS5A_E.fig','Figure1H.fig','Figure3.fig','Figure2E.fig']
summary={}
for name in figs:
    matches=list(root.rglob(name))
    if not matches:
        summary[name]={'error':'not found'}; continue
    p=matches[0]
    try:
        d=loadmat(p,squeeze_me=True,struct_as_record=False)
        acc=[]
        for k,v in d.items():
            if not k.startswith('__'): walk(v,k,acc)
        summary[name]=acc
        print('\n###',name,'objects',len(acc))
        for r in acc:
            if 'XData' in r or 'YData' in r or 'String' in r:
                print(json.dumps(r)[:2000])
    except Exception as e:
        summary[name]={'error':repr(e)}
(out/'figure_objects.json').write_text(json.dumps(summary,indent=2))
