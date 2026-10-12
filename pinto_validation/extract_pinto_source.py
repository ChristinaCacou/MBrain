import hashlib, json, os, urllib.request, zipfile
from pathlib import Path
import h5py, numpy as np

ROOT=Path('analysis/pinto_external_source'); ROOT.mkdir(parents=True,exist_ok=True)
URLS={
 'fig2':'https://cdn.elifesciences.org/articles/70263/elife-70263-fig2-data1-v1.zip',
 'fig4':'https://cdn.elifesciences.org/articles/70263/elife-70263-fig4-data1-v1.zip',
}

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def conv(v):
 if isinstance(v,(bytes,np.bytes_)): return v.decode(errors='replace')
 a=np.asarray(v)
 if a.dtype.kind=='S': return [x.decode(errors='replace') for x in a.ravel().tolist()]
 if a.dtype.kind=='O':
  out=[]
  for x in a.ravel().tolist(): out.append(conv(x))
  return out
 if a.ndim==0: return a.item()
 return a.tolist()

def walk_h5(path):
 out={}
 with h5py.File(path,'r') as f:
  def visit(name,obj):
   if isinstance(obj,h5py.Dataset):
    rec={'shape':list(obj.shape),'dtype':str(obj.dtype)}
    try:
     v=obj[()]
     if np.asarray(v).size<=20000: rec['value']=conv(v)
    except Exception as e: rec['error']=repr(e)
    out[name]=rec
  f.visititems(visit)
 return out

meta={}
for key,url in URLS.items():
 z=ROOT/f'{key}.zip'; urllib.request.urlretrieve(url,z)
 d=ROOT/key; d.mkdir(exist_ok=True)
 with zipfile.ZipFile(z) as zz: zz.extractall(d); names=zz.namelist()
 meta[key]={'url':url,'sha256':sha(z),'zip_names':names}
 h5s=list(d.rglob('*.hdf5'))+list(d.rglob('*.h5'))
 dumps={str(p.relative_to(ROOT)):walk_h5(p) for p in h5s}
 (ROOT/f'{key}_h5dump.json').write_text(json.dumps(dumps,indent=2))
 meta[key]['h5_files']=[str(p.relative_to(ROOT)) for p in h5s]

(ROOT/'source_manifest.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2))
for key in URLS:
 p=ROOT/f'{key}_h5dump.json'
 print('\n====',key,'DATASET PATHS ====')
 d=json.loads(p.read_text())
 for fn,items in d.items():
  print('FILE',fn)
  for k,v in items.items(): print(k,v.get('shape'),v.get('dtype'))
