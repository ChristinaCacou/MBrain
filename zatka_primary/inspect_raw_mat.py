import io,bisect,json,requests,zipfile,h5py,numpy as np
from pathlib import Path
META='https://api.figshare.com/v2/articles/13008038'; OUT=Path('analysis/zatka_primary/raw_schema'); OUT.mkdir(parents=True,exist_ok=True)
meta=requests.get(META,timeout=60).json(); parts=sorted([f for f in meta['files'] if f['name'].startswith('Zatka-Haas_et_al_Dataset.zip.')],key=lambda f:f['name']); urls=[x['download_url'] for x in parts]; sizes=[int(x['size']) for x in parts]; starts=[]; total=0
for z in sizes: starts.append(total); total+=z
class RC(io.RawIOBase):
 def __init__(self): self.pos=0; self.s=requests.Session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,o,w=0): self.pos=max(0,min(total,o if w==0 else self.pos+o if w==1 else total+o)); return self.pos
 def read(self,n=-1):
  if self.pos>=total:return b''
  if n<0:n=total-self.pos
  n=min(n,total-self.pos); out=[]
  while n:
   i=bisect.bisect_right(starts,self.pos)-1; local=self.pos-starts[i]; take=min(n,sizes[i]-local); r=self.s.get(urls[i],headers={'Range':f'bytes={local}-{local+take-1}'},timeout=120); r.raise_for_status(); b=r.content
   if len(b)!=take: raise RuntimeError((i,take,len(b),r.status_code))
   out.append(b); self.pos+=take; n-=take
  return b''.join(out)
with zipfile.ZipFile(RC()) as z:
 for member in ['OptogeneticInactivation/26CoordSet.mat','OptogeneticInactivation/Inactivation_52Coord.mat']:
  p=OUT/Path(member).name
  with z.open(member) as src, open(p,'wb') as dst:
   while True:
    b=src.read(1<<20)
    if not b:break
    dst.write(b)
  print(member,p.stat().st_size)
coord=OUT/'26CoordSet.mat'
try:
 from scipy.io import loadmat
 d=loadmat(coord,squeeze_me=True,struct_as_record=False); print('COORD KEYS',d.keys()); print({k:np.asarray(v).tolist() for k,v in d.items() if not k.startswith('__')})
except Exception as e: print('coord scipy error',repr(e))
p=OUT/'Inactivation_52Coord.mat'
print('HDF5?',h5py.is_hdf5(p))
rows=[]
with h5py.File(p,'r') as f:
 def visit(name,obj):
  x={'name':name,'type':type(obj).__name__}
  if isinstance(obj,h5py.Dataset):
   x.update(shape=list(obj.shape),dtype=str(obj.dtype));
   try:
    if obj.size<=50 and obj.dtype.kind not in 'OV': x['sample']=np.asarray(obj).reshape(-1).tolist()
   except: pass
  rows.append(x)
 f.visititems(visit)
 print('TOP',list(f.keys()))
(OUT/'schema.json').write_text(json.dumps(rows,indent=2,default=str))
for x in rows[:400]: print(x)
