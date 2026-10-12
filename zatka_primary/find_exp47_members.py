import io,bisect,requests,zipfile,json,re
from pathlib import Path
meta=requests.get('https://api.figshare.com/v2/articles/13008038',timeout=60).json(); parts=sorted([f for f in meta['files'] if f['name'].startswith('Zatka-Haas_et_al_Dataset.zip.')],key=lambda f:f['name']); urls=[x['download_url'] for x in parts]; sizes=[int(x['size']) for x in parts]; starts=[]; total=0
for z in sizes: starts.append(total); total+=z
class RC(io.RawIOBase):
 def __init__(self): self.pos=0; self.s=requests.Session()
 def readable(self): return True
 def seekable(self): return True
 def tell(self): return self.pos
 def seek(self,o,w=0): self.pos=max(0,min(total,o if w==0 else self.pos+o if w==1 else total+o)); return self.pos
 def read(self,n=-1):
  if self.pos>=total:return b''
  if n<0:n=total-self.pos
  n=min(n,total-self.pos); out=[]
  while n:
   i=bisect.bisect_right(starts,self.pos)-1; local=self.pos-starts[i]; take=min(n,sizes[i]-local); r=self.s.get(urls[i],headers={'Range':f'bytes={local}-{local+take-1}'},timeout=120); r.raise_for_status(); b=r.content
   if len(b)!=take: raise RuntimeError((i,take,len(b)))
   out.append(b); self.pos+=take; n-=take
  return b''.join(out)
with zipfile.ZipFile(RC()) as z:
 rows=[{'filename':x.filename,'file_size':x.file_size,'compress_size':x.compress_size} for x in z.infolist() if re.search(r'exp47|47b|inact.*47|47.*inact',x.filename,re.I)]
Path('analysis/zatka_primary').mkdir(parents=True,exist_ok=True); Path('analysis/zatka_primary/exp47_members.json').write_text(json.dumps(rows,indent=2)); print(json.dumps(rows,indent=2))
