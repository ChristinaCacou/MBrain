import io, json, requests, zipfile, bisect, re
from pathlib import Path
API='https://api.figshare.com/v2/articles/13008038'
OUT=Path('analysis/zatka_primary/splitzip'); OUT.mkdir(parents=True,exist_ok=True)
meta=requests.get(API,timeout=60).json()
parts=[f for f in meta['files'] if f['name'].startswith('Zatka-Haas_et_al_Dataset.zip.')]
parts=sorted(parts,key=lambda f:f['name'])
urls=[p['download_url'] for p in parts]; sizes=[int(p['size']) for p in parts]
starts=[]; s=0
for z in sizes: starts.append(s); s+=z
class RemoteConcat(io.RawIOBase):
    def __init__(self): self.pos=0; self.session=requests.Session(); self.total=s
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self,offset,whence=io.SEEK_SET):
        if whence==io.SEEK_SET: p=offset
        elif whence==io.SEEK_CUR: p=self.pos+offset
        elif whence==io.SEEK_END: p=self.total+offset
        else: raise ValueError
        self.pos=max(0,min(self.total,p)); return self.pos
    def read(self,n=-1):
        if self.pos>=self.total:return b''
        if n is None or n<0:n=self.total-self.pos
        n=min(n,self.total-self.pos); chunks=[]; remain=n
        while remain>0:
            i=bisect.bisect_right(starts,self.pos)-1
            local=self.pos-starts[i]
            take=min(remain,sizes[i]-local)
            h={'Range':f'bytes={local}-{local+take-1}'}
            r=self.session.get(urls[i],headers=h,timeout=120); r.raise_for_status()
            data=r.content
            if len(data)!=take:
                raise RuntimeError(f'range mismatch part {i} requested {take}, got {len(data)}, status {r.status_code}')
            chunks.append(data); self.pos+=take; remain-=take
        return b''.join(chunks)
rf=RemoteConcat()
with zipfile.ZipFile(rf) as z:
    infos=z.infolist()
    rows=[]
    for x in infos:
        if re.search(r'(inactivation|26CoordSet|ctxOutlines|sessionList)',x.filename,re.I):
            rows.append({'filename':x.filename,'file_size':x.file_size,'compress_size':x.compress_size,'header_offset':x.header_offset,'compress_type':x.compress_type})
    (OUT/'matching_members.json').write_text(json.dumps(rows,indent=2))
    print('total members',len(infos),'matches',len(rows),'archive bytes',s)
    print(json.dumps(rows,indent=2))
