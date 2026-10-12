import json, struct, urllib.request
from pathlib import Path

PARTS=[
 ('001','https://ndownloader.figshare.com/files/24786056',5033164800),
 ('002','https://ndownloader.figshare.com/files/24786080',5033164800),
 ('003','https://ndownloader.figshare.com/files/24786128',5033164800),
 ('004','https://ndownloader.figshare.com/files/24786167',4524251185),
]

def get_range(url,start,end):
    req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}','User-Agent':'Mozilla/5.0'})
    with urllib.request.urlopen(req,timeout=120) as r:
        return r.read()

last_url=PARTS[-1][1]; last_size=PARTS[-1][2]
tail_len=min(4_000_000,last_size)
tail=get_range(last_url,last_size-tail_len,last_size-1)
pos=tail.rfind(b'PK\x05\x06')
if pos<0: raise RuntimeError('EOCD not found')
eocd=tail[pos:pos+22]
(sig,disk_no,cd_disk,n_disk,n_total,cd_size,cd_off,comment_len)=struct.unpack('<4sHHHHIIH',eocd)
print('EOCD',disk_no,cd_disk,n_disk,n_total,cd_size,cd_off,comment_len)
# central directory is usually wholly on one disk; support spill forward if needed
remaining=cd_size; disk=cd_disk; off=cd_off; chunks=[]
while remaining>0:
    if disk>=len(PARTS): raise RuntimeError('central dir disk overflow')
    size=PARTS[disk][2]; take=min(remaining,size-off)
    chunks.append(get_range(PARTS[disk][1],off,off+take-1))
    remaining-=take; disk+=1; off=0
cd=b''.join(chunks)
entries=[]; i=0
while i+46<=len(cd):
    if cd[i:i+4]!=b'PK\x01\x02':
        break
    vals=struct.unpack('<4s6H3I5H2I',cd[i:i+46])
    _,vermade,verneed,flag,method,mtime,mdate,crc,csize,usize,nlen,xlen,clen,diskstart,intattr,extattr,loff=vals
    name=cd[i+46:i+46+nlen].decode('utf-8','replace')
    extra=cd[i+46+nlen:i+46+nlen+xlen]
    comment=cd[i+46+nlen+xlen:i+46+nlen+xlen+clen]
    entries.append(dict(name=name,compressed_size=csize,uncompressed_size=usize,method=method,flag=flag,disk_start=diskstart,local_offset=loff,crc32=crc))
    i+=46+nlen+xlen+clen
out=Path('analysis/zatka_primary'); out.mkdir(parents=True,exist_ok=True)
(out/'split_zip_index.json').write_text(json.dumps({'eocd':dict(disk_no=disk_no,cd_disk=cd_disk,n_total=n_total,cd_size=cd_size,cd_off=cd_off),'entries':entries},indent=2))
for e in entries:
    n=e['name'].lower()
    if any(k in n for k in ['inactiv','opto','exp47','sessionlist','coordset','ctxoutline']):
        print(json.dumps(e))
print('entries',len(entries))
