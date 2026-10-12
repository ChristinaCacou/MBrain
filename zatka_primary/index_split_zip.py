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

def read_across(disk,off,n):
    chunks=[]; remain=n
    while remain>0:
        if disk>=len(PARTS): raise RuntimeError(f'read overflow disk={disk} off={off} remain={remain}')
        size=PARTS[disk][2]
        if off>=size:
            disk+=1; off=0; continue
        take=min(remain,size-off)
        chunks.append(get_range(PARTS[disk][1],off,off+take-1))
        remain-=take; disk+=1; off=0
    return b''.join(chunks)

def zip64_extra(extra, need_usize, need_csize, need_loff, need_disk):
    i=0; vals=[]
    while i+4<=len(extra):
        hid,sz=struct.unpack('<HH',extra[i:i+4]); body=extra[i+4:i+4+sz]; i+=4+sz
        if hid!=0x0001: continue
        j=0; out={}
        if need_usize: out['usize']=struct.unpack('<Q',body[j:j+8])[0]; j+=8
        if need_csize: out['csize']=struct.unpack('<Q',body[j:j+8])[0]; j+=8
        if need_loff: out['loff']=struct.unpack('<Q',body[j:j+8])[0]; j+=8
        if need_disk: out['disk']=struct.unpack('<I',body[j:j+4])[0]; j+=4
        return out
    return {}

last_url=PARTS[-1][1]; last_size=PARTS[-1][2]
tail_len=min(4_000_000,last_size)
tail_start=last_size-tail_len
tail=get_range(last_url,tail_start,last_size-1)
pos=tail.rfind(b'PK\x05\x06')
if pos<0: raise RuntimeError('EOCD not found')
eocd=tail[pos:pos+22]
(sig,disk_no,cd_disk,n_disk,n_total,cd_size32,cd_off32,comment_len)=struct.unpack('<4sHHHHIIH',eocd)
meta={'classic_eocd':dict(disk_no=disk_no,cd_disk=cd_disk,n_disk=n_disk,n_total=n_total,cd_size=cd_size32,cd_off=cd_off32,comment_len=comment_len)}
print('CLASSIC EOCD',meta['classic_eocd'])

# ZIP64 locator should precede classic EOCD when 32-bit fields overflow.
locpos=tail.rfind(b'PK\x06\x07',max(0,pos-256),pos)
if locpos>=0:
    sig,disk64,off64,total_disks=struct.unpack('<4sIQI',tail[locpos:locpos+20])
    zhdr=read_across(disk64,off64,56)
    if zhdr[:4]!=b'PK\x06\x06': raise RuntimeError('ZIP64 EOCD signature missing')
    sig,recsize,vermade,verneed,zdisk,zcd_disk,zn_disk,zn_total,zcd_size,zcd_off=struct.unpack('<4sQHHIIQQQQ',zhdr[:56])
    cd_disk=int(zcd_disk); cd_size=int(zcd_size); cd_off=int(zcd_off); n_total=int(zn_total)
    meta['zip64']=dict(locator_disk=disk64,locator_offset=off64,total_disks=total_disks,disk_no=zdisk,cd_disk=zcd_disk,n_disk=zn_disk,n_total=zn_total,cd_size=zcd_size,cd_off=zcd_off)
    print('ZIP64',meta['zip64'])
else:
    cd_disk=cd_disk; cd_size=cd_size32; cd_off=cd_off32

cd=read_across(cd_disk,cd_off,cd_size)
entries=[]; i=0
while i+46<=len(cd):
    if cd[i:i+4]!=b'PK\x01\x02':
        print('central directory stopped at',i,'of',len(cd),'sig',cd[i:i+4].hex()); break
    vals=struct.unpack('<4s6H3I5H2I',cd[i:i+46])
    _,vermade,verneed,flag,method,mtime,mdate,crc,csize32,usize32,nlen,xlen,clen,diskstart16,intattr,extattr,loff32=vals
    name=cd[i+46:i+46+nlen].decode('utf-8','replace')
    extra=cd[i+46+nlen:i+46+nlen+xlen]
    z=zip64_extra(extra,usize32==0xffffffff,csize32==0xffffffff,loff32==0xffffffff,diskstart16==0xffff)
    usize=int(z.get('usize',usize32)); csize=int(z.get('csize',csize32)); loff=int(z.get('loff',loff32)); diskstart=int(z.get('disk',diskstart16))
    entries.append(dict(name=name,compressed_size=csize,uncompressed_size=usize,method=method,flag=flag,disk_start=diskstart,local_offset=loff,crc32=crc))
    i+=46+nlen+xlen+clen
meta['parsed_entries']=len(entries)
out=Path('analysis/zatka_primary'); out.mkdir(parents=True,exist_ok=True)
(out/'split_zip_index.json').write_text(json.dumps({'metadata':meta,'entries':entries},indent=2))
for e in entries:
    n=e['name'].lower()
    if any(k in n for k in ['inactiv','opto','exp47','sessionlist','coordset','ctxoutline']):
        print(json.dumps(e))
print('entries',len(entries))
