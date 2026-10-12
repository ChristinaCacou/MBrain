import json, zipfile, hashlib
from pathlib import Path
import requests, h5py, numpy as np

URLS = {
    'fig2': 'https://cdn.elifesciences.org/articles/70263/elife-70263-fig2-data1-v1.zip',
    'fig4': 'https://cdn.elifesciences.org/articles/70263/elife-70263-fig4-data1-v1.zip',
}
ROOT = Path('analysis/pinto_external')
CACHE = ROOT/'source'
CACHE.mkdir(parents=True, exist_ok=True)

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def safe_sample(x, n=12):
    try:
        a=np.asarray(x)
        if a.dtype.kind in 'SOU':
            vals=a.reshape(-1)[:n]
            out=[]
            for v in vals:
                if isinstance(v,(bytes,np.bytes_)): out.append(v.decode(errors='replace'))
                elif isinstance(v,np.generic): out.append(v.item())
                else: out.append(str(v))
            return out
        vals=a.reshape(-1)[:n]
        return [v.item() if isinstance(v,np.generic) else v for v in vals]
    except Exception as e:
        return [f'<sample_error:{e}>']

def h5_inventory(path):
    out={}
    with h5py.File(path,'r') as f:
        def visitor(name,obj):
            if isinstance(obj,h5py.Dataset):
                rec={'type':'dataset','shape':list(obj.shape),'dtype':str(obj.dtype)}
                try:
                    if obj.size <= 5000:
                        rec['sample']=safe_sample(obj[()])
                except Exception as e:
                    rec['sample_error']=str(e)
                out['/'+name]=rec
            else:
                out['/'+name]={'type':'group'}
        f.visititems(visitor)
    return out

def main():
    inv={'source_urls':URLS,'downloads':{}}
    for key,url in URLS.items():
        zp=CACHE/f'{key}.zip'
        r=requests.get(url, timeout=120)
        r.raise_for_status(); zp.write_bytes(r.content)
        target=CACHE/key; target.mkdir(exist_ok=True)
        with zipfile.ZipFile(zp) as z:
            z.extractall(target)
            names=z.namelist()
        rec={'zip_sha256':sha256(zp),'zip_size':zp.stat().st_size,'members':names,'files':{}}
        for fp in sorted(target.rglob('*')):
            if not fp.is_file(): continue
            frec={'size':fp.stat().st_size,'sha256':sha256(fp)}
            suf=fp.suffix.lower()
            if suf in {'.h5','.hdf5'}:
                frec['hdf5']=h5_inventory(fp)
            elif suf in {'.txt','.csv','.tsv','.json'} and fp.stat().st_size < 1_000_000:
                frec['text_head']=fp.read_text(errors='replace')[:12000]
            rec['files'][str(fp.relative_to(target))]=frec
        inv['downloads'][key]=rec
    (ROOT/'source_inventory.json').write_text(json.dumps(inv, indent=2, default=str))
    print(json.dumps(inv, indent=2, default=str)[:50000])

if __name__=='__main__': main()
