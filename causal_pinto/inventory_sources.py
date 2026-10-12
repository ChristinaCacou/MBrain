import json, zipfile, hashlib
from pathlib import Path
import requests, h5py, numpy as np
from scipy.io import loadmat

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

def safe_sample(x, n=20):
    try:
        a=np.asarray(x)
        vals=a.reshape(-1)[:n]
        out=[]
        for v in vals:
            if isinstance(v,(bytes,np.bytes_)): out.append(v.decode(errors='replace'))
            elif isinstance(v,np.generic): out.append(v.item())
            else: out.append(v if isinstance(v,(str,int,float,bool)) else str(v))
        return out
    except Exception as e:
        return [f'<sample_error:{e}>']

def summarize_obj(obj, depth=0):
    if depth>8: return {'type':'depth_limit'}
    if isinstance(obj, dict):
        return {'type':'dict','keys':list(obj.keys()),'items':{str(k):summarize_obj(v,depth+1) for k,v in obj.items() if not str(k).startswith('__')}}
    if isinstance(obj, np.ndarray):
        rec={'type':'ndarray','shape':list(obj.shape),'dtype':str(obj.dtype)}
        if obj.dtype==object:
            flat=obj.reshape(-1)
            rec['items']=[summarize_obj(v,depth+1) for v in flat[:min(len(flat),25)]]
        elif obj.size <= 20000:
            rec['sample']=safe_sample(obj)
        return rec
    if isinstance(obj, (list,tuple)):
        return {'type':type(obj).__name__,'length':len(obj),'items':[summarize_obj(v,depth+1) for v in obj[:25]]}
    if isinstance(obj, np.generic): obj=obj.item()
    if isinstance(obj,(str,int,float,bool)) or obj is None: return {'type':type(obj).__name__,'value':obj}
    return {'type':type(obj).__name__,'repr':repr(obj)[:1000]}

def h5_inventory(path):
    out={}
    with h5py.File(path,'r') as f:
        def visitor(name,obj):
            if isinstance(obj,h5py.Dataset):
                rec={'type':'dataset','shape':list(obj.shape),'dtype':str(obj.dtype)}
                try:
                    if obj.size <= 5000: rec['sample']=safe_sample(obj[()])
                except Exception as e: rec['sample_error']=str(e)
                out['/'+name]=rec
            else: out['/'+name]={'type':'group'}
        f.visititems(visitor)
    return out

def mat_inventory(path):
    try:
        d=loadmat(path, simplify_cells=True)
        return summarize_obj(d)
    except NotImplementedError:
        return {'hdf5_fallback':h5_inventory(path)}

def main():
    inv={'source_urls':URLS,'downloads':{}}
    for key,url in URLS.items():
        zp=CACHE/f'{key}.zip'
        r=requests.get(url, timeout=120); r.raise_for_status(); zp.write_bytes(r.content)
        target=CACHE/key; target.mkdir(exist_ok=True)
        with zipfile.ZipFile(zp) as z:
            z.extractall(target); names=z.namelist()
        rec={'zip_sha256':sha256(zp),'zip_size':zp.stat().st_size,'members':names,'files':{}}
        for fp in sorted(target.rglob('*')):
            if not fp.is_file(): continue
            frec={'size':fp.stat().st_size,'sha256':sha256(fp)}
            suf=fp.suffix.lower()
            if suf in {'.h5','.hdf5'}: frec['hdf5']=h5_inventory(fp)
            elif suf=='.mat': frec['mat']=mat_inventory(fp)
            elif suf in {'.txt','.csv','.tsv','.json'} and fp.stat().st_size < 1_000_000:
                frec['text_head']=fp.read_text(errors='replace')[:12000]
            rec['files'][str(fp.relative_to(target))]=frec
        inv['downloads'][key]=rec
    (ROOT/'source_inventory.json').write_text(json.dumps(inv, indent=2, default=str))
    print(json.dumps(inv, indent=2, default=str)[:100000])

if __name__=='__main__': main()
