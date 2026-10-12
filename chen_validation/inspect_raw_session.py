from pathlib import Path
import json, urllib.request, scipy.io, h5py
OUT=Path('analysis/chen_raw_schema'); OUT.mkdir(parents=True,exist_ok=True)
record=json.load(urllib.request.urlopen('https://zenodo.org/api/records/6713616'))
files={f['key']:f for f in record['files']}
want_base='BAYLORGC95_2019_12_19'
report={'record_file_count':len(files),'matches':[k for k in files if want_base in k]}
for prefix in ['data_structure_','meta_data_']:
    key=prefix+want_base+'.mat'
    if key not in files:
        report[key]='NOT_FOUND'; continue
    url=files[key]['links']['content']
    dest=Path('chen_validation')/key
    if not dest.exists(): urllib.request.urlretrieve(url,dest)
    report[key]={'size':dest.stat().st_size,'url':url}
    with open(dest,'rb') as fh: hdr=fh.read(128)
    report[key]['header']=hdr[:80].decode('latin1','replace')
    try:
        report[key]['whosmat']=[(a,b,list(c) if isinstance(c,tuple) else c) for a,b,c in scipy.io.whosmat(dest)]
        d=scipy.io.loadmat(dest,squeeze_me=True,struct_as_record=False)
        top={}
        for k,v in d.items():
            if k.startswith('__'): continue
            rec={'python_type':type(v).__name__}
            if hasattr(v,'_fieldnames'): rec['fields']=list(v._fieldnames)
            elif hasattr(v,'shape'): rec['shape']=list(v.shape); rec['dtype']=str(getattr(v,'dtype',''))
            top[k]=rec
        report[key]['top']=top
    except NotImplementedError:
        with h5py.File(dest,'r') as h:
            report[key]['h5_keys']=list(h.keys())
            report[key]['h5_shapes']={k:list(h[k].shape) if hasattr(h[k],'shape') else None for k in h.keys()}
(OUT/'schema.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
