from pathlib import Path
import urllib.request, h5py, json
OUT=Path('analysis/alm3_schema'); OUT.mkdir(parents=True,exist_ok=True)
URL='https://dandiarchive.s3.amazonaws.com/blobs/a1d/c8d/a1dc8d9a-bbb6-4260-b56e-78844dbd957e'
F=Path('alm3_replication/example.nwb')
if not F.exists(): urllib.request.urlretrieve(URL,F)
rows=[]
with h5py.File(F,'r') as h:
    def visit(name,obj):
        rec={'path':name,'type':'group' if isinstance(obj,h5py.Group) else 'dataset'}
        if isinstance(obj,h5py.Dataset):
            rec['shape']=list(obj.shape); rec['dtype']=str(obj.dtype)
        # include only structural text attributes, never numeric dataset values
        attrs={}
        for k,v in obj.attrs.items():
            if isinstance(v,(bytes,str)):
                try: attrs[k]=v.decode() if isinstance(v,bytes) else v
                except: attrs[k]=str(v)
        if attrs: rec['text_attrs']=attrs
        rows.append(rec)
    h.visititems(visit)
    # print dynamic-table column names structurally, but not any cell values
    table_info={}
    for p in ['/intervals/trials','/units','/general/optogenetics']:
        if p in h:
            g=h[p]
            table_info[p]={'keys':list(g.keys()),'attrs':{k:(v.decode() if isinstance(v,bytes) else str(v)) for k,v in g.attrs.items() if isinstance(v,(bytes,str))}}
(OUT/'schema.json').write_text(json.dumps({'source_url':URL,'objects':rows,'tables':table_info},indent=2))
print(json.dumps({'tables':table_info,'n_objects':len(rows)},indent=2))
