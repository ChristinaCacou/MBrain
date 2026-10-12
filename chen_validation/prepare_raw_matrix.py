import json, os
from remotezip import RemoteZip
standard=['BAYLORGC4','BAYLORGC12','BAYLORGC13','BAYLORGC15','BAYLORGC17','BAYLORGC18','BAYLORGC19','BAYLORGC21','BAYLORGC22','BAYLORGC25','BAYLORGC26','BAYLORGC50','BAYLORGC72','BAYLORGC75','BAYLORGC86','BAYLORGC87','BAYLORGC93','BAYLORGC95']
items=[]
for idx in [1,2,3,5]:
    url=f'https://zenodo.org/records/6713616/files/datafiles_{idx}_of_8.zip?download=1'
    with RemoteZip(url) as rz:
        for n in rz.namelist():
            base=n.split('/')[-1]
            if base.startswith('data_structure_') and any(base.startswith('data_structure_'+a+'_') for a in standard):
                items.append({'archive':str(idx),'member':n})
# three GC95 raw data files are exposed individually, while their metadata are in archive 6
for date in ['2019_12_13','2019_12_16','2019_12_19']:
    items.append({'archive':'direct','member':f'data_structure_BAYLORGC95_{date}.mat'})
# dedupe by basename in case a direct file also appears in a zip
seen=set(); out=[]
for x in items:
    b=x['member'].split('/')[-1]
    if b not in seen: seen.add(b); out.append(x)
print(json.dumps(out,separators=(',',':')))
