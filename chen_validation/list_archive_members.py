import os, json
from remotezip import RemoteZip
idx=int(os.environ.get('ARCHIVE_INDEX','1'))
url=f'https://zenodo.org/records/6713616/files/datafiles_{idx}_of_8.zip?download=1'
standard=['BAYLORGC4','BAYLORGC12','BAYLORGC13','BAYLORGC15','BAYLORGC17','BAYLORGC18','BAYLORGC19','BAYLORGC21','BAYLORGC22','BAYLORGC25','BAYLORGC26','BAYLORGC50','BAYLORGC72','BAYLORGC75','BAYLORGC86','BAYLORGC87','BAYLORGC93','BAYLORGC95']
with RemoteZip(url) as rz:
    names=rz.namelist()
sel=[n for n in names if any(('data_structure_'+a+'_') in n or ('meta_data_'+a+'_') in n for a in standard)]
out={'archive':idx,'url':url,'n_all':len(names),'n_standard_members':len(sel),'members':sel}
os.makedirs('analysis/chen_raw_manifest',exist_ok=True)
open(f'analysis/chen_raw_manifest/archive_{idx}.json','w').write(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
