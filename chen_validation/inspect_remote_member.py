import io, json, os, numpy as np
from remotezip import RemoteZip
from scipy.io import loadmat, whosmat
url='https://zenodo.org/records/6713616/files/datafiles_1_of_8.zip?download=1'
member='datafiles_1_of_8/data_structure_BAYLORGC12_2018_04_03.mat'
with RemoteZip(url) as rz:
    info=rz.getinfo(member)
    print('member',member,'compressed',info.compress_size,'uncompressed',info.file_size,flush=True)
    raw=rz.read(member)
print('downloaded bytes',len(raw),flush=True)
bio=io.BytesIO(raw)
print('whosmat',whosmat(bio),flush=True)
bio.seek(0); d=loadmat(bio,squeeze_me=True,struct_as_record=False)
obj=d.get('obj')
rep={'obj_type':type(obj).__name__,'fields':list(getattr(obj,'_fieldnames',[]))}
for field in rep['fields']:
    v=getattr(obj,field)
    x={'type':type(v).__name__}
    if hasattr(v,'shape'): x['shape']=list(v.shape)
    if field in ['trialTypeStr']:
        try:x['value']=[str(z) for z in np.atleast_1d(v).tolist()]
        except:pass
    if field in ['trialIds','trialStartTimes']:
        a=np.atleast_1d(v); x['preview']=a[:10].tolist(); x['n']=len(a)
    if field in ['trialPropertiesHash','eventSeriesHash']:
        x['subfields']=list(getattr(v,'_fieldnames',[]))
        if hasattr(v,'keyNames'):
            try:x['keyNames']=[str(z) for z in np.atleast_1d(v.keyNames).tolist()]
            except:pass
        if hasattr(v,'value'):
            vals=np.atleast_1d(v.value); x['n_values']=len(vals); x['value_types']=[type(z).__name__ for z in vals[:10]]
    rep[field]=x
# first unit fields
try:
    units=np.atleast_1d(obj.eventSeriesHash.value)
    u=units[0]
    rep['first_unit_fields']=list(u._fieldnames)
    for f in u._fieldnames:
        v=getattr(u,f); z={'type':type(v).__name__}
        if hasattr(v,'shape'): z['shape']=list(v.shape)
        if f in ['hemisphere','cellType','manual_quality_score','eventTrials']:
            try:
                a=np.atleast_1d(v); z['preview']=[str(q) for q in a[:10].tolist()]
            except: pass
        rep['unit_'+f]=z
except Exception as e: rep['unit_error']=repr(e)
os.makedirs('analysis/chen_raw_schema_remote',exist_ok=True)
open('analysis/chen_raw_schema_remote/schema.json','w').write(json.dumps(rep,indent=2,default=str))
print(json.dumps(rep,indent=2,default=str),flush=True)
