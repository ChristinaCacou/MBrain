import json
from remotezip import RemoteZip
URL='https://zenodo.org/records/17973874/files/2026_Ghaderi_NCOMMS_Zenodo_data_code.zip?download=1'
with RemoteZip(URL) as rz:
    members=sorted(z.filename for z in rz.infolist() if z.filename.startswith('data_electrophysiology/') and z.filename.endswith('.nwb'))
print(json.dumps([{'member':m} for m in members],separators=(',',':')))
