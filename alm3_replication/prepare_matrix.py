import json
from dandi.dandiapi import DandiAPIClient
VERSION='0.220126.1907'
with DandiAPIClient() as client:
    ds=client.get_dandiset('000011',VERSION)
    paths=sorted(a.path for a in ds.get_assets() if 'behavior+ecephys+ogen.nwb' in a.path)
print(json.dumps({'include':[{'asset':p,'id':f'{i:03d}'} for i,p in enumerate(paths)]},separators=(',',':')))
