from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO
import requests, os, json
from pathlib import Path

out=Path('analysis/alm3_schema'); out.mkdir(parents=True,exist_ok=True)
with DandiAPIClient() as client:
    ds=client.get_dandiset('000011','draft')
    assets=[a for a in ds.get_assets() if 'behavior+ecephys+ogen.nwb' in a.path]
    if not assets:
        raise RuntimeError('No behavior+ecephys+ogen NWB assets found')
    a=assets[0]
    url=a.get_content_url(follow_redirects=1, strip_query=False)
    path=Path('alm3_replication/tmp.nwb')
    with requests.get(url,stream=True,timeout=120) as r:
        r.raise_for_status()
        with open(path,'wb') as f:
            for chunk in r.iter_content(8*1024*1024):
                if chunk: f.write(chunk)

with NWBHDF5IO(str(path),'r',load_namespaces=True) as io:
    nwb=io.read()
    info={
      'asset_path': a.path,
      'identifier': nwb.identifier,
      'session_description': nwb.session_description,
      'subject_id': getattr(getattr(nwb,'subject',None),'subject_id',None),
      'trials_columns': list(nwb.trials.colnames) if nwb.trials is not None else [],
      'units_columns': list(nwb.units.colnames) if nwb.units is not None else [],
      'acquisition_keys': list(nwb.acquisition.keys()),
      'stimulus_keys': list(nwb.stimulus.keys()),
      'processing_keys': list(nwb.processing.keys()),
      'interval_keys': list(nwb.intervals.keys()),
      'ogen_sites': list(nwb.ogen_sites.keys()) if getattr(nwb,'ogen_sites',None) is not None else [],
      'n_trials': len(nwb.trials) if nwb.trials is not None else 0,
      'n_units': len(nwb.units) if nwb.units is not None else 0,
    }
    if nwb.trials is not None:
      df=nwb.trials.to_dataframe()
      info['trial_dtypes']={c:str(df[c].dtype) for c in df.columns}
      info['trial_samples']={c:[str(x) for x in df[c].dropna().head(8).tolist()] for c in df.columns}
    if nwb.units is not None:
      udf=nwb.units.to_dataframe()
      info['unit_dtypes']={c:str(udf[c].dtype) for c in udf.columns}
      info['unit_samples']={c:[str(x) for x in udf[c].dropna().head(5).tolist()] for c in udf.columns if c!='spike_times'}
    out.joinpath('schema.json').write_text(json.dumps(info,indent=2,default=str))
    print(json.dumps(info,indent=2,default=str))
path.unlink(missing_ok=True)
