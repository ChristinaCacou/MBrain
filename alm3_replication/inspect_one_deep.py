from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO
import requests, json, numpy as np
from pathlib import Path

out=Path('analysis/alm3_deep_schema'); out.mkdir(parents=True,exist_ok=True)
with DandiAPIClient() as client:
    ds=client.get_dandiset('000011','draft')
    assets=[a for a in ds.get_assets() if 'behavior+ecephys+ogen.nwb' in a.path][:3]
    rows=[]
    for ai,a in enumerate(assets):
        url=a.get_content_url(follow_redirects=1, strip_query=False)
        path=Path(f'alm3_replication/tmp_{ai}.nwb')
        with requests.get(url,stream=True,timeout=180) as r:
            r.raise_for_status()
            with open(path,'wb') as f:
                for chunk in r.iter_content(8*1024*1024):
                    if chunk: f.write(chunk)
        with NWBHDF5IO(str(path),'r',load_namespaces=True) as io:
            nwb=io.read(); df=nwb.trials.to_dataframe(); udf=nwb.units.to_dataframe()
            info={'asset_path':a.path,'subject':getattr(nwb.subject,'subject_id',None),'n_trials':len(df),'n_units':len(udf),
                  'trial_uniques':{},'acquisition':{},'ogen_sites':{},'unit_locations':[]}
            for c in df.columns:
                vals=df[c].astype(str).value_counts(dropna=False).head(30)
                info['trial_uniques'][c]={str(k):int(v) for k,v in vals.items()}
            for k,obj in nwb.acquisition.items():
                d={'class':obj.__class__.__name__,'description':getattr(obj,'description','')}
                if hasattr(obj,'time_series'):
                    d['time_series']={}
                    for n,ts in obj.time_series.items():
                        x={'class':ts.__class__.__name__,'description':getattr(ts,'description',''),'unit':getattr(ts,'unit',None)}
                        try:
                            data=np.asarray(ts.data[:]); x['data_shape']=list(data.shape); x['data_unique']=[str(z) for z in np.unique(data)[:30]]
                        except Exception as e: x['data_error']=repr(e)
                        try:
                            tt=np.asarray(ts.timestamps[:]); x['n_timestamps']=int(len(tt)); x['first_timestamps']=[float(z) for z in tt[:20]]
                        except Exception as e: x['timestamp_error']=repr(e)
                        d['time_series'][n]=x
                info['acquisition'][k]=d
            for k,site in nwb.ogen_sites.items():
                info['ogen_sites'][k]={'description':getattr(site,'description',''),'location':getattr(site,'location',''),'excitation_lambda':getattr(site,'excitation_lambda',None)}
            for eg in nwb.electrode_groups.values():
                info['unit_locations'].append({'name':eg.name,'location':eg.location,'description':eg.description})
            rows.append(info)
        path.unlink(missing_ok=True)
    out.joinpath('deep_schema.json').write_text(json.dumps(rows,indent=2,default=str)); print(json.dumps(rows,indent=2,default=str))
