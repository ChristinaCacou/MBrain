import json,shutil,tempfile,os
from pathlib import Path
import numpy as np,h5py
from remotezip import RemoteZip
URL='https://zenodo.org/records/17973874/files/2026_Ghaderi_NCOMMS_Zenodo_data_code.zip?download=1'
REGIONS=['A1','wS1','wS2','wM2','ALM'];EPOCHS=['audio','delay','whisker']
OUT=Path('analysis/ghaderi_raw_opto');OUT.mkdir(parents=True,exist_ok=True)
def rstr(ds):
    try:return np.asarray(ds.asstr()[()],dtype=str).ravel()
    except:
        a=ds[()];return np.asarray([z.decode(errors='replace') if isinstance(z,bytes) else str(z) for z in np.asarray(a).ravel()],dtype=str)
def main():
    rows=[];files=[]
    with RemoteZip(URL) as rz:
        members=sorted(z.filename for z in rz.infolist() if z.filename.startswith('data_optogenetics/') and z.filename.endswith('.nwb'))
        for ii,m in enumerate(members):
            local=f'/tmp/opto_{ii}.nwb'
            with rz.open(m) as src,open(local,'wb') as dst:shutil.copyfileobj(src,dst,8<<20)
            try:
                with h5py.File(local,'r') as f:
                    tr=f['/intervals/trials/trial_type'][()].astype(float).ravel();lick=f['/intervals/trials/lick_flag'][()].astype(bool).ravel();early=f['/intervals/trials/early_lick'][()].astype(bool).ravel();win=rstr(f['/intervals/trials/opto_window']);area=rstr(f['/intervals/trials/opto_area'])
                    n=min(len(tr),len(lick),len(early),len(win),len(area));tr=tr[:n];lick=lick[:n];early=early[:n];win=win[:n];area=area[:n];completed=~early;session=Path(m).stem
                    regs=sorted(set(a.strip() for a in area if a.strip() in REGIONS))
                    for reg in regs:
                        # Source sessions are area-specific; retain only trials carrying the session's opto-area label or nolight trials.
                        ar=(np.char.strip(area.astype(str))==reg)
                        base=completed&(tr==1)&(np.char.strip(win.astype(str))=='nolight')
                        if base.sum()<1:continue
                        b=float(lick[base].mean()*100.0)
                        for ep in EPOCHS:
                            cur=completed&(tr==1)&(np.char.strip(win.astype(str))==ep)
                            # opto_area can be nan/blank on no-light but should equal region on light trials.
                            if np.any(ar):cur &= ar
                            if cur.sum()<1:continue
                            p=float(lick[cur].mean()*100.0)
                            rows.append({'session':session,'member':m,'region':reg,'epoch':ep,'nolight_plick':b,'light_plick':p,'causal_impairment_pp':b-p,'n_nolight':int(base.sum()),'n_light':int(cur.sum())})
                    files.append({'member':m,'regions':regs,'n_trials':int(n)})
            finally:
                try:os.remove(local)
                except:pass
    out={'dataset':'Ghaderi raw optogenetic NWB trial tables','n_files':len(files),'rows':rows,'files':files,'method_note':'Direct NWB trial-table reconstruction of completed trial_type 1 lick probability; impairment is nolight minus matched light-window lick probability.'}
    (OUT/'opto_raw.json').write_text(json.dumps(out,indent=2));print(json.dumps({'n_files':len(files),'n_rows':len(rows)},indent=2))
if __name__=='__main__':main()
