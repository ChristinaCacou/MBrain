import json,itertools,csv
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
ROOT=Path('all-ghaderi-ephys');OP=Path('opto-artifact');OUT=Path('analysis/ghaderi_raw_final');OUT.mkdir(parents=True,exist_ok=True)
REG=['A1','wS1','wS2','wM2','ALM'];EPS=['audio','delay','whisker']
rows=[];files=[]
for p in ROOT.rglob('*.json'):
    try:d=json.loads(p.read_text())
    except:continue
    files.append({'session':d.get('session'),'member':d.get('member'),'n_rows':len(d.get('rows',[]))})
    for r in d.get('rows',[]):
        z=dict(r);z['session']=d.get('session');rows.append(z)
ephys=pd.DataFrame(rows)
ops=[]
for p in OP.rglob('opto_raw.json'):
    d=json.loads(p.read_text());ops.extend(d.get('rows',[]))
opto=pd.DataFrame(ops)
if ephys.empty:raise RuntimeError('no ephys rows')
if opto.empty:raise RuntimeError('no opto rows')
dec=ephys.groupby(['region','epoch'],as_index=False).agg(decoder_accuracy=('decoder_accuracy','mean'),n_decoder_sessions=('session','nunique'),median_n_units=('n_units','median'))
ca=opto.groupby(['region','epoch'],as_index=False).agg(causal_impairment_pp=('causal_impairment_pp','mean'),n_opto_sessions=('session','nunique'),n_light_trials=('n_light','sum'),n_nolight_trials=('n_nolight','sum'))
tab=dec.merge(ca,on=['region','epoch'])
tab=tab[tab.region.isin(REG)&tab.epoch.isin(EPS)].copy()
def stat(t,regs=REG,eps=EPS):
    xs=[];ys=[]
    for ep in eps:
        q=t[t.epoch==ep].set_index('region').reindex(regs).dropna(subset=['decoder_accuracy','causal_impairment_pp'])
        if len(q)<3:return np.nan
        rx=rankdata(q.decoder_accuracy.to_numpy(float));ry=rankdata(q.causal_impairment_pp.to_numpy(float));rx-=rx.mean();ry-=ry.mean();xs.append(rx);ys.append(ry)
    x=np.concatenate(xs);y=np.concatenate(ys);return float((x@y)/np.sqrt((x@x)*(y@y)))
def exact(t,regs=REG,eps=EPS):
    obs=stat(t,regs,eps);vals=[]
    for perm in itertools.permutations(regs):
        mp=dict(zip(regs,perm));ca=t[['region','epoch','causal_impairment_pp']].copy();ca['region']=ca.region.map(mp);u=t.drop(columns=['causal_impairment_pp']).merge(ca,on=['region','epoch'],how='left');vals.append(stat(u,regs,eps))
    vals=np.asarray(vals,float);return obs,float(np.mean(vals>=obs-1e-12)),vals.tolist()
obs,p,perms=exact(tab)
looR={}
for r in REG:
    regs=[x for x in REG if x!=r];q=tab[tab.region.isin(regs)];s,pp,_=exact(q,regs,EPS);looR[r]={'statistic':s,'exact_one_sided_p':pp}
looE={}
for e in EPS:
    eps=[x for x in EPS if x!=e];q=tab[tab.epoch.isin(eps)];s,pp,_=exact(q,REG,eps);looE[e]={'statistic':s,'exact_one_sided_p':pp}
result={'analysis_status':'POST_HOC_RAW_NWB_GHADERI_RECONSTRUCTION','dataset':'Ghaderi et al. Nature Communications 2026; Zenodo 17973874','n_ephys_files_parsed':len(files),'n_ephys_session_region_epoch_rows':len(ephys),'n_opto_trialtable_rows':len(opto),'n_region_epoch_cells':len(tab),'primary_raw_reconstruction_statistic':obs,'exact_common_region_label_one_sided_p':p,'leave_one_region_out':looR,'leave_one_epoch_out':looE,'method_note':'This is a post hoc raw-NWB reconstruction sensitivity, not a replacement for the prospectively specified processed-data analysis. It uses source task classes and fixed 50-neuron population size but fold-wise standardization and epoch-integrated raw spike counts.'}
tab.to_csv(OUT/'region_epoch_table.csv',index=False);ephys.to_csv(OUT/'ephys_session_metrics.csv',index=False);opto.to_csv(OUT/'opto_session_metrics.csv',index=False);(OUT/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));print(tab.to_string(index=False))
