#!/usr/bin/env python3
"""Independent numerical audit of the frozen common-driver positive case.
Not preregistered confirmatory inference. Operates on generated *synthetic* data.
Must not be used as a biological validation result.
"""
from pathlib import Path
from itertools import product
import json, hashlib
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

BASE=Path(__file__).resolve().parent/'replication'
source=BASE/'common_driver/locked_outputs/cell_state_metrics.csv'
c=pd.read_csv(source)
feats=['effective_dimensionality','lz76_norm','compressed_complexity']
q=c.query("treatment_code=='SHAM' and modality_code=='ELEC'").pivot(index=['cell_id','replicate_day'],columns='state_code',values='effective_dimensionality').dropna()
d=(q['HAB']-q['BASE']).to_numpy()
dz=float(d.mean()/d.std(ddof=1))
tp=float(stats.ttest_1samp(d,0).pvalue)
day_d=q.assign(diff=q['HAB']-q['BASE']).reset_index().groupby('replicate_day')['diff'].mean().to_numpy()
null=[abs(np.mean(day_d*np.array(s))) for s in product([-1,1],repeat=len(day_d))]
pday=sum(x>=abs(day_d.mean())-1e-12 for x in null)/len(null)

mech=c.query("treatment_code=='SHAM' and modality_code=='MECH' and state_code in ['BASE','HAB']")
elec=c.query("treatment_code=='SHAM' and modality_code=='ELEC' and state_code in ['BASE','HAB']")
s=StandardScaler().fit(mech[feats]); clf=LogisticRegression(C=1.0,solver='liblinear',random_state=260416,max_iter=1000).fit(s.transform(mech[feats]),(mech.state_code=='HAB').astype(int))
auc=float(roc_auc_score((elec.state_code=='HAB').astype(int),clf.predict_proba(s.transform(elec[feats]))[:,1]))

# Held-out predictive score, independently implemented from the model description.
log_scores=[]
for day in sorted(c.replicate_day.unique()):
 train=c[c.replicate_day!=day]; test=c[c.replicate_day==day]
 m=train[train.modality_code=='MECH'].pivot(index='cell_id',columns='state_code',values=feats)
 shift=np.column_stack([(m[(feat,'HAB')]-m[(feat,'BASE')]).to_numpy() for feat in feats]); shift=shift[np.isfinite(shift).all(1)]
 mu=shift.mean(0)
 cov=np.cov(shift,rowvar=False)+np.eye(3)*1e-6
 inv=np.linalg.pinv(cov)
 e=test[test.modality_code=='ELEC'].pivot(index='cell_id',columns='state_code',values=feats)
 ds=np.column_stack([(e[(feat,'HAB')]-e[(feat,'BASE')]).to_numpy() for feat in feats]); ds=ds[np.isfinite(ds).all(1)]
 for v in ds:
  log_scores.append({'day':day,'score':float(-.5*(v-mu)@inv@(v-mu)+.5*v@inv@v)})
pred=float(np.mean([x['score'] for x in log_scores]))
day_scores=pd.DataFrame(log_scores).groupby('day')['score'].mean().to_numpy(); pday_predict=sum(abs(np.mean(day_scores*np.array(s)))>=abs(day_scores.mean())-1e-12 for s in product([-1,1],repeat=len(day_scores)))/2**len(day_scores)

# Independently recompute amplitude-normalized electrical *dimensionality* for both cases.
import csv
r={}
for group in ['common_driver','no_change']:
 path=BASE/group
 meta=pd.read_csv(path/'SYNTHETIC_metadata.csv').query("modality_code=='ELEC'")
 rows=[]
 for rec in meta.itertuples(index=False):
  with np.load(path/'synthetic_trajectories'/rec.trajectory_file) as z: x=z['trajectory']
  y=x-x[:20].mean(0); y=y/np.sqrt((y*y).sum()); eig=np.linalg.eigvalsh(np.cov(y,rowvar=False));eig=eig[eig>1e-12]
  dim=float(eig.sum()**2/(eig*eig).sum())
  rows.append((rec.cell_id,rec.state_code,dim))
 pp=pd.DataFrame(rows,columns=['id','state','dimension']).pivot(index='id',columns='state',values='dimension');delta=pp.HAB-pp.BASE
 r[group]={'base':float(pp.BASE.mean()),'hab':float(pp.HAB.mean()),'paired_delta':float(delta.mean()),'paired_t_p':float(stats.ttest_rel(pp.HAB,pp.BASE).pvalue),'n':len(pp)}

result={'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'status':'synthetic numerical audit','positive_case':{'n_cells':len(q),'baseline_dim':float(q.BASE.mean()),'post_dim':float(q.HAB.mean()),'paired_dz':dz,'paired_t_p':tp,'mechanical_to_electrical_auc':auc,'heldout_mean_delta_log_score':pred,'day_means_electrical_dim':day_d.tolist(),'day_sign_flip_p_two_sided':pday,'day_means_predictive':day_scores.tolist(),'day_predictive_sign_flip_p_two_sided':pday_predict},'geometry_independent_trajectory_recomputation':r,'note':'P values of frozen pipeline are cell-level; day-level exact test differs. No consciousness or biological claim is made.'}
(Path(__file__).resolve().parent/'independent_audit_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))