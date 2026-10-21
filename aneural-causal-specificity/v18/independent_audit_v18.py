#!/usr/bin/env python3
"""Independent post-hoc diagnostic of the stored synthetic experiment.
Does not call or edit the frozen analysis pipeline. Checks its input-independent
properties and the consequence of treating recording day as unit of replication.
"""
from pathlib import Path
import argparse, hashlib, json, itertools
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

P='common_driver'; C='no_change'
FEATS=['effective_dimensionality','lz76_norm','compressed_complexity']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def auc_dayout(df):
    days=sorted(df.replicate_day.unique());out=[]
    for day in days:
        train=df[(df.replicate_day!=day)&(df.modality_code=='MECH')]
        test=df[(df.replicate_day==day)&(df.modality_code=='ELEC')]
        scaler=StandardScaler().fit(train[FEATS].to_numpy())
        model=LogisticRegression(C=1,solver='liblinear',max_iter=1000,random_state=260416)
        model.fit(scaler.transform(train[FEATS]),(train.state_code=='HAB').astype(int))
        y=(test.state_code=='HAB').astype(int)
        pred=model.predict_proba(scaler.transform(test[FEATS]))[:,1]
        out.append({'day':day,'n_electrical_rows':len(test),'auc':roc_auc_score(y,pred)})
    return out

def score_dayout(df):
    # Independent reconstruction of mechanical delta models using electrical baseline and heldout HAB.
    from scipy.stats import multivariate_normal
    rows=[];days=sorted(df.replicate_day.unique())
    for day in days:
        tr=df[df.replicate_day!=day];te=df[df.replicate_day==day]
        m0=tr[(tr.modality_code=='MECH')&(tr.state_code=='BASE')].set_index('cell_id')[FEATS]
        m1=tr[(tr.modality_code=='MECH')&(tr.state_code=='HAB')].set_index('cell_id')[FEATS]
        ids=m0.index.intersection(m1.index)
        delta=(m1.loc[ids]-m0.loc[ids]).mean(axis=0).to_numpy()
        resid=(m1.loc[ids]-m0.loc[ids]).to_numpy()-delta
        cov=np.cov(resid,rowvar=False)+np.eye(len(FEATS))*1e-6
        e0=te[(te.modality_code=='ELEC')&(te.state_code=='BASE')].set_index('cell_id')[FEATS]
        e1=te[(te.modality_code=='ELEC')&(te.state_code=='HAB')].set_index('cell_id')[FEATS]
        ds=[]
        for cid in e0.index.intersection(e1.index):
            b=e0.loc[cid].to_numpy(float); o=e1.loc[cid].to_numpy(float)
            ds.append(float(multivariate_normal.logpdf(o,mean=b+delta,cov=cov)-multivariate_normal.logpdf(o,mean=b,cov=cov)))
        rows.append({'day':day,'n_cells':len(ds),'mean_delta_log_score':float(np.mean(ds))})
    return rows

def exact_signflip(values, alternative='greater'):
    """Unstudentized symmetric sign-flip reference, conditional on four observed day values."""
    v=np.asarray(values,float);obs=np.mean(v)
    null=np.array([np.mean(v*np.asarray(z)) for z in itertools.product([-1,1],repeat=len(v))])
    if alternative=='greater': return float(np.mean(null>=obs-1e-12)),len(null)
    return float(np.mean(np.abs(null)>=abs(obs)-1e-12)),len(null)

def run(root, case):
    p=Path(root)/case/'locked_outputs'/'cell_state_metrics.csv'
    d=pd.read_csv(p)
    electrical=d[(d.modality_code=='ELEC')&(d.treatment_code=='SHAM')]
    pv=electrical.pivot(index=['cell_id','replicate_day'],columns='state_code',values='effective_dimensionality').dropna()
    pv['delta']=pv.HAB-pv.BASE
    day=pv.groupby('replicate_day').delta.agg(['mean','count']).reset_index()
    primary= float(pv.delta.mean()); day_vals=day['mean'].to_list()
    day_p_two, n=exact_signflip(day_vals,'two')
    day_p_directional,_=exact_signflip([-v for v in day_vals],'greater')
    decode=auc_dayout(d)
    predicted=score_dayout(d)
    pred_exact,n2=exact_signflip([x['mean_delta_log_score'] for x in predicted],'greater')
    return {'case':case,'n_cells':len(pv),'day_geometry':day.to_dict(orient='records'),
            'primary_mean_geometry_delta':primary,'day_signflip_two_sided':day_p_two,
            'day_signflip_directional_decrease':day_p_directional,
            'day_signflip_support_n':n,'crossmodal_dayout':decode,
            'crossmodal_dayout_pooled_auc_mean':float(np.mean([x['auc'] for x in decode])),
            'predictive_dayout':predicted,'predictive_day_signflip_one_sided':pred_exact,
            'predictive_day_signflip_support_n':n2}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run-dir',required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args()
    out={'scope':'synthetic, post-hoc diagnostic, not preregistered; no cells acquired',
         'cases':[run(a.run_dir,P),run(a.run_dir,C)]}
    Path(a.out).write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))