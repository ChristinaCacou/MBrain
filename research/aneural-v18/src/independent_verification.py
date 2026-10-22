#!/usr/bin/env python3
"""Separate numerical verification of the central v18 synthetic result.
Does NOT import the preregistered analysis functions; reproduces response dimensionality
from saved trajectories and applies a four-day exact sign-flip sensitivity test.
"""
from pathlib import Path
import hashlib, json, itertools
import numpy as np
import pandas as pd
from scipy import stats
ROOT=Path(__file__).parent
BASE=ROOT/'frozen_audit'

def dim(x):
    y=x-x[:20].mean(0,keepdims=True)
    y=y/np.sqrt(np.sum(y*y))
    c=np.cov(y,rowvar=False)
    eig=np.linalg.eigvalsh(c)
    eig=eig[eig>1e-12]
    return float(eig.sum()**2/np.dot(eig,eig))

out={}
for case in ['common_driver','no_change']:
    folder=BASE/case
    meta=pd.read_csv(folder/'SYNTHETIC_metadata.csv')
    means={};records=[]
    for _,r in meta.iterrows():
        if r.modality_code!='ELEC':continue
        path=folder/'synthetic_trajectories'/r.trajectory_file
        with np.load(path) as f: x=f['trajectory']
        records.append({'day':r.replicate_day,'cell':r.cell_id,'state':r.state_code,'dim':dim(x)})
    t=pd.DataFrame(records)
    piv=t.pivot_table(index=['day','cell'],columns='state',values='dim').reset_index()
    d=piv.HAB-piv.BASE
    day=piv.assign(diff=d).groupby('day')['diff'].mean().to_dict()
    dayvals=np.array(list(day.values()))
    sign_means=np.array([np.mean(dayvals*np.array(s)) for s in itertools.product([-1.,1.],repeat=len(dayvals))])
    exact_p=float(np.mean(np.abs(sign_means)>=abs(dayvals.mean())-1e-12))
    tt=stats.ttest_rel(piv.HAB,piv.BASE)
    reference=pd.read_csv(folder/'locked_outputs'/'primary_paired_differences.csv')
    deviation=float(np.max(np.abs(piv.sort_values(['day','cell'])['dim'].to_numpy()[:0]))) if False else float(np.max(np.abs(piv.BASE.to_numpy()-reference.BASE.to_numpy())))
    output={
        'n_virtual_cells':len(piv),'baseline_mean_independent':float(piv.BASE.mean()),
        'habituated_mean_independent':float(piv.HAB.mean()),
        'paired_mean_difference':float(d.mean()),'cell_level_t_p':float(tt.pvalue),
        'recomputed_vs_locked_BASE_max_absolute_difference':deviation,
        'day_level_mean_differences':day,'exact_two_sided_day_signflip_p':exact_p,
        'n_exact_day_sign_assignments':2**len(dayvals)
    }
    assert deviation < 1e-8, deviation
    out[case]=output
rho=.6;gamma=.25;c=(rho-gamma)/rho;resid=1-(gamma**2+c*c*rho+2*gamma*c*rho)
cov=gamma+c*rho
variance=gamma**2+c*c*rho+2*gamma*c*rho+resid
out['analytic_structural_model']={'rho_target':rho,'gamma':gamma,'common_cause_covariance_ME':rho,'directed_model_covariance_ME':cov,'directed_model_var_E':variance,'residual_variance':resid,'directed_intervention_derivative':gamma,'common_cause_intervention_derivative':0.0}
ROOT.joinpath('independent_verification_results.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))