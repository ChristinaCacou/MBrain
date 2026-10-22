#!/usr/bin/env python3
"""Independent numerical cross-check of synthetic benchmark output.

Independent computations of eigenvalue participation ratios, paired t statistic,
Gaussian causal covariance identities, AUC rank statistic, and Monte Carlo
rejection fractions. Does not implement the locked scientific decision rule.
"""
from pathlib import Path
import json, numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
ROOT=Path(__file__).resolve().parent.parent
SYN=ROOT/'central_reproduction'
CAUSAL=ROOT/'causal_reproduction'
results={}
for scenario in ['common_driver','no_change']:
    d=SYN/scenario
    tr=pd.read_csv(d/'locked_outputs/trial_metrics.csv')
    meta=pd.read_csv(d/'SYNTHETIC_metadata.csv')
    sampled=meta.sample(n=min(48,len(meta)),random_state=477)
    errors=[]
    for r in sampled.itertuples(index=False):
        with np.load(d/'synthetic_trajectories'/r.trajectory_file) as f: x=f['trajectory']
        y=x-x[:20].mean(axis=0)
        norm=y/np.linalg.norm(y)
        cov=np.cov(norm.T)
        e=np.linalg.eigvalsh(cov);e=e[e>1e-12]
        dim=float(np.square(e.sum()) /np.square(e).sum())
        locked=float(tr.loc[(tr.cell_id==r.cell_id)&(tr.trial_id==r.trial_id),'effective_dimensionality'].iloc[0])
        errors.append(abs(dim-locked))
    t=pd.read_csv(d/'locked_outputs/primary_paired_differences.csv')
    diff=t['difference_hab_minus_base'].to_numpy(float)
    p=stats.ttest_1samp(diff,0).pvalue
    with open(d/'locked_outputs/confirmatory_results.json') as f: reported=json.load(f)
    results[scenario]={
      'num_trials':len(tr),'num_virtual_cells':tr.cell_id.nunique(),'n_days':tr.replicate_day.nunique(),
      'independent_eigenvalue_metric_max_abs_error':max(errors),
      'independent_paired_t_p':float(p),
      'reported_paired_t_p':reported['primary']['paired_t_p'],
      'reported_auc':reported['cross_modal']['auc'],
      'reported_mean_logscore_delta':reported['model_score']['mean_delta_log_score'],
      'reported_three_gate_decision':reported['decision']['local_scalar_gain_null_rejected']}
    # Recompute AUC with an independently assembled training/testing table,
    # using the exact preregistered estimator and a separately computed rank AUC.
    cell=pd.read_csv(d/'locked_outputs/cell_state_metrics.csv')
    ftrs=['effective_dimensionality','lz76_norm','compressed_complexity']
    trn=cell[(cell.treatment_code=='SHAM')&(cell.modality_code=='MECH')&(cell.state_code.isin(['BASE','HAB']))]
    tst=cell[(cell.treatment_code=='SHAM')&(cell.modality_code=='ELEC')&(cell.state_code.isin(['BASE','HAB']))]
    sc=StandardScaler().fit(trn[ftrs].values)
    clf=LogisticRegression(C=1.0,solver='liblinear',random_state=260416,max_iter=1000)
    clf.fit(sc.transform(trn[ftrs].values),(trn.state_code=='HAB').astype(int).values)
    y=(tst.state_code=='HAB').astype(int).values
    scores=clf.predict_proba(sc.transform(tst[ftrs].values))[:,1]
    auc_by_pairwise=float((scores[y==1,None]>scores[y==0]).mean()+.5*(scores[y==1,None]==scores[y==0]).mean())
    results[scenario]['independent_auc_by_rank']=auc_by_pairwise
    results[scenario]['max_auc_difference']=abs(auc_by_pairwise-reported['cross_modal']['auc'])
# Structural equation independent algebra audit
rho=.60; gamma=.25
c=(rho-gamma)/np.sqrt(rho)
v=1-gamma**2-c**2-2*gamma*c*np.sqrt(rho)
var_E=gamma**2+c**2+2*gamma*c*np.sqrt(rho)+v
cov_ME=gamma+c*np.sqrt(rho)
results['causal_algebra']={'common_cause_covariance':rho,'coupled_covariance':cov_ME,'coupled_E_variance':var_E,'residual_variance':v,'intervention_delta_E_common_cause':0.,'intervention_delta_E_coupled':gamma}
summary=pd.read_csv(CAUSAL/'ExtendedData_Table9b_intervention_assurance.csv')
raw=pd.read_csv(CAUSAL/'ExtendedData_Table9c_intervention_trials.csv')
results['campaigns']={}
for row in summary.itertuples():
    z=raw[raw.scenario==row.scenario]
    rate=float((z.pvalue<.05).mean())
    results['campaigns'][row.scenario]={'n_campaigns':len(z),'calculated_rejection_rate':rate,'reported_rate':float(row.rejection_rate),'difference':abs(rate-float(row.rejection_rate))}

assert max(results[x]['independent_eigenvalue_metric_max_abs_error'] for x in ['common_driver','no_change'])<1e-9
assert max(results[x]['max_auc_difference'] for x in ['common_driver','no_change'])<1e-12
assert abs(cov_ME-rho)<1e-12 and abs(var_E-1)<1e-12 and v>0
assert all(v['difference']<1e-15 for v in results['campaigns'].values())
(ROOT/'AUDIT_NUMERICAL_RESULTS.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
