import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, theilslopes, pearsonr

IN=Path('chen-final')
OUT=Path('analysis/chen_posthoc_robustness')
OUT.mkdir(parents=True,exist_ok=True)
mouse=pd.read_csv(IN/'mouse_rows.csv')
sessions=pd.DataFrame(json.loads((IN/'session_rows.json').read_text()))
primary=json.loads((IN/'result.json').read_text())

endpoints=['right_neuronal_recovery','right_behavior_recovery','left_neuronal_recovery','left_behavior_recovery']
pvals=np.array([primary[e]['one_sided_permutation_p'] for e in endpoints],float)
bonf=np.minimum(1,pvals*len(pvals))
order=np.argsort(pvals); bh=np.empty_like(pvals); prev=1.0
for pos in range(len(pvals)-1,-1,-1):
    idx=order[pos]; rank=pos+1; val=min(prev,pvals[idx]*len(pvals)/rank); bh[idx]=val; prev=val

cov=sessions.groupby('animal').agg(total_n_control=('n_control','sum'),total_n_rightstim=('n_rightstim','sum'),total_n_leftstim=('n_leftstim','sum')).reset_index()
df=mouse.merge(cov,on='animal',how='left')

def partial_rank(data,x,y,covs):
    z=data[[x,y]+covs].dropna().copy(); n=len(z)
    xr=rankdata(z[x].to_numpy(float)); yr=rankdata(z[y].to_numpy(float))
    C=np.column_stack([np.ones(n)]+[rankdata(z[c].to_numpy(float)) for c in covs])
    rx=xr-C@np.linalg.lstsq(C,xr,rcond=None)[0]
    ry=yr-C@np.linalg.lstsq(C,yr,rcond=None)[0]
    return float(pearsonr(rx,ry).statistic)

rows=[]
for j,e in enumerate(endpoints):
    x=df['coupling'].to_numpy(float); y=df[e].to_numpy(float)
    rho=float(spearmanr(x,y).statistic)
    loo=[]
    for i in range(len(df)):
        loo.append(float(spearmanr(np.delete(x,i),np.delete(y,i)).statistic))
    ts=theilslopes(y,x,alpha=.95)
    rows.append(dict(endpoint=e,n=len(df),rho=rho,raw_one_sided_permutation_p=float(pvals[j]),bonferroni_p=float(bonf[j]),bh_q=float(bh[j]),drop_one_mouse_rho_min=float(np.nanmin(loo)),drop_one_mouse_rho_max=float(np.nanmax(loo)),drop_one_mouse_rho_median=float(np.nanmedian(loo)),theil_sen_slope=float(ts.slope),theil_sen_95_low=float(ts.low_slope),theil_sen_95_high=float(ts.high_slope)))

right={
  'adjust_n_sessions': partial_rank(df,'coupling','right_behavior_recovery',['n_sessions']),
  'adjust_unit_counts': partial_rank(df,'coupling','right_behavior_recovery',['n_units_left','n_units_right']),
  'adjust_local_code_strength': partial_rank(df,'coupling','right_behavior_recovery',['transferable_local_cv_auc']),
  'adjust_trial_totals': partial_rank(df,'coupling','right_behavior_recovery',['total_n_control','total_n_rightstim']),
}
result={
 'analysis_status':'POST_HOC_CHEN_ROBUSTNESS_AFTER_PRIMARY_RESULT_KNOWN',
 'family':'four direct side-by-endpoint coupling tests',
 'multiplicity_note':'Bonferroni and Benjamini-Hochberg adjustments are applied to the four original permutation P values; these are post hoc robustness summaries, not a redefinition of the primary analysis.',
 'right_behavior_drop_one_mouse_note':'Each mouse is omitted in turn and the raw Spearman association is recomputed.',
 'partial_rank_note':'Descriptive rank-residualization sensitivities; no confirmatory P values assigned because they were specified after the primary result was known.',
 'endpoint_rows':rows,
 'right_behavior_partial_rank_sensitivities':right
}
pd.DataFrame(rows).to_csv(OUT/'endpoint_robustness.csv',index=False)
(OUT/'result.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
