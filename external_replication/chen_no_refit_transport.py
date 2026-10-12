"""Post hoc zero-recalibration transport audit: Chen 2021 -> Chen 2024.

The 2021 raw-scale linear map is fit once, then applied to the 2024 right-side
behavioral endpoint without any slope/intercept refitting. This audit was conceived
after the Chen 2024 outcomes were known and is therefore descriptive/post hoc, not
confirmatory.
"""
from pathlib import Path
from itertools import permutations
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED=20261012
B=100000
OUT=Path('analysis/chen_no_refit_transport');OUT.mkdir(parents=True,exist_ok=True)
SRC=pd.read_csv('analysis_inputs/chen_mouse_rows.csv')
TGT=pd.read_csv('analysis_inputs/chen2024_mouse_endpoint_table.csv')

src=SRC[['coupling','right_behavior_recovery']].dropna().copy()
tgt=TGT[['subject','coupling_for_right_recovery','right_recovery']].dropna().copy()
x=src.coupling.to_numpy(float); y=src.right_behavior_recovery.to_numpy(float)
slope=float(np.cov(x,y,ddof=0)[0,1]/np.var(x)); intercept=float(y.mean()-slope*x.mean())
xt=tgt.coupling_for_right_recovery.to_numpy(float); yt=tgt.right_recovery.to_numpy(float); pred=intercept+slope*xt
n=len(yt)

def r2(a,b):
 d=np.sum((a-a.mean())**2)
 return float(1-np.sum((a-b)**2)/d) if d>0 else float('nan')

rho=float(spearmanr(pred,yt).statistic)
# Exact 8! one-sided label permutation for the fixed predictions.
ge=0;tot=0
for yp in permutations(yt):
 rr=float(spearmanr(pred,np.asarray(yp,float)).statistic)
 ge+=int(rr>=rho-1e-12);tot+=1
p_exact=ge/tot

cal_slope=float(np.cov(pred,yt,ddof=0)[0,1]/np.var(pred));cal_intercept=float(yt.mean()-cal_slope*pred.mean())
mae=float(np.mean(np.abs(yt-pred)));rmse=float(np.sqrt(np.mean((yt-pred)**2)))
source_mean=float(y.mean());mae_source_mean=float(np.mean(np.abs(yt-source_mean)))

rng=np.random.default_rng(SEED);idx=rng.integers(0,n,size=(B,n));yy=yt[idx];pp=pred[idx]
ym=yy.mean(axis=1);pm=pp.mean(axis=1);den=((yy-ym[:,None])**2).sum(axis=1);sse=((yy-pp)**2).sum(axis=1)
r2b=np.where(den>0,1-sse/den,np.nan);maeb=np.mean(np.abs(yy-pp),axis=1)
vp=np.mean((pp-pm[:,None])**2,axis=1);cv=np.mean((pp-pm[:,None])*(yy-ym[:,None]),axis=1)
cs=np.where(vp>0,cv/vp,np.nan);ci=ym-cs*pm

out={
 'status':'POST_HOC_ZERO_RECALIBRATION_TRANSPORT_AUDIT',
 'source':'Chen et al. 2021 raw mouse table','target':'Chen et al. 2024 DANDI frozen-analysis mouse table',
 'n_source_mice':int(len(src)),'n_target_mice':int(n),
 'source_fit':{'slope':slope,'intercept':intercept},
 'target_fixed_prediction':{
  'spearman_rho':rho,'exact_one_sided_permutation_p':float(p_exact),
  'external_r2_vs_target_mean':r2(yt,pred),'mae':mae,'rmse':rmse,
  'source_mean_baseline_mae':mae_source_mean,'mae_fractional_improvement_vs_source_mean':float(1-mae/mae_source_mean),
  'calibration_slope':cal_slope,'calibration_intercept':cal_intercept,
  'bootstrap95_r2':[float(x) for x in np.nanpercentile(r2b,[2.5,97.5])],
  'bootstrap95_mae':[float(x) for x in np.nanpercentile(maeb,[2.5,97.5])],
  'bootstrap95_calibration_slope':[float(x) for x in np.nanpercentile(cs,[2.5,97.5])],
  'bootstrap95_calibration_intercept':[float(x) for x in np.nanpercentile(ci,[2.5,97.5])],
 },
 'interpretation':'The raw-scale 2021 relation was transferred without refitting. Because this audit was specified after 2024 outcomes were known, it strengthens calibration/transport plausibility but cannot be treated as prospective confirmation.'
}
pd.DataFrame({'subject':tgt.subject,'coupling':xt,'observed_recovery':yt,'fixed_2021_prediction':pred}).to_csv(OUT/'predictions.csv',index=False)
(OUT/'result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
