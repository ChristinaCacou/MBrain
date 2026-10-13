import json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr

SEED=20261012; NPERM=200000; NBOOT=50000
OUT=Path('analysis/chen_zero_refit_prediction'); OUT.mkdir(parents=True,exist_ok=True)
tr=pd.read_csv('analysis_inputs/chen_mouse_rows.csv')
te=pd.read_csv('analysis_inputs/chen2024_right_behavior_mouse_rows.csv')
model=LinearRegression().fit(tr[['coupling']],tr['right_behavior_recovery'])
pred=model.predict(te[['coupling']]); obs=te['right_behavior_recovery'].to_numpy(float)
r2=float(1-np.sum((obs-pred)**2)/np.sum((obs-obs.mean())**2)); mae=float(np.mean(np.abs(obs-pred))); rho=float(spearmanr(pred,obs).statistic)
rng=np.random.default_rng(SEED); ge=0; den=np.sum((obs-obs.mean())**2)
for _ in range(NPERM):
    y=rng.permutation(obs); rr=1-np.sum((y-pred)**2)/den; ge+=rr>=r2-1e-15
p=float((ge+1)/(NPERM+1))
rng=np.random.default_rng(SEED+1); r2b=[]; maeb=[]
for _ in range(NBOOT):
    ix=rng.integers(0,len(obs),len(obs)); y=obs[ix]; q=pred[ix]; d=np.sum((y-y.mean())**2)
    if d>0:r2b.append(1-np.sum((y-q)**2)/d)
    maeb.append(np.mean(np.abs(y-q)))
res={'status':'POST_HOC_ZERO_REFIT_EXTERNAL_PREDICTION','training':'Chen et al. 2021 mouse-level right-behavior model','test':'Chen et al. 2024 frozen-analysis right-behavior secondary mice','n_train':int(len(tr)),'n_test':int(len(te)),'slope':float(model.coef_[0]),'intercept':float(model.intercept_),'external_r2':r2,'mae':mae,'spearman_rho':rho,'permutation_p_r2':p,'bootstrap95_r2':[float(np.percentile(r2b,2.5)),float(np.percentile(r2b,97.5))],'bootstrap95_mae':[float(np.percentile(maeb,2.5)),float(np.percentile(maeb,97.5))],'interpretation':'Post hoc diagnostic only; does not convert the prespecified Chen 2024 secondary endpoint into a confirmatory replication.'}
(OUT/'result.json').write_text(json.dumps(res,indent=2,sort_keys=True)); te.assign(predicted=pred).to_csv(OUT/'predictions.csv',index=False); print(json.dumps(res,indent=2))
