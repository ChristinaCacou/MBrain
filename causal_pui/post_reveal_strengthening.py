from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

SEED=20261011
ROOT=Path('analysis/pui_blind_fixedC')
OUT=ROOT/'post_reveal_strengthening'; OUT.mkdir(parents=True,exist_ok=True)
sess=pd.read_csv(ROOT/'blind_features_session.csv')
rev=pd.read_csv(ROOT/'reveal_vectorized/revealed_region_table.csv')
regions=rev.region.tolist()
outcome=rev.set_index('region').loc[regions,'silencing_effect'].to_numpy(float)

# Equal-mouse weighting: first average sessions within mouse and region, then mice within region.
mm=(sess.groupby(['region','mouse'],as_index=False)
    .agg(pui=('pui_nats_per_trial','mean'), marginal=('marginal_gain_nats_per_trial','mean'),
         auc=('target_only_auc','mean'), bacc=('target_only_bacc','mean')))
eq=mm.groupby('region',as_index=False).mean(numeric_only=True).set_index('region').loc[regions]

def perm_p(x,y,n=500000,seed=SEED):
    rng=np.random.default_rng(seed); obs=float(spearmanr(x,y).statistic)
    ry=rankdata(y); rx=rankdata(x); rx=rx-rx.mean(); denx=np.sqrt(np.sum(rx**2)); vals=np.empty(n)
    for i in range(n):
        yp=ry[rng.permutation(len(ry))]; yp=yp-yp.mean(); vals[i]=np.dot(rx,yp)/(denx*np.sqrt(np.sum(yp**2)))
    return obs,float((np.sum(vals>=obs-1e-12)+1)/(n+1))

def loo(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); pred=np.full(len(y),np.nan)
    for i in range(len(y)):
        tr=np.arange(len(y))!=i
        X=np.column_stack([np.ones(tr.sum()),x[tr]])
        b=np.linalg.lstsq(X,y[tr],rcond=None)[0]
        pred[i]=b[0]+b[1]*x[i]
    rho=float(spearmanr(pred,y).statistic); mse=float(np.mean((pred-y)**2)); var=float(np.mean((y-y.mean())**2)); r2=float(1-mse/var)
    return rho,r2,pred

metrics={}
for col in ['pui','marginal','auc','bacc']:
    x=eq[col].to_numpy(float); rho,p=perm_p(x,outcome); lr,lr2,pred=loo(x,outcome)
    metrics[col]={'rho':rho,'one_sided_perm_p':p,'loo_rho':lr,'loo_r2':lr2,'values':x.tolist(),'loo_predictions':pred.tolist()}

# Reproduce split-half stability using frozen session features only. Regions with >=4 units.
def split_stability(df, unit, value='pui_nats_per_trial', n=3000):
    rng=np.random.default_rng(SEED)
    if unit=='mouse':
        base=df.groupby(['region','mouse'],as_index=False)[value].mean()
        idcol='mouse'
    else:
        base=df[['region','session',value]].copy(); idcol='session'
    counts=base.groupby('region')[idcol].nunique(); keep=counts[counts>=4].index.tolist()
    vals=[]
    for _ in range(n):
        a=[]; b=[]
        for r in keep:
            q=base[base.region==r]
            ids=q[idcol].unique().copy(); rng.shuffle(ids)
            h=len(ids)//2
            if h<1: continue
            ia=ids[:h]; ib=ids[h:2*h]
            a.append(q[q[idcol].isin(ia)][value].mean()); b.append(q[q[idcol].isin(ib)][value].mean())
        vals.append(float(spearmanr(a,b).statistic))
    arr=np.asarray(vals,float)
    return {'n_regions':len(keep),'median_rho':float(np.nanmedian(arr)),'split_2p5':float(np.nanpercentile(arr,2.5)),'split_97p5':float(np.nanpercentile(arr,97.5))}

rel_session=split_stability(sess,'session')
rel_mouse=split_stability(sess,'mouse')
# Spearman-Brown full-length reliability and an upper-bound sensitivity if the outcome were perfectly reliable.
def atten(rel):
    r=rel['median_rho']; sb=2*r/(1+r) if r>-1 else np.nan
    return {'half_split_median_rho':r,'spearman_brown_full_reliability':float(sb),'max_observable_abs_correlation_if_outcome_reliability_1':float(np.sqrt(max(sb,0)))}
atten_session=atten(rel_session); atten_mouse=atten(rel_mouse)

result={'equal_mouse_weighting':metrics,'session_split':rel_session,'mouse_split':rel_mouse,
        'attenuation_sensitivity_session':atten_session,'attenuation_sensitivity_mouse':atten_mouse,
        'scope_note':'All analyses are post-reveal diagnostics. Spearman-Brown/attenuation values are heuristic sensitivities because split-half rank stability is not a classical measurement reliability coefficient and silencing-outcome reliability is unknown.'}
(OUT/'result.json').write_text(json.dumps(result,indent=2))
pd.DataFrame({'region':regions,'silencing_effect':outcome,'pui_equal_mouse':eq.pui.values,'marginal_equal_mouse':eq.marginal.values,'auc_equal_mouse':eq.auc.values,'bacc_equal_mouse':eq.bacc.values}).to_csv(OUT/'equal_mouse_region_table.csv',index=False)
print(json.dumps(result,indent=2))
