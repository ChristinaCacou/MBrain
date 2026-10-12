import json,csv
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,rankdata,theilslopes

SEED=20261012;NPERM=500000;NBOOT=50000
ROOT=Path('all-robust-results');OUT=Path('analysis/chen_robust_final');OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for p in ROOT.rglob('*.json'):
    try:r=json.loads(p.read_text())
    except:continue
    if 'error' not in r:rows.append(r)
def finite(v):
    try:return v is not None and np.isfinite(float(v))
    except:return False
animals=sorted(set(r['animal'] for r in rows));mice=[]
keys=['control_bilateral_coupling','control_bilateral_within_condition_coupling','control_bilateral_coupling_fixed10','control_bilateral_coupling_fixed20','transferable_local_cv_auc','left_neuronal_recovery','right_neuronal_recovery','left_behavior_recovery','right_behavior_recovery','left_behavior_delta','right_behavior_delta','left_neuronal_delta','right_neuronal_delta','control_behavior_accuracy','left_stim_behavior_accuracy','right_stim_behavior_accuracy','n_units_left','n_units_right','n_control','n_leftstim','n_rightstim']
for a in animals:
    q=[r for r in rows if r['animal']==a];d={'animal':a,'n_sessions':len(q)}
    for k in keys:
        z=[float(r[k]) for r in q if finite(r.get(k))];d[k]=float(np.mean(z)) if z else None
        if k in ['control_bilateral_coupling_fixed10','control_bilateral_coupling_fixed20']:d[k+'_n_sessions']=len(z)
    if finite(d.get('n_units_left')) and finite(d.get('n_units_right')):d['n_units_total']=d['n_units_left']+d['n_units_right']
    if finite(d.get('n_control')) and finite(d.get('n_rightstim')):d['right_trials_total']=d['n_control']+d['n_rightstim']
    mice.append(d)
def perm(x,y,nperm=NPERM):
    x=np.asarray(x,float);y=np.asarray(y,float);obs=float(spearmanr(x,y).statistic);rx=rankdata(x);ry=rankdata(y);rx-=rx.mean();den=np.sqrt(np.sum(rx*rx));rng=np.random.default_rng(SEED);ge=0;done=0
    while done<nperm:
        b=min(10000,nperm-done);idx=np.argsort(rng.random((b,len(y))),axis=1);yp=ry[idx];yp=yp-yp.mean(1,keepdims=True);v=(yp@rx)/(den*np.sqrt(np.sum(yp*yp,axis=1)));ge+=int(np.sum(v>=obs-1e-12));done+=b
    return obs,float((ge+1)/(nperm+1))
def boot(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);rng=np.random.default_rng(SEED);v=[]
    for _ in range(NBOOT):
        idx=rng.integers(0,len(x),len(x))
        if len(np.unique(x[idx]))<2 or len(np.unique(y[idx]))<2:continue
        z=spearmanr(x[idx],y[idx]).statistic
        if np.isfinite(z):v.append(float(z))
    return [float(np.percentile(v,2.5)),float(np.percentile(v,97.5))] if v else [None,None]
def loo(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);pred=np.full(len(y),np.nan)
    if len(x)<4:return None,None,[]
    for i in range(len(y)):
        tr=np.arange(len(y))!=i;A=np.column_stack([np.ones(tr.sum()),x[tr]]);b=np.linalg.lstsq(A,y[tr],rcond=None)[0];pred[i]=b[0]+b[1]*x[i]
    rho=float(spearmanr(pred,y).statistic);var=np.mean((y-y.mean())**2);r2=float(1-np.mean((pred-y)**2)/var) if var>0 else None;return rho,r2,pred.tolist()
def assoc(ykey,xkey='control_bilateral_coupling'):
    q=[m for m in mice if finite(m.get(xkey)) and finite(m.get(ykey))];x=[m[xkey] for m in q];y=[m[ykey] for m in q]
    if len(q)<3:return {'n_mice':len(q)}
    rho,p=perm(x,y);lr,lr2,_=loo(x,y);ts=theilslopes(y,x,alpha=.95)
    return {'n_mice':len(q),'spearman_rho':rho,'one_sided_permutation_p':p,'bootstrap_95_ci':boot(x,y),'loo_spearman_rho':lr,'loo_r2':lr2,'theil_sen_slope':float(ts.slope),'theil_sen_95_ci':[float(ts.low_slope),float(ts.high_slope)]}
def partial_rank(ykey,covkeys,xkey='control_bilateral_coupling',nperm=200000):
    q=[m for m in mice if finite(m.get(xkey)) and finite(m.get(ykey)) and all(finite(m.get(k)) for k in covkeys)]
    if len(q)<5:return {'n_mice':len(q)}
    rx=rankdata([m[xkey] for m in q]).astype(float);ry=rankdata([m[ykey] for m in q]).astype(float);C=np.column_stack([rankdata([m[k] for m in q]) for k in covkeys]);Xd=np.column_stack([np.ones(len(q)),C]);M=np.eye(len(q))-Xd@np.linalg.pinv(Xd);ex=M@rx;ey=M@ry;obs=float((ex@ey)/np.sqrt((ex@ex)*(ey@ey)));rng=np.random.default_rng(SEED+17);ge=0;done=0
    while done<nperm:
        b=min(5000,nperm-done);idx=np.argsort(rng.random((b,len(q))),axis=1);py=ry[idx];epy=py@M.T;vals=(epy@ex)/(np.sqrt(ex@ex)*np.sqrt(np.sum(epy*epy,axis=1)));ge+=int(np.sum(vals>=obs-1e-12));done+=b
    return {'n_mice':len(q),'partial_rank_r':obs,'one_sided_permutation_p':float((ge+1)/(nperm+1)),'covariates':covkeys}
def jackknife(ykey,xkey='control_bilateral_coupling'):
    q=[m for m in mice if finite(m.get(xkey)) and finite(m.get(ykey))];out=[]
    for i,m in enumerate(q):
        z=[u for j,u in enumerate(q) if j!=i];r=spearmanr([u[xkey] for u in z],[u[ykey] for u in z]);out.append({'omitted_mouse':m['animal'],'rho':float(r.statistic),'two_sided_asymptotic_p':float(r.pvalue)})
    return {'min_rho':float(min(z['rho'] for z in out)),'max_rho':float(max(z['rho'] for z in out)),'all_positive':all(z['rho']>0 for z in out),'all_two_sided_p_lt_0_05':all(z['two_sided_asymptotic_p']<.05 for z in out),'rows':out}
def holm(pdict):
    items=sorted(pdict.items(),key=lambda x:x[1]);m=len(items);out={};prev=0
    for i,(k,p) in enumerate(items):a=max(prev,min(1,(m-i)*p));out[k]=a;prev=a
    return out
primary4={k:assoc(k) for k in ['right_behavior_recovery','right_neuronal_recovery','left_behavior_recovery','left_neuronal_recovery']};hadj=holm({k:v.get('one_sided_permutation_p',1) for k,v in primary4.items()})
for k in primary4:primary4[k]['holm_adjusted_p_across_four_side_endpoint_tests']=hadj[k]
result={'analysis_status':'POST_HOC_RAW_CHEN_ROBUSTNESS_FROM_SESSION_SPIKES','dataset':'Chen et al. Cell 2021 / Zenodo 6713616','n_raw_sessions_parsed':len(rows),'n_mice_any':len(mice),'original_side_endpoint_tests':primary4,'right_behavior_additive_endpoint':assoc('right_behavior_delta'),'right_behavior_fixed10_coupling_ratio':assoc('right_behavior_recovery','control_bilateral_coupling_fixed10'),'right_behavior_fixed10_coupling_delta':assoc('right_behavior_delta','control_bilateral_coupling_fixed10'),'right_behavior_fixed20_coupling_ratio':assoc('right_behavior_recovery','control_bilateral_coupling_fixed20'),'right_behavior_fixed20_coupling_delta':assoc('right_behavior_delta','control_bilateral_coupling_fixed20'),'jackknife_right_behavior_original':jackknife('right_behavior_recovery'),'partial_right_behavior_adjust_total_units':partial_rank('right_behavior_recovery',['n_units_total']),'partial_right_behavior_adjust_trial_count':partial_rank('right_behavior_recovery',['right_trials_total']),'partial_right_behavior_adjust_control_accuracy':partial_rank('right_behavior_recovery',['control_behavior_accuracy']),'partial_right_behavior_adjust_local_auc':partial_rank('right_behavior_recovery',['transferable_local_cv_auc']),'partial_right_behavior_adjust_units_and_control_accuracy':partial_rank('right_behavior_recovery',['n_units_total','control_behavior_accuracy']),'full_coupling_vs_unit_count':assoc('n_units_total'),'fixed20_coupling_vs_unit_count':assoc('n_units_total','control_bilateral_coupling_fixed20'),'full_coupling_vs_control_accuracy':assoc('control_behavior_accuracy'),'note':'All analyses beyond the original locked raw Chen analysis are explicitly post hoc robustness analyses. Fixed-unit coupling controls population yield; additive endpoints test ratio-denominator sensitivity; partial-rank sensitivities test unit yield, trial count, control accuracy and local code strength.'}
(OUT/'result.json').write_text(json.dumps(result,indent=2));(OUT/'mouse_rows.json').write_text(json.dumps(mice,indent=2));(OUT/'session_rows.json').write_text(json.dumps(rows,indent=2))
with open(OUT/'mouse_rows.csv','w',newline='') as f:
    if mice:w=csv.DictWriter(f,fieldnames=list(mice[0].keys()));w.writeheader();w.writerows(mice)
print(json.dumps(result,indent=2))
