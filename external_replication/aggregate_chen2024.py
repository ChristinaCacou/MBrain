import json, os
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
from sklearn.linear_model import LinearRegression

SEED=20261012; NPERM=500000; NBOOT=50000; BATCH=10000
ROOT=Path(os.environ.get('INPUT_DIR','all-session-results'))
OUT=Path('analysis/chen2024_external_final'); OUT.mkdir(parents=True,exist_ok=True)

def corr(x,y):
    r=spearmanr(x,y).statistic
    return float(r) if np.isfinite(r) else float('nan')

def rankcorr_rows(X,Y):
    RX=rankdata(X,axis=1); RY=rankdata(Y,axis=1)
    AX=RX-RX.mean(axis=1,keepdims=True); AY=RY-RY.mean(axis=1,keepdims=True)
    den=np.sqrt(np.sum(AX*AX,axis=1)*np.sum(AY*AY,axis=1))
    out=np.full(len(X),np.nan); ok=den>0; out[ok]=np.sum(AX[ok]*AY[ok],axis=1)/den[ok]; return out

def inference(df, outcome, seed):
    d=df[np.isfinite(df['coupling']) & np.isfinite(df[outcome])].copy()
    x=d['coupling'].to_numpy(float); y=d[outcome].to_numpy(float); n=len(d)
    if n<3:return {'n_mice':n,'mice':d.subject.tolist(),'rho':None,'p_one_sided':None,'bootstrap95':None,'loo_rho':None,'loo_r2':None,'loo_mae':None}
    obs=corr(x,y)
    rx=rankdata(x); ry=rankdata(y); xc=rx-rx.mean(); denx=np.sqrt(xc@xc)
    rng=np.random.default_rng(seed); ge=done=0
    while done<NPERM:
        k=min(BATCH,NPERM-done); order=np.argsort(rng.random((k,n)),axis=1); P=ry[order]; P-=P.mean(axis=1,keepdims=True)
        r=(P@xc)/(np.sqrt(np.sum(P*P,axis=1))*denx); ge+=int(np.sum(r>=obs-1e-15)); done+=k
    p=(ge+1)/(NPERM+1)
    rng=np.random.default_rng(seed+1); boots=[]; done=0
    while done<NBOOT:
        k=min(BATCH,NBOOT-done); idx=rng.integers(0,n,size=(k,n)); r=rankcorr_rows(x[idx],y[idx]); boots.append(r[np.isfinite(r)]);done+=k
    bv=np.concatenate(boots) if boots else np.array([]); ci=[float(np.percentile(bv,2.5)),float(np.percentile(bv,97.5))] if len(bv) else [None,None]
    pred=np.zeros(n)
    for i in range(n):
        tr=np.arange(n)!=i; m=LinearRegression().fit(x[tr,None],y[tr]);pred[i]=m.predict([[x[i]]])[0]
    loo_r=corr(pred,y); mse=float(np.mean((pred-y)**2)); var=float(np.mean((y-y.mean())**2)); r2=float(1-mse/var) if var>0 else None; mae=float(np.mean(np.abs(pred-y)))
    return {'n_mice':n,'mice':d.subject.tolist(),'rho':obs,'p_one_sided':float(p),'bootstrap95':ci,'loo_rho':loo_r,'loo_r2':r2,'loo_mae':mae}

def main():
    files=list(ROOT.rglob('*.json')); rows=[]
    for fp in files:
        try:r=json.loads(fp.read_text())
        except:continue
        if 'control_bilateral_coupling' not in r:continue
        rows.append({'asset_id':r.get('asset_id'),'subject':r.get('subject'),'identifier':r.get('identifier'),'coupling':r.get('control_bilateral_coupling'),'within_choice_coupling':r.get('control_bilateral_within_choice_coupling'),'left_recovery':(r.get('behavior_recovery') or {}).get('left'),'right_recovery':(r.get('behavior_recovery') or {}).get('right'),'side_average_recovery':(r.get('behavior_recovery') or {}).get('side_average'),'left_neuronal':(r.get('neuronal_recovery') or {}).get('left'),'right_neuronal':(r.get('neuronal_recovery') or {}).get('right'),'n_left_units':r.get('n_good_units_left_ALM'),'n_right_units':r.get('n_good_units_right_ALM'),'control_trials':r.get('control_trial_count'),'left_stim_trials':(r.get('stim_trial_counts') or {}).get('left',0),'right_stim_trials':(r.get('stim_trial_counts') or {}).get('right',0),'bilateral_stim_trials':(r.get('stim_trial_counts') or {}).get('bilateral',0)})
    s=pd.DataFrame(rows)
    if s.empty: raise RuntimeError('no session result JSON files found')
    for c in ['coupling','within_choice_coupling','left_recovery','right_recovery','side_average_recovery','left_neuronal','right_neuronal']:
        s[c]=pd.to_numeric(s[c],errors='coerce')
    s.to_csv(OUT/'session_table.csv',index=False)
    mouse_base=sorted(s.subject.dropna().unique()); mrows=[]
    for subj in mouse_base:
        q=s[s.subject==subj]; z={'subject':subj,'n_sessions_total':len(q)}
        for ep in ['side_average_recovery','left_recovery','right_recovery','left_neuronal','right_neuronal']:
            qq=q[np.isfinite(q.coupling)&np.isfinite(q[ep])]
            z[ep]=float(qq[ep].mean()) if len(qq) else np.nan; z[f'coupling_for_{ep}']=float(qq.coupling.mean()) if len(qq) else np.nan; z[f'within_for_{ep}']=float(qq.within_choice_coupling.mean()) if len(qq) and qq.within_choice_coupling.notna().any() else np.nan; z[f'n_sessions_{ep}']=int(len(qq))
        mrows.append(z)
    wide=pd.DataFrame(mrows); wide.to_csv(OUT/'mouse_endpoint_table.csv',index=False)
    results={}
    eps=['side_average_recovery','left_recovery','right_recovery','left_neuronal','right_neuronal']
    for i,ep in enumerate(eps):
        d=pd.DataFrame({'subject':wide.subject,'coupling':wide[f'coupling_for_{ep}'],ep:wide[ep]});results[ep]=inference(d,ep,SEED+100*i)
    controls={}
    for i,ep in enumerate(['left_recovery','right_recovery']):
        d=pd.DataFrame({'subject':wide.subject,'coupling':wide[f'within_for_{ep}'],ep:wide[ep]});controls[ep]=inference(d,ep,SEED+1000+100*i)
    for ep in ['left_recovery','right_recovery']:
        p=results[ep].get('p_one_sided');results[ep]['bonferroni_two_side_secondary']=min(1.0,2*p) if p is not None else None
    primary=results['side_average_recovery'];primary['locked_success']=bool(primary.get('n_mice',0)>=8 and primary.get('rho') is not None and primary['rho']>0 and primary.get('p_one_sided',1)<.05 and primary.get('loo_r2') is not None and primary['loo_r2']>0);primary['status']=('STRUCTURALLY_UNDERPOWERED_OR_INELIGIBLE_UNDER_FROZEN_PRIMARY' if primary.get('n_mice',0)<8 else ('PASS' if primary['locked_success'] else 'FAIL'))
    result={'analysis_status':'FROZEN_EXTERNAL_CHEN2024_BILATERAL_COUPLING_TEST','dataset':'Chen et al. Cell 2024 / DANDI 000363 v0.230822.0128','n_session_json':len(s),'n_subjects_any_session':len(mouse_base),'primary_endpoint':'side_average_recovery','primary':primary,'secondary':{k:results[k] for k in ['left_recovery','right_recovery','left_neuronal','right_neuronal']},'negative_control_within_choice':controls,'lock':'external_replication/CHEN2024_EXTERNAL_LOCK.md','note':'Primary and secondary definitions were frozen before aggregate perturbation-outcome extraction. A schema-validation session was run after the lock and before this aggregate analysis.'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True));print(json.dumps(result,indent=2,sort_keys=True));print(s.to_string(index=False));print(wide.to_string(index=False))
if __name__=='__main__':main()
