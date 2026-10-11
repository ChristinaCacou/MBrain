import ast, hashlib, json
from pathlib import Path
from urllib.request import urlopen
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LinearRegression

SEED=20261011; N_PERM=500000; N_BOOT=50000; BATCH=10000
ROOT=Path('analysis/pui_blind_fixedC'); FEATURE=ROOT/'blind_features_region.csv'; LOCK=ROOT/'blind_lock.json'; OUT=ROOT/'reveal_vectorized'; OUT.mkdir(parents=True,exist_ok=True)
COMMIT='a17677ef2ab38c198761ad66f8a8548e4b857901'; URL=f'https://raw.githubusercontent.com/elliottower/neural-geometry-reliability/{COMMIT}/experiments/shared_bundle.py'

def sha(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def outcomes():
 t=urlopen(URL).read().decode(); tree=ast.parse(t)
 for node in tree.body:
  if isinstance(node,ast.Assign) and any(isinstance(z,ast.Name) and z.id=='SILENCING_EFFECTS' for z in node.targets):
   return {str(k):float(v) for k,v in ast.literal_eval(node.value).items()},hashlib.sha256(t.encode()).hexdigest()
 raise RuntimeError('missing outcomes')

def corr_ranked(rx,ry):
 rx=np.asarray(rx,float); ry=np.asarray(ry,float); a=rx-rx.mean(); b=ry-ry.mean(); return float(a@b/np.sqrt((a@a)*(b@b)))

def observed(x,y): return corr_ranked(rankdata(x),rankdata(y))

def permutation_stats(x,y,rng,n=N_PERM):
 rx=rankdata(np.asarray(x,float)); ry=rankdata(np.asarray(y,float)); obs= corr_ranked(rx,ry); xc=rx-rx.mean(); denx=np.sqrt(xc@xc); one=two=0; done=0
 while done<n:
  k=min(BATCH,n-done); keys=rng.random((k,len(ry))); order=np.argsort(keys,axis=1); P=ry[order]; Pc=P-P.mean(axis=1,keepdims=True); r=(Pc@xc)/(np.sqrt(np.sum(Pc*Pc,axis=1))*denx); one += int(np.sum(r>=obs-1e-15)); two += int(np.sum(np.abs(r)>=abs(obs)-1e-15)); done+=k
 return obs,(one+1)/(n+1),(two+1)/(n+1)

def bootstrap_ci(x,y,rng,n=N_BOOT):
 x=np.asarray(x,float); y=np.asarray(y,float); vals=[]; done=0
 while done<n:
  k=min(BATCH,n-done); idx=rng.integers(0,len(x),size=(k,len(x))); X=x[idx]; Y=y[idx]; RX=rankdata(X,axis=1); RY=rankdata(Y,axis=1); AX=RX-RX.mean(axis=1,keepdims=True); AY=RY-RY.mean(axis=1,keepdims=True); den=np.sqrt(np.sum(AX*AX,axis=1)*np.sum(AY*AY,axis=1)); ok=den>0; vals.append(np.sum(AX[ok]*AY[ok],axis=1)/den[ok]); done+=k
 v=np.concatenate(vals); return [float(np.percentile(v,2.5)),float(np.percentile(v,97.5))]

def delta_perm(x1,x2,y,obsdelta,rng,n=N_PERM):
 r1=rankdata(x1); r2=rankdata(x2); ry=rankdata(y); a1=r1-r1.mean(); a2=r2-r2.mean(); d1=np.sqrt(a1@a1); d2=np.sqrt(a2@a2); count=done=0
 while done<n:
  k=min(BATCH,n-done); P=ry[np.argsort(rng.random((k,len(ry))),axis=1)]; P-=P.mean(axis=1,keepdims=True); dy=np.sqrt(np.sum(P*P,axis=1)); q=(P@a1)/(dy*d1)-(P@a2)/(dy*d2); count += int(np.sum(q>=obsdelta-1e-15)); done+=k
 return (count+1)/(n+1)

def partial(x,y,controls):
 xr=rankdata(x); yr=rankdata(y); Z=np.column_stack([rankdata(z) for z in controls]); A=np.column_stack([np.ones(len(xr)),Z]); rx=xr-A@np.linalg.lstsq(A,xr,rcond=None)[0]; ry=yr-A@np.linalg.lstsq(A,yr,rcond=None)[0]; return float(spearmanr(rx,ry).statistic),float(spearmanr(rx,ry).pvalue)

def loo(x,y):
 x=np.asarray(x,float).reshape(-1,1); y=np.asarray(y,float); pred=np.zeros(len(y))
 for i in range(len(y)):
  k=np.arange(len(y))!=i; m=LinearRegression().fit(x[k],y[k]); pred[i]=m.predict(x[i:i+1])[0]
 rho,p=spearmanr(pred,y); mse=float(np.mean((pred-y)**2)); var=float(np.mean((y-y.mean())**2)); return float(rho),float(p),float(1-mse/var),mse

def main():
 lock=json.loads(LOCK.read_text()); actual=sha(FEATURE); exp=lock['feature_sha256']['blind_features_region.csv'];
 if actual!=exp or lock.get('perturbation_labels_loaded') is not False: raise RuntimeError('blind lock violation')
 f=pd.read_csv(FEATURE); o,osha=outcomes(); od=pd.DataFrame({'region':list(o),'silencing_effect':list(o.values())}); m=od.merge(f,on='region').sort_values('region').reset_index(drop=True)
 if len(m)!=len(o): raise RuntimeError('outcome overlap incomplete')
 y=m.silencing_effect.to_numpy(float); metrics=['pui_mean','pui_median','marginal_mean','local_bacc_mean','local_auc_mean','snr_mean','mean_window_spikes','neurons_recorded_mean','n_sessions']; labels=['PUI (primary)','PUI median','marginal information','local decoder balanced accuracy','local decoder AUC','condition SNR','mean window spikes','neuron yield','session coverage']; rows=[]
 for j,(key,label) in enumerate(zip(metrics,labels)):
  x=m[key].to_numpy(float); ok=np.isfinite(x); xx=x[ok]; yy=y[ok]; rho,po,pt=permutation_stats(xx,yy,np.random.default_rng(SEED+100*j)); ci=bootstrap_ci(xx,yy,np.random.default_rng(SEED+100*j+1)); lr,lp,r2,mse=loo(xx,yy)
  rows.append(dict(metric=key,label=label,n=len(xx),spearman_rho=rho,permutation_p_one_sided_positive=po,permutation_p_two_sided=pt,bootstrap_ci_low=ci[0],bootstrap_ci_high=ci[1],loo_spearman_rho=lr,loo_spearman_p=lp,loo_r2=r2,loo_mse=mse))
 r=pd.DataFrame(rows); pr=float(r.loc[r.metric=='pui_mean','spearman_rho'].iloc[0]); mr=float(r.loc[r.metric=='marginal_mean','spearman_rho'].iloc[0]); delta=pr-mr; dp=delta_perm(m.pui_mean.to_numpy(),m.marginal_mean.to_numpy(),y,delta,np.random.default_rng(SEED+999)); par,parp=partial(m.pui_mean.to_numpy(),y,[m.neurons_recorded_mean.to_numpy(),m.snr_mean.to_numpy()]); loo_r=[]
 for i in range(len(m)):
  k=np.arange(len(m))!=i; loo_r.append(observed(m.pui_mean.to_numpy()[k],y[k]))
 primary={'feature_sha256_verified':actual,'blind_lock_script_sha256':lock['script_sha256'],'outcome_source_commit':COMMIT,'outcome_source_sha256':osha,'n_regions':len(m),'regions':m.region.tolist(),'primary_metric':'pui_mean','primary_spearman_rho':pr,'primary_permutation_p_one_sided_positive':float(r.loc[r.metric=='pui_mean','permutation_p_one_sided_positive'].iloc[0]),'primary_permutation_p_two_sided':float(r.loc[r.metric=='pui_mean','permutation_p_two_sided'].iloc[0]),'primary_bootstrap_95_ci':[float(r.loc[r.metric=='pui_mean','bootstrap_ci_low'].iloc[0]),float(r.loc[r.metric=='pui_mean','bootstrap_ci_high'].iloc[0])],'partial_rank_rho_controlling_neuron_yield_and_snr':par,'partial_rank_p':parp,'marginal_rho':mr,'delta_rho_pui_minus_marginal':delta,'delta_rho_permutation_p_one_sided':float(dp),'leave_one_region_out_rho_min':float(np.min(loo_r)),'leave_one_region_out_rho_median':float(np.median(loo_r)),'leave_one_region_out_rho_max':float(np.max(loo_r)),'success_criterion':'positive primary association with one-sided permutation p<0.05 AND PUI rho > marginal rho'}; primary['success']=bool(primary['primary_permutation_p_one_sided_positive']<.05 and delta>0)
 m.to_csv(OUT/'revealed_region_table.csv',index=False); r.to_csv(OUT/'metric_comparison.csv',index=False); (OUT/'primary_result.json').write_text(json.dumps(primary,indent=2,sort_keys=True)); (OUT/'RESULT.md').write_text(f"# Exact blinded causal-prediction test\n\n**Result: {'PASS' if primary['success'] else 'FAIL'}**\n\nPrimary PUI rho = {pr:.6f}; one-sided permutation p = {primary['primary_permutation_p_one_sided_positive']:.6g}; two-sided p = {primary['primary_permutation_p_two_sided']:.6g}.\n\nBootstrap 95% CI = [{primary['primary_bootstrap_95_ci'][0]:.4f}, {primary['primary_bootstrap_95_ci'][1]:.4f}].\n\nMarginal rho = {mr:.6f}; delta = {delta:.6f}; delta permutation p = {dp:.6g}.\n\nPartial rho controlling neuron yield and SNR = {par:.6f}.\n\nFrozen feature SHA-256: `{actual}`.\n")
 print(json.dumps(primary,indent=2)); print(r.to_string(index=False))
if __name__=='__main__': main()
