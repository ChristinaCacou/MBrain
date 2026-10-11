import hashlib, json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, log_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

SEED=20261011
N_OUTER=4
N_DRAWS=3
FIXED_NEURONS=5
TARGET_PCS=2
OTHER_PCS=1
C_FIXED=1.0
WINDOW_START=0.0
WINDOW_END=0.2
MIN_RESPONDED_TRIALS=80
MIN_CLASS=20
URLS={
 'steinmetz_part0.npz':'https://raw.githubusercontent.com/tipsyfermion/NMA_Mapping_Brain_Networks_2020/master/code/steinmetz_part0.npz',
 'steinmetz_part1.npz':'https://raw.githubusercontent.com/tipsyfermion/NMA_Mapping_Brain_Networks_2020/master/code/steinmetz_part1.npz',
 'steinmetz_part2.npz':'https://raw.githubusercontent.com/tipsyfermion/NMA_Mapping_Brain_Networks_2020/master/code/steinmetz_part2.npz'}
OUT=Path('analysis/pui_blind_fixedC'); CACHE=Path('.cache/steinmetz')
OUT.mkdir(parents=True,exist_ok=True); CACHE.mkdir(parents=True,exist_ok=True)

def hsh(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def alias(a):
 a=str(a).strip()
 if a.startswith('SSp'): return 'SSp'
 if a.startswith('ACA'): return 'ACA'
 if a.startswith('ORB'): return 'ORB'
 return a

def load():
 ss=[]; sums={}
 for name,url in URLS.items():
  p=CACHE/name
  if not p.exists(): urlretrieve(url,p)
  sums[name]=hsh(p)
  with np.load(p,allow_pickle=True) as z: part=z['dat'].tolist()
  ss.extend(part); print(name,len(part),flush=True)
 if len(ss)!=39: raise RuntimeError(len(ss))
 return ss,sums

def win(dat):
 dt=float(np.asarray(dat.get('bin_size',.01)).reshape(-1)[0]); t0=float(np.asarray(dat.get('stim_onset',.5)).reshape(-1)[0])
 n=np.asarray(dat['spks']).shape[2]
 a=max(0,min(n-1,int(round((t0+WINDOW_START)/dt)))); b=max(a+1,min(n,int(round((t0+WINDOW_END)/dt))))
 return slice(a,b),dt

def eligible(dat,n):
 r=np.asarray(dat['response'],float).reshape(-1); ok=np.isfinite(r)&(r!=0)
 try:
  rt=np.asarray(dat['response_time'],float).reshape(-1)
  if len(rt)==n: ok &= np.isfinite(rt)&(rt>WINDOW_END)
 except Exception: pass
 idx=np.flatnonzero(ok); rr=r[idx]; vals=np.unique(rr)
 if vals.size!=2: return np.array([],int),np.array([],int)
 return idx,(rr==vals.max()).astype(int)

def Hmat(dat,n):
 cl=np.asarray(dat['contrast_left'],float).reshape(-1); cr=np.asarray(dat['contrast_right'],float).reshape(-1); rsp=np.asarray(dat['response'],float).reshape(-1)
 fb=np.asarray(dat.get('feedback_type',np.zeros(n)),float).reshape(-1)
 H=np.column_stack([cl,cr,cr-cl,cr+cl]); prev=np.zeros((n,4))
 prev[1:,0]=rsp[:-1]
 if len(fb)==n: prev[1:,1]=fb[:-1]
 prev[1:,2]=cl[:-1]; prev[1:,3]=cr[:-1]
 return np.column_stack([H,prev])

def project(X,tr,te,k):
 sc=StandardScaler(); A=sc.fit_transform(X[tr]); B=sc.transform(X[te]); kk=min(k,A.shape[1],max(1,A.shape[0]-1))
 pc=PCA(n_components=kk,svd_solver='full',random_state=SEED); return pc.fit_transform(A),pc.transform(B)

def eval_lr(A,y,B,z):
 sc=StandardScaler(); At=sc.fit_transform(A); Bt=sc.transform(B)
 m=LogisticRegression(C=C_FIXED,solver='lbfgs',max_iter=1000,random_state=SEED); m.fit(At,y); p=m.predict_proba(Bt)[:,1]; pred=(p>=.5).astype(int)
 return float(log_loss(z,p,labels=[0,1])),float(balanced_accuracy_score(z,pred)),float(roc_auc_score(z,p)) if np.unique(z).size==2 else np.nan

def snr(signal,cl,cr):
 keys=list(zip(cl.tolist(),cr.tolist())); means=[]; within=[]
 for k in sorted(set(keys)):
  ii=np.array([i for i,x in enumerate(keys) if x==k]);
  if len(ii)>=2: means.append(np.mean(signal[ii])); within.append(np.var(signal[ii],ddof=1))
 return float(np.var(means,ddof=1)/(np.mean(within)+1e-12)) if len(means)>=2 else np.nan

def process(si,dat):
 sp=np.asarray(dat['spks'],np.float32); nn,nt,_=sp.shape; idx,y=eligible(dat,nt)
 if len(idx)<MIN_RESPONDED_TRIALS or min(np.bincount(y,minlength=2))<MIN_CLASS: return []
 H=Hmat(dat,nt)[idx]; sl,dt=win(dat); raw=np.asarray(dat['brain_area']).astype(str); aa=np.asarray([alias(x) for x in raw])
 regs={a:np.flatnonzero(aa==a) for a in sorted(set(aa)) if a.lower() not in {'root','nan','','none'} and np.sum(aa==a)>=FIXED_NEURONS}
 if len(regs)<2:return []
 acts={a:sp[ii][:,idx,sl].sum(axis=2).T.astype(float) for a,ii in regs.items()}
 cl=np.asarray(dat['contrast_left'],float).reshape(-1)[idx]; cr=np.asarray(dat['contrast_right'],float).reshape(-1)[idx]
 mouse=str(np.asarray(dat.get('mouse_name',f'session{si}')).reshape(-1)[0]); date=str(np.asarray(dat.get('date_exp','')).reshape(-1)[0])
 folds=list(StratifiedKFold(n_splits=N_OUTER,shuffle=True,random_state=SEED+si).split(H,y)); out=[]
 for d in range(N_DRAWS):
  rng=np.random.default_rng(SEED+si*1009+d*9176); sampled={a:acts[a][:,rng.choice(acts[a].shape[1],FIXED_NEURONS,replace=False)] for a in regs}
  foldcache=[]
  for fi,(tr,te) in enumerate(folds):
   projections={}
   for a,X in sampled.items():
    # two PCs retained for target use; first PC is used for conditioning when this is another region
    projections[a]=project(X,tr,te,TARGET_PCS)
   hres=eval_lr(H[tr],y[tr],H[te],y[te])
   foldcache.append((tr,te,projections,hres))
  for target in sorted(sampled):
   st=[]
   for tr,te,P,hres in foldcache:
    Ttr,Tte=P[target]
    otr=[]; ote=[]
    for a in sorted(P):
     if a==target: continue
     otr.append(P[a][0][:,:OTHER_PCS]); ote.append(P[a][1][:,:OTHER_PCS])
    Otr=np.column_stack(otr); Ote=np.column_stack(ote)
    baseA=np.column_stack([H[tr],Otr]); baseB=np.column_stack([H[te],Ote])
    fullA=np.column_stack([H[tr],Otr,Ttr]); fullB=np.column_stack([H[te],Ote,Tte])
    margA=np.column_stack([H[tr],Ttr]); margB=np.column_stack([H[te],Tte])
    bl=eval_lr(baseA,y[tr],baseB,y[te])[0]; fl=eval_lr(fullA,y[tr],fullB,y[te])[0]; ml=eval_lr(margA,y[tr],margB,y[te])[0]; tl,tb,ta=eval_lr(Ttr,y[tr],Tte,y[te]); hl=hres[0]
    st.append([hl,bl,fl,ml,tl,tb,ta])
   m=np.nanmean(np.asarray(st),axis=0); signal=sampled[target].mean(axis=1)
   out.append(dict(session=si,mouse=mouse,date=date,region=target,draw=d,n_trials=len(idx),n_regions_conditioned=len(sampled)-1,n_neurons_recorded=len(regs[target]),n_neurons_matched=FIXED_NEURONS,bin_size_s=dt,window_start_s=WINDOW_START,window_end_s=WINDOW_END,pui_nats_per_trial=m[1]-m[2],marginal_gain_nats_per_trial=m[0]-m[3],target_only_loss=m[4],target_only_bacc=m[5],target_only_auc=m[6],condition_snr=snr(signal,cl,cr),mean_window_spikes=float(np.mean(signal))))
 return out

def main():
 sessions,raws=load(); rows=[]
 for i,d in enumerate(sessions):
  rr=process(i,d); rows+=rr; print(f'{i+1}/39 rows={len(rr)} total={len(rows)}',flush=True)
 df=pd.DataFrame(rows).sort_values(['region','session','draw']).reset_index(drop=True)
 sess=(df.groupby(['region','session'],as_index=False).agg(pui_nats_per_trial=('pui_nats_per_trial','mean'),marginal_gain_nats_per_trial=('marginal_gain_nats_per_trial','mean'),target_only_bacc=('target_only_bacc','mean'),target_only_auc=('target_only_auc','mean'),condition_snr=('condition_snr','mean'),mean_window_spikes=('mean_window_spikes','mean'),n_neurons_recorded=('n_neurons_recorded','first'),n_trials=('n_trials','first'),n_regions_conditioned=('n_regions_conditioned','first'),mouse=('mouse','first')))
 reg=(sess.groupby('region',as_index=False).agg(pui_mean=('pui_nats_per_trial','mean'),pui_median=('pui_nats_per_trial','median'),pui_sd=('pui_nats_per_trial','std'),marginal_mean=('marginal_gain_nats_per_trial','mean'),local_bacc_mean=('target_only_bacc','mean'),local_auc_mean=('target_only_auc','mean'),snr_mean=('condition_snr','mean'),mean_window_spikes=('mean_window_spikes','mean'),neurons_recorded_mean=('n_neurons_recorded','mean'),conditioned_regions_mean=('n_regions_conditioned','mean'),n_sessions=('session','nunique'),n_mice=('mouse','nunique')))
 files={'blind_features_draws.csv':df,'blind_features_session.csv':sess,'blind_features_region.csv':reg}
 for n,x in files.items(): x.to_csv(OUT/n,index=False)
 config={'seed':SEED,'outcome_blind':True,'perturbation_labels_loaded':False,'primary_metric':'pui_mean','target':'responded-trial choice direction','window_s':[WINDOW_START,WINDOW_END],'response_before_window_excluded':True,'fixed_neurons_per_region':FIXED_NEURONS,'target_pcs':TARGET_PCS,'other_region_pcs':OTHER_PCS,'outer_folds':N_OUTER,'draws':N_DRAWS,'logistic_C_fixed':C_FIXED,'history_covariates':['contrast_left','contrast_right','signed_evidence','total_contrast','prev_response','prev_feedback','prev_contrast_left','prev_contrast_right'],'raw_sha256':raws,'script_sha256':hsh(Path(__file__))}
 config['feature_sha256']={n:hsh(OUT/n) for n in files}
 (OUT/'blind_lock.json').write_text(json.dumps(config,indent=2,sort_keys=True))
 (OUT/'BLIND_LOCK.md').write_text('# Blinded PUI feature lock\n\nNo perturbation outcomes were loaded or referenced in this phase.\n\nPrimary predictor: `pui_mean`.\n\n```json\n'+json.dumps(config['feature_sha256'],indent=2,sort_keys=True)+'\n```\n')
 print(reg.sort_values('pui_mean',ascending=False).head(30).to_string(index=False),flush=True)
if __name__=='__main__': main()
