import json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
IN=Path('analysis_inputs/chen_mouse_rows.csv'); OUT=Path('analysis/chen_specificity'); OUT.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(IN); x=df.coupling.to_numpy(float); SEED=20261012

def rowcorr(A,B):
 A=rankdata(A,axis=1); B=rankdata(B,axis=1); A=A-A.mean(axis=1,keepdims=True); B=B-B.mean(axis=1,keepdims=True)
 return np.sum(A*B,axis=1)/np.sqrt(np.sum(A*A,axis=1)*np.sum(B*B,axis=1))
def rho_vec(xv,yv):
 X=rankdata(xv); Y=rankdata(yv); X=X-X.mean(); Y=Y-Y.mean(); return float(np.dot(X,Y)/np.sqrt(np.dot(X,X)*np.dot(Y,Y)))
def exact_swap(a,b):
 a=np.asarray(a,float); b=np.asarray(b,float); n=len(a); obs=rho_vec(x,a)-rho_vec(x,b); ge=tw=0
 batch=8192; total=1<<n
 for start in range(0,total,batch):
  masks=np.arange(start,min(start+batch,total),dtype=np.uint32)[:,None]; bits=((masks>>np.arange(n,dtype=np.uint32))&1).astype(bool)
  A=np.where(bits,b[None,:],a[None,:]); B=np.where(bits,a[None,:],b[None,:]); X=np.broadcast_to(x,A.shape)
  d=rowcorr(X,A)-rowcorr(X,B); ge+=int(np.sum(d>=obs-1e-15)); tw+=int(np.sum(np.abs(d)>=abs(obs)-1e-15))
 rng=np.random.default_rng(SEED); boots=[]
 for _ in range(10):
  idx=rng.integers(0,n,size=(5000,n)); X=x[idx]; A=a[idx]; B=b[idx]; boots.append(rowcorr(X,A)-rowcorr(X,B))
 boots=np.concatenate(boots)
 return {'delta_rho':obs,'exact_one_sided_p':ge/total,'exact_two_sided_p':tw/total,'bootstrap_95_ci':[float(np.nanpercentile(boots,2.5)),float(np.nanpercentile(boots,97.5))],'n_mice':n,'bootstrap_resamples':50000}
res={
 'right_behavior_vs_left_behavior':exact_swap(df.right_behavior_recovery,df.left_behavior_recovery),
 'right_behavior_vs_right_neuronal':exact_swap(df.right_behavior_recovery,df.right_neuronal_recovery),
 'right_neuronal_vs_left_neuronal':exact_swap(df.right_neuronal_recovery,df.left_neuronal_recovery),
 'status':'POST_HOC_PAIRED_SPECIFICITY_AUDIT',
 'interpretation':'Post hoc paired specificity contrasts. A significant association in one endpoint and nonsignificant association in another does not itself establish a significant difference; these exact within-mouse swap tests directly test the correlation differences.'
}
(OUT/'paired_contrasts.json').write_text(json.dumps(res,indent=2)); print(json.dumps(res,indent=2))
