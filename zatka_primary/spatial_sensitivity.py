import json,urllib.request
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr,t as tdist
from scipy.spatial.distance import cdist
SEED=20261012; NPERM=500000; BATCH=10000
URLS=[f'https://raw.githubusercontent.com/tipsyfermion/NMA_Mapping_Brain_Networks_2020/master/code/steinmetz_part{i}.npz' for i in range(3)]
OUT=Path('analysis/zatka_spatial'); CACHE=Path('.cache/steinmetz'); OUT.mkdir(parents=True,exist_ok=True); CACHE.mkdir(parents=True,exist_ok=True)
def alias(a):
 a=str(a).strip()
 if a.startswith('SSp'): return 'SSp'
 if a.startswith('ACA'): return 'ACA'
 if a.startswith('ORB'): return 'ORB'
 return a
def residual(v,Z):
 A=np.column_stack([np.ones(len(v)),Z]); return v-A@np.linalg.lstsq(A,v,rcond=None)[0]
def partial_rank_perm(x,y,Z):
 rx=rankdata(x); ry=rankdata(y); xr=residual(rx,Z); yr=residual(ry,Z); obs=float(np.corrcoef(xr,yr)[0,1]); rng=np.random.default_rng(SEED); ge=tw=done=0
 # fixed x residual; each permuted y rank is residualized on the same spatial basis
 while done<NPERM:
  k=min(BATCH,NPERM-done); keys=rng.random((k,len(ry))); ord=np.argsort(keys,axis=1); P=ry[ord]; A=np.column_stack([np.ones(len(ry)),Z]); H=A@np.linalg.pinv(A); R=np.eye(len(ry))-H; PY=P@R.T; num=PY@xr; den=np.sqrt(np.sum(PY*PY,axis=1)*np.sum(xr*xr)); q=num/den; ge+=int(np.sum(q>=obs-1e-15)); tw+=int(np.sum(np.abs(q)>=abs(obs)-1e-15)); done+=k
 return {'rho':obs,'permutation_p_one_sided':(ge+1)/(NPERM+1),'permutation_p_two_sided':(tw+1)/(NPERM+1)}
def moran(v,W):
 z=v-np.mean(v); return float(len(v)/W.sum()*(z@W@z)/(z@z))
def gls_rank(x,y,D,length):
 X=np.column_stack([np.ones(len(x)),rankdata(x)]); yy=rankdata(y); S=np.exp(-D/length)+np.eye(len(x))*1e-6; Si=np.linalg.inv(S); XtSi=X.T@Si; beta=np.linalg.solve(XtSi@X,XtSi@yy); r=yy-X@beta; df=len(x)-2; s2=float(r@Si@r/df); cov=s2*np.linalg.inv(XtSi@X); se=float(np.sqrt(cov[1,1])); tt=float(beta[1]/se); p=float(1-tdist.cdf(tt,df)); return {'length_um':float(length),'slope_rank':float(beta[1]),'se':se,'t':tt,'df':df,'one_sided_p':p}
# load revealed outcome/features
m=pd.read_csv('analysis/pui_blind_fixedC/reveal_vectorized/revealed_region_table.csv'); regs=set(m.region)
coords={r:[] for r in regs}
for i,u in enumerate(URLS):
 p=CACHE/f'part{i}.npz'
 if not p.exists(): urllib.request.urlretrieve(u,p)
 with np.load(p,allow_pickle=True) as z: ss=z['dat'].tolist()
 for d in ss:
  aa=np.asarray([alias(x) for x in np.asarray(d['brain_area']).astype(str)]); c=np.asarray(d['ccf'],float)
  for r in regs:
   q=c[aa==r]
   if len(q): coords[r].append(q.mean(axis=0))
cent={r:np.mean(coords[r],axis=0) for r in regs if coords[r]}
C=np.vstack([cent[r] for r in m.region]); m[['ccf_0','ccf_1','ccf_2']]=C
x=m.pui_mean.to_numpy(float); y=m.silencing_effect.to_numpy(float); Z=(C-C.mean(axis=0))/C.std(axis=0,ddof=1)
D=cdist(C,C); nn=np.median(np.partition(D+np.eye(len(D))*1e99,1,axis=1)[:,1]); med=np.median(D[np.triu_indices(len(D),1)]); W=np.exp(-D/max(nn,1)); np.fill_diagonal(W,0)
# centered Moran operator; retain two strongest positive spatial modes
M=np.eye(len(W))-np.ones_like(W)/len(W); vals,vecs=np.linalg.eigh(M@W@M); order=np.argsort(vals)[::-1]; pos=[i for i in order if vals[i]>1e-10]; E=vecs[:,pos[:min(2,len(pos))]]
res={'n_regions':len(m),'regions':m.region.tolist(),'centroid_coordinate_units':'source Steinmetz CCF units (micrometres)','median_nearest_neighbor_um':float(nn),'median_pairwise_um':float(med),'raw_spearman_rho':float(spearmanr(x,y).statistic),'moran_I_pui':moran(rankdata(x),W),'moran_I_outcome':moran(rankdata(y),W),'coordinate_adjusted_partial_rank':partial_rank_perm(x,y,Z),'moran_eigenvector_adjusted_partial_rank':partial_rank_perm(x,y,E) if E.shape[1] else None,'n_positive_moran_modes':len(pos),'n_moran_modes_adjusted':int(E.shape[1]),'gls_rank_sensitivities':[gls_rank(x,y,D,l) for l in [nn,med,2*med]],'status':'POST_REVEAL_SPATIAL_SENSITIVITY','interpretation_note':'These analyses are descriptive post-reveal spatial robustness checks, not replacements for the locked primary region-label permutation test.'}
m.to_csv(OUT/'region_table_with_ccf_centroids.csv',index=False); (OUT/'spatial_sensitivity.json').write_text(json.dumps(res,indent=2)); print(json.dumps(res,indent=2)); print(m[['region','pui_mean','silencing_effect','ccf_0','ccf_1','ccf_2']].to_string(index=False))
