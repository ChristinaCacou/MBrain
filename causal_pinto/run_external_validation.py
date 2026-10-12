import json, zipfile, hashlib, itertools
from pathlib import Path
import requests
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import rankdata

SEED=20261011
N_BOOT=50000
URLS={
 'fig2':'https://cdn.elifesciences.org/articles/70263/elife-70263-fig2-data1-v1.zip',
 'fig4':'https://cdn.elifesciences.org/articles/70263/elife-70263-fig4-data1-v1.zip',
}
REGIONS=['V1','mV2','PPC','RSC','mM2','aM2','M1']
ROOT=Path('analysis/pinto_external')
CACHE=ROOT/'validation_source'; CACHE.mkdir(parents=True,exist_ok=True)

def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def download_extract(key,url):
 zp=CACHE/f'{key}.zip'; r=requests.get(url,timeout=120); r.raise_for_status(); zp.write_bytes(r.content)
 td=CACHE/key; td.mkdir(exist_ok=True)
 with zipfile.ZipFile(zp) as z: z.extractall(td)
 mats=list(td.glob('*.mat'))
 if len(mats)!=1: raise RuntimeError(f'{key}: expected one .mat, got {mats}')
 return zp,mats[0]

def recursive_find(obj,key,path='root'):
 hits=[]
 if isinstance(obj,dict):
  for k,v in obj.items():
   p=f'{path}/{k}'
   if k==key: hits.append((p,v))
   hits.extend(recursive_find(v,key,p))
 elif isinstance(obj,(list,tuple)):
  for i,v in enumerate(obj): hits.extend(recursive_find(v,key,f'{path}[{i}]'))
 elif isinstance(obj,np.ndarray) and obj.dtype==object:
  for i,v in enumerate(obj.reshape(-1)): hits.extend(recursive_find(v,key,f'{path}[{i}]'))
 return hits

def get_unique(obj,key):
 hits=recursive_find(obj,key)
 if len(hits)!=1: raise RuntimeError(f'{key}: {len(hits)} hits {[h[0] for h in hits]}')
 return hits[0]

def clean_labels(x):
 a=np.asarray(x).reshape(-1)
 return [str(v).strip() for v in a]

def spearman(x,y):
 rx=rankdata(np.asarray(x,float)); ry=rankdata(np.asarray(y,float))
 rx=rx-rx.mean(); ry=ry-ry.mean(); den=np.sqrt(np.sum(rx*rx)*np.sum(ry*ry))
 return float(np.sum(rx*ry)/den) if den>0 else np.nan

def exact_perm_positive(x,y):
 x=np.asarray(x,float); y=np.asarray(y,float); obs=spearman(x,y); vals=[]
 for perm in itertools.permutations(y.tolist()): vals.append(spearman(x,np.asarray(perm,float)))
 vals=np.asarray(vals,float); p=float(np.mean(vals>=obs-1e-15))
 return obs,p,int(np.sum(vals>=obs-1e-15)),int(len(vals))

def bootstrap_ci(x,y):
 rng=np.random.default_rng(SEED); x=np.asarray(x,float); y=np.asarray(y,float); vals=[]
 for _ in range(N_BOOT):
  ii=rng.integers(0,len(x),len(x)); r=spearman(x[ii],y[ii])
  if np.isfinite(r): vals.append(r)
 return [float(np.percentile(vals,2.5)),float(np.percentile(vals,97.5))]

def main():
 lock=json.loads((ROOT/'EXTERNAL_TEST_LOCK.json').read_text())
 sources={}; mats={}
 for k,u in URLS.items():
  z,m=download_extract(k,u); mats[k]=loadmat(m,simplify_cells=True)
  sources[k]={'url':u,'zip_sha256':sha256(z),'mat_sha256':sha256(m),'mat_file':m.name}

 # perturbation fields
 _,fig2=get_unique(mats['fig2'],'logisticRegression_coefficients')
 labels=clean_labels(fig2['inactivation_condition_labels'])
 t=np.asarray(fig2['time_axis'],float).reshape(-1)
 beta=np.asarray(fig2['coefficient_mean'],float)
 fdr=np.asarray(fig2['coefficient_FDR_isSig'],float)
 if beta.shape[0]!=len(labels) and beta.shape[1]==len(labels): beta=beta.T; fdr=fdr.T
 if beta.shape!=(len(labels),len(t)): raise RuntimeError((beta.shape,len(labels),len(t)))
 nanidx=np.flatnonzero(~np.isfinite(t))
 if len(nanidx)<2: raise RuntimeError(f'expected two NaN separators in time axis, got {t}')
 pre_idx=np.arange(0,nanidx[0]); during_idx=np.arange(nanidx[0]+1,nanidx[1]); post_idx=np.arange(nanidx[1]+1,len(t))
 if len(pre_idx)<2 or len(during_idx)<1: raise RuntimeError((t,pre_idx,during_idx,post_idx))

 # recording fields
 path_area,areas_raw=get_unique(mats['fig4'],'expfit_areaList')
 path_tau,tau_raw=get_unique(mats['fig4'],'expFit_tau_median')
 path_tau_mouse,tau_mouse_raw=get_unique(mats['fig4'],'expFit_tau_bymouse')
 areas=clean_labels(areas_raw); tau=np.asarray(tau_raw,float).reshape(-1)
 if len(areas)!=len(tau): raise RuntimeError((areas,tau))
 tau_map=dict(zip(areas,tau))

 rows=[]
 pre_t=t[pre_idx]
 spacing=float(np.nanmedian(np.abs(np.diff(pre_t)))) if len(pre_t)>1 else 1.0
 for reg in REGIONS:
  if reg not in labels: raise RuntimeError(f'{reg} missing from perturb labels {labels}')
  if reg not in tau_map: raise RuntimeError(f'{reg} missing from tau labels {areas}')
  i=labels.index(reg); b=beta[i]; sig=fdr[i]
  w=np.maximum(0.0,-b[pre_idx])
  horizon=float(np.sum(w*np.abs(pre_t))/np.sum(w)) if np.sum(w)>0 else np.nan
  pre_auc=float(np.sum(w)*spacing)
  siglags=np.abs(pre_t[(sig[pre_idx]>0.5)&(b[pre_idx]<0)])
  oldest=float(np.max(siglags)) if len(siglags) else 0.0
  immediate=float(np.nanmean(np.maximum(0.0,-b[during_idx])))
  postdef=float(np.nanmean(np.maximum(0.0,-b[post_idx]))) if len(post_idx) else np.nan
  rec={'region':reg,'tau_median_s':float(tau_map[reg]),'causal_memory_horizon_s':horizon,'prelaser_deficit_auc':pre_auc,'oldest_FDR_significant_pre_laser_lag_s':oldest,'immediate_during_laser_deficit':immediate,'mean_postlaser_deficit':postdef}
  for j in range(len(t)):
   rec[f'beta_{j}']=float(b[j]) if np.isfinite(b[j]) else np.nan
   rec[f'time_{j}_s']=float(t[j]) if np.isfinite(t[j]) else np.nan
   rec[f'fdrsig_{j}']=float(sig[j]) if np.isfinite(sig[j]) else np.nan
  rows.append(rec)
 df=pd.DataFrame(rows)
 ok=np.isfinite(df.tau_median_s)&np.isfinite(df.causal_memory_horizon_s)
 x=df.loc[ok,'tau_median_s'].to_numpy(); y=df.loc[ok,'causal_memory_horizon_s'].to_numpy()
 rho,p,count,total=exact_perm_positive(x,y)
 ci=bootstrap_ci(x,y)
 loo=[]
 for i in range(len(x)):
  k=np.arange(len(x))!=i; loo.append(spearman(x[k],y[k]))
 secondary={}
 for key in ['prelaser_deficit_auc','oldest_FDR_significant_pre_laser_lag_s','immediate_during_laser_deficit','mean_postlaser_deficit']:
  yy=df.loc[ok,key].to_numpy(float); secondary[key]={'spearman_rho':spearman(x,yy)}
 result={
  'analysis_status':'external post-publication replication/stress test; specification locked before statistic calculation but source study and qualitative hypothesis were public',
  'dataset_doi':'10.7554/eLife.70263',
  'n_regions':int(np.sum(ok)),'regions':df.loc[ok,'region'].tolist(),
  'source_hashes':sources,
  'recording_field_paths':{'area':path_area,'tau':path_tau,'tau_bymouse':path_tau_mouse},
  'perturbation_labels':labels,'time_axis_s':[None if not np.isfinite(v) else float(v) for v in t],
  'segment_indices':{'pre':pre_idx.tolist(),'during':during_idx.tolist(),'post':post_idx.tolist()},
  'primary_predictor':'tau_median_s','primary_outcome':'causal_memory_horizon_s',
  'primary_spearman_rho':rho,'exact_one_sided_permutation_p':p,'permutations_ge_observed':count,'n_exact_permutations':total,
  'bootstrap_95_ci':ci,'leave_one_region_out_rho':loo,'leave_one_region_out_min':float(np.nanmin(loo)),'leave_one_region_out_median':float(np.nanmedian(loo)),'leave_one_region_out_max':float(np.nanmax(loo)),
  'success_criterion':'rho > 0 and exact one-sided permutation p < 0.05','success':bool(rho>0 and p<0.05),
  'secondary_correlations':secondary,
  'guardrail':'This result does not validate PUI. It tests alignment between an observational temporal-timescale feature and the temporal horizon of an independent perturbation effect in a different task/lab.'
 }
 df.to_csv(ROOT/'external_region_table.csv',index=False)
 (ROOT/'external_result.json').write_text(json.dumps(result,indent=2,sort_keys=True))
 (ROOT/'EXTERNAL_RESULT.md').write_text(
  '# Pinto 2022 external perturbation stress test\n\n'
  f"**{'PASS' if result['success'] else 'FAIL'}** under the locked external-test criterion.\n\n"
  f"Primary Spearman rho = {rho:.6f}; exact one-sided permutation p = {p:.6g} ({count}/{total} permutations >= observed).\n\n"
  f"Bootstrap 95% CI = [{ci[0]:.4f}, {ci[1]:.4f}]. Leave-one-region-out rho range = [{result['leave_one_region_out_min']:.4f}, {result['leave_one_region_out_max']:.4f}].\n\n"
  'This is a post-publication external replication/stress test, not a preregistered or blinded test, and it does not rescue or validate the failed PUI result.\n')
 print(df.to_string(index=False)); print(json.dumps(result,indent=2))

if __name__=='__main__': main()
