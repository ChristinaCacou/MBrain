import io,bisect,json,requests,zipfile
from pathlib import Path
import numpy as np,pandas as pd,h5py
from scipy.io import loadmat
from scipy.stats import spearmanr
META='https://api.figshare.com/v2/articles/13008038'; OUT=Path('analysis/zatka_primary_reconstruction'); OUT.mkdir(parents=True,exist_ok=True)
TARGET=['ACA','MOp','MOs','ORB','PL','RSP','SSp','SSs','VISam','VISl','VISp','VISpm']
SECOND={'PL':.3333,'ORB':.3085,'VISpm':.2248,'MOp':.1869,'SSp':.1862,'VISl':.1722,'SSs':.1600,'MOs':.1529,'ACA':.1451,'RSP':.1421,'VISp':.1414,'VISam':.0818}
# Range-readable concatenation of the four Figshare split-ZIP parts.
meta=requests.get(META,timeout=60).json(); parts=sorted([f for f in meta['files'] if f['name'].startswith('Zatka-Haas_et_al_Dataset.zip.')],key=lambda f:f['name']); urls=[x['download_url'] for x in parts]; sizes=[int(x['size']) for x in parts]; starts=[]; total=0
for z in sizes: starts.append(total); total+=z
class RC(io.RawIOBase):
 def __init__(self): self.pos=0; self.s=requests.Session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,o,w=0): self.pos=max(0,min(total,o if w==0 else self.pos+o if w==1 else total+o)); return self.pos
 def read(self,n=-1):
  if self.pos>=total:return b''
  if n<0:n=total-self.pos
  n=min(n,total-self.pos); out=[]
  while n:
   i=bisect.bisect_right(starts,self.pos)-1; local=self.pos-starts[i]; take=min(n,sizes[i]-local); r=self.s.get(urls[i],headers={'Range':f'bytes={local}-{local+take-1}'},timeout=120); r.raise_for_status(); b=r.content
   if len(b)!=take: raise RuntimeError((i,take,len(b),r.status_code))
   out.append(b); self.pos+=take; n-=take
  return b''.join(out)
with zipfile.ZipFile(RC()) as z:
 for member in ['OptogeneticInactivation/26CoordSet.mat','OptogeneticInactivation/Inactivation_52Coord.mat']:
  p=OUT/Path(member).name
  with z.open(member) as src,open(p,'wb') as dst:
   for b in iter(lambda:src.read(1<<20),b''): dst.write(b)
# Decode MATLAB table via its MCOS variable-name/data cell arrays.
with h5py.File(OUT/'Inactivation_52Coord.mat','r') as f:
 names=[''.join(chr(int(x)) for x in np.asarray(f[r]).reshape(-1)) for r in f['#refs#/1'][()].reshape(-1)]; refs=f['#refs#/y'][()].reshape(-1); need={'stimulus','response','laserIdx','sessionID'}; d={n:np.asarray(f[r]).T for n,r in zip(names,refs) if n in need}
stim=d['stimulus']; response=d['response'].ravel().astype(int); laser=d['laserIdx'].ravel().astype(int); session=d['sessionID'].ravel().astype(int); keep=(stim[:,0]==stim[:,1])&(stim[:,0]>0); response=response[keep]; laser=laser[keep]; session=session[keep]
sessions=np.unique(session); delta=np.full((3,52,len(sessions)),np.nan); ntr=np.zeros((52,len(sessions)),int)
for si,s in enumerate(sessions):
 m=session==s; rr=response[m]; ll=laser[m]; bc=np.array([(rr[ll==0]==r).sum() for r in (1,2,3)],float)
 if bc.sum()==0: continue
 base=bc/bc.sum()
 for j in range(1,53):
  cc=np.array([(rr[ll==j]==r).sum() for r in (1,2,3)],float); ntr[j-1,si]=int(cc.sum())
  if cc.sum(): delta[:,j-1,si]=cc/cc.sum()-base
delta=np.nanmean(delta,axis=2); effect=np.sqrt(np.sum(delta**2,axis=0)); cs=loadmat(OUT/'26CoordSet.mat')['coordSet']; coords=np.vstack([cs,np.column_stack([-cs[:,0],cs[:,1]])]); coord=pd.DataFrame({'laserIdx':np.arange(1,53),'ML_mm':coords[:,0],'AP_mm':coords[:,1],'delta_left':delta[0],'delta_right':delta[1],'delta_nogo':delta[2],'total_effect':effect,'n_sessions':np.sum(ntr>0,axis=1)})
# Anatomical mapping fixed independently of outcome values: exact laser center -> first target Beryl cortical label when descending from +0.5 to -2.0 mm DV.
from iblatlas.atlas import AllenAtlas
ba=AllenAtlas(res_um=50); br=ba.regions; zs=np.arange(.0005,-.00205,-.00005); mapped=[]; traces=[]
for ml,ap in coords:
 xyz=np.column_stack([np.full(len(zs),ml/1000),np.full(len(zs),ap/1000),zs]); ids=ba.get_labels(xyz,mapping='Beryl'); ac=np.asarray(br.id2acronym(ids,mapping='Beryl')).astype(str); hit=next((a for a in ac if a in TARGET),None); mapped.append(hit); traces.append('|'.join(dict.fromkeys(ac.tolist())))
coord['region']=mapped; coord['dv_trace_unique']=traces; coord.to_csv(OUT/'coordinate_effects_atlas_mapped.csv',index=False)
reg=coord.dropna(subset=['region']).groupby('region',as_index=False).agg(primary_raw_effect=('total_effect','mean'),n_coordinates=('laserIdx','size'),mean_sessions_per_coordinate=('n_sessions','mean'))
# Frozen feature table from the inherited outcome-blind branch.
f=pd.read_csv('analysis/pui_blind_fixedC/blind_features_region.csv')[['region','pui_mean','marginal_mean','local_auc_mean','local_bacc_mean']]; tab=f.merge(reg,on='region'); tab['secondary_effect']=tab.region.map(SECOND); tab.to_csv(OUT/'primary_region_table.csv',index=False)
def perm_p(x,y,n=500000,seed=20261012):
 obs=float(spearmanr(x,y).statistic); rng=np.random.default_rng(seed); ge=0
 for _ in range(n//10000):
  P=np.array([rng.permutation(y) for __ in range(10000)]); rx=np.argsort(np.argsort(x)).astype(float); rx-=rx.mean(); ry=np.argsort(np.argsort(P,axis=1),axis=1).astype(float); ry-=ry.mean(axis=1,keepdims=True); q=(ry@rx)/np.sqrt(np.sum(ry*ry,axis=1)*np.sum(rx*rx)); ge+=int(np.sum(q>=obs-1e-15))
 return obs,(ge+1)/(n+1)
res={'status':'POST_REVEAL_PRIMARY_ZATKA_RECONSTRUCTION','mapping_rule':'First frozen-target Beryl label encountered at exact laser ML/AP center while descending DV +0.5 to -2.0 mm, AllenAtlas 50um. Mapping uses anatomy only, not outcome magnitude.','n_raw_trials_total':int(len(stim)),'n_equal_nonzero_contrast_trials':int(len(response)),'n_sessions':int(len(sessions)),'n_coordinates':52,'n_coordinates_mapped_to_frozen_targets':int(coord.region.notna().sum()),'mapped_regions':reg.region.tolist(),'n_frozen_regions_recovered':int(len(tab))}
if len(tab)>=3:
 res['cpg_vs_primary_raw_effect']={}; res['cpg_vs_primary_raw_effect']['rho'],res['cpg_vs_primary_raw_effect']['one_sided_permutation_p']=perm_p(tab.pui_mean.to_numpy(),tab.primary_raw_effect.to_numpy())
 res['primary_vs_secondary_effect_spearman']=float(spearmanr(tab.primary_raw_effect,tab.secondary_effect).statistic); res['primary_vs_secondary_effect_mae']=float(np.mean(np.abs(tab.primary_raw_effect-tab.secondary_effect))); res['primary_vs_secondary_effect_max_abs_error']=float(np.max(np.abs(tab.primary_raw_effect-tab.secondary_effect)))
(OUT/'result.json').write_text(json.dumps(res,indent=2)); print(json.dumps(res,indent=2)); print(coord[['laserIdx','ML_mm','AP_mm','total_effect','region']].to_string(index=False)); print(tab.to_string(index=False))
