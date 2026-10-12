import json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr

SEED=20261012; NDRAW=100000; NPERM=500000; BATCH=10000
ROOT=Path('analysis/zatka_primary_reconstruction'); OUT=Path('analysis/zatka_mapping_uncertainty');OUT.mkdir(parents=True,exist_ok=True)
TARGET=['ACA','MOp','MOs','ORB','PL','RSP','SSp','SSs','VISam','VISl','VISp','VISpm']
SECOND={'PL':.3333,'ORB':.3085,'VISpm':.2248,'MOp':.1869,'SSp':.1862,'VISl':.1722,'SSs':.1600,'MOs':.1529,'ACA':.1451,'RSP':.1421,'VISp':.1414,'VISam':.0818}

def canon(a):
    a=str(a)
    if a in ['MOp','MOs','SSs','VISam','VISl','VISp','VISpm','PL']:return a
    for root in ['ACA','ORB','RSP','SSp']:
        if a.startswith(root):return root
    return None

def perm_positive(x,y,seed):
    x=np.asarray(x,float);y=np.asarray(y,float);obs=float(spearmanr(x,y).statistic);rx=rankdata(x);ry=rankdata(y);xc=rx-rx.mean();denx=np.sqrt(xc@xc);rng=np.random.default_rng(seed);ge=done=0
    while done<NPERM:
        k=min(BATCH,NPERM-done);order=np.argsort(rng.random((k,len(y))),axis=1);P=ry[order];P-=P.mean(axis=1,keepdims=True);r=(P@xc)/(np.sqrt(np.sum(P*P,axis=1))*denx);ge+=int(np.sum(r>=obs-1e-15));done+=k
    return obs,float((ge+1)/(NPERM+1))

def aggregate(coord,cands,rule):
    sm={r:0. for r in TARGET};wt={r:0. for r in TARGET}
    for e,arr in zip(coord.total_effect,cands):
        if not arr:continue
        if rule=='first':alloc=[(arr[0],1.)]
        elif rule=='last':alloc=[(arr[-1],1.)]
        elif rule=='fractional':alloc=[(r,1/len(arr)) for r in arr]
        else:raise ValueError(rule)
        for r,w in alloc:sm[r]+=w*float(e);wt[r]+=w
    return np.array([sm[r]/wt[r] if wt[r]>0 else np.nan for r in TARGET]),wt

def rowcorr(RV,rfix):
    A=RV-RV.mean(axis=1,keepdims=True);b=rfix-rfix.mean();return (A@b)/(np.sqrt(np.sum(A*A,axis=1))*np.sqrt(np.sum(b*b)))

def main():
    coord=pd.read_csv(ROOT/'coordinate_effects_atlas_mapped.csv');feat=pd.read_csv('analysis/pui_blind_fixedC/blind_features_region.csv').set_index('region');pui=np.array([feat.loc[r,'pui_mean'] for r in TARGET],float);sec=np.array([SECOND[r] for r in TARGET],float)
    cands=[]
    for s in coord.dv_trace_unique.fillna(''):
        q=[]
        for a in str(s).split('|'):
            r=canon(a)
            if r and r not in q:q.append(r)
        cands.append(q)
    candidate_table=pd.DataFrame({'laserIdx':coord.laserIdx,'ML_mm':coord.ML_mm,'AP_mm':coord.AP_mm,'total_effect':coord.total_effect,'candidate_targets':['|'.join(x) for x in cands],'n_candidate_targets':[len(x) for x in cands]});candidate_table.to_csv(OUT/'coordinate_candidate_targets.csv',index=False)
    deterministic={}
    for rule in ['first','last','fractional']:
        v,w=aggregate(coord,cands,rule);ok=np.isfinite(v);rho,p=perm_positive(pui[ok],v[ok],SEED+len(rule));rho2=float(spearmanr(sec[ok],v[ok]).statistic);deterministic[rule]={'n_regions':int(ok.sum()),'effects':{r:(float(v[i]) if np.isfinite(v[i]) else None) for i,r in enumerate(TARGET)},'effective_coordinate_weights':{r:float(w[r]) for r in TARGET},'cpg_rho':rho,'cpg_positive_one_sided_permutation_p':p,'raw_vs_pinned_secondary_rho':rho2}
    rng=np.random.default_rng(SEED);nr=len(TARGET);sums=np.zeros((NDRAW,nr));counts=np.zeros((NDRAW,nr),dtype=np.int16);ri={r:i for i,r in enumerate(TARGET)}
    for e,arr in zip(coord.total_effect,cands):
        if not arr:continue
        z=rng.integers(0,len(arr),size=NDRAW)
        for j,r in enumerate(arr):
            m=z==j;sums[m,ri[r]]+=float(e);counts[m,ri[r]]+=1
    mapped=np.divide(sums,counts,where=counts>0,out=np.full_like(sums,np.nan));present=np.isfinite(mapped);complete=np.all(present,axis=1);V=mapped[complete]
    rp=rs=np.array([])
    if len(V):
        RV=rankdata(V,axis=1);rp=rowcorr(RV,rankdata(pui));rs=rowcorr(RV,rankdata(sec))
    topidx=np.argsort(np.where(present,mapped,-np.inf),axis=1)[:,-3:];topcounts=np.array([np.sum(topidx==j) for j in range(nr)]);pres=present.sum(axis=0)
    ensemble_regions=[]
    for j,r in enumerate(TARGET):
        vv=mapped[present[:,j],j]
        ensemble_regions.append({'region':r,'candidate_coordinate_count':int(sum(r in q for q in cands)),'mapping_presence_fraction':float(pres[j]/NDRAW),'effect_median':float(np.median(vv)) if len(vv) else None,'effect_2_5':float(np.percentile(vv,2.5)) if len(vv) else None,'effect_97_5':float(np.percentile(vv,97.5)) if len(vv) else None,'top3_fraction_given_present':float(topcounts[j]/pres[j]) if pres[j] else None})
    pd.DataFrame(ensemble_regions).to_csv(OUT/'region_mapping_ensemble.csv',index=False)
    result={'analysis_status':'POST_REVEAL_ATLAS_MAPPING_SENSITIVITY','n_coordinates':int(len(coord)),'n_random_mapping_draws':NDRAW,'candidate_rule':'Each coordinate may map to any frozen target whose Beryl acronym or collapsed parent appears anywhere along the exact ML/AP vertical atlas trace. Random draws choose uniformly among those admissible labels. This is a stress-test ensemble, not a posterior anatomical probability model.','n_complete_12_region_draws':int(complete.sum()),'complete_draw_fraction':float(complete.mean()),'complete_draw_cpg_rho_median':float(np.median(rp)) if len(rp) else None,'complete_draw_cpg_rho_2_5':float(np.percentile(rp,2.5)) if len(rp) else None,'complete_draw_cpg_rho_97_5':float(np.percentile(rp,97.5)) if len(rp) else None,'complete_draw_fraction_cpg_rho_positive':float(np.mean(rp>0)) if len(rp) else None,'complete_draw_raw_vs_pinned_secondary_rho_median':float(np.median(rs)) if len(rs) else None,'complete_draw_raw_vs_pinned_secondary_rho_2_5':float(np.percentile(rs,2.5)) if len(rs) else None,'complete_draw_raw_vs_pinned_secondary_rho_97_5':float(np.percentile(rs,97.5)) if len(rs) else None,'deterministic_rules':deterministic,'region_ensemble':ensemble_regions,'interpretation':'The raw coordinate-level outcome is highly sensitive to coordinate-to-region assignment, but the failed CPG prediction is not rescued: across complete 12-region admissible mapping draws the CPG rank association remains negative.'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True));print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
