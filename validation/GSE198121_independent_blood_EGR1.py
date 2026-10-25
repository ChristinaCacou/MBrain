#!/usr/bin/env python3
"""GSE198121 independent living clinical peripheral-blood EGR1 check; not brain replication."""
import re,json,sys
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
df=pd.read_csv(sys.argv[1],compression="infer",low_memory=False)
out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
p={}
for col in df.columns:
 z=re.findall(r'(?<![A-Za-z0-9])([HN][1-7])(?![0-9])',str(col),re.I)
 if len(z)!=1:continue
 code=z[0].upper()
 if code in p:raise ValueError('Duplicate subject '+code)
 p[code]=col
expected={f'H{i}' for i in range(1,8)}|{f'N{i}' for i in range(1,8)}
if set(p)!=expected:raise ValueError('Bad sample headers; found '+str(p))
cc=[p[f'H{i}'] for i in range(1,8)];nc=[p[f'N{i}'] for i in range(1,8)]
val=df[cc+nc].apply(pd.to_numeric,errors='raise').to_numpy(dtype=float)
assert val.shape[1]==14 and np.isfinite(val).all()
gfields=[x for x in df.columns if x not in set(cc+nc)]
m=np.zeros(len(df),dtype=bool)
for col in gfields:
 s=df[col].astype(str)
 m |= s.str.strip().str.upper().isin(['EGR1','ENSG00000120738']).to_numpy()
 m |= s.str.contains(r'(?<![A-Za-z0-9])ENSG00000120738(?![0-9])',regex=True,case=False,na=False).to_numpy()
 m |= s.str.contains(r'(?<![A-Za-z0-9])EGR1(?![A-Za-z0-9])',regex=True,case=False,na=False).to_numpy()
case=val[:,:7];control=val[:,7:]
t,pval=ttest_ind(case,control,axis=1,equal_var=False,nan_policy='raise')
valid=np.isfinite(pval);q=np.full(len(df),np.nan)
if valid.any():q[valid]=multipletests(pval[valid],method='fdr_bh')[1]
delta=case.mean(axis=1)-control.mean(axis=1)
res=df[gfields].copy();res['mean_OUD_minus_control_VST']=delta;res['Welch_t']=t;res['Welch_p']=pval;res['BH_q']=q
res.to_csv(out/'GSE198121_blood_all_genes_Welch.csv',index=False)
hit=res.loc[m].copy();hit.to_csv(out/'GSE198121_EGR1_results.csv',index=False)
summary={'GEO':'GSE198121','n_OUD':7,'n_control':7,'sex':'male','tissue':'peripheral blood','input':'published VST normalized expression, NOT original counts','test':'exploratory Welch t, BH across finite rows','genes_with_finite_p':int(valid.sum()),'EGR1_matching_rows':int(m.sum()),'EGR1_all_negative':bool((delta[m]<0).all()) if m.any() else None,'EGR1_values':hit.to_dict('records'),'limitation':'Independent living clinical patient cohort, different tissue; not brain or transcriptome-wide model-matched replication, and cannot establish causality.'}
(out/'SUMMARY.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
print(json.dumps(summary,indent=2,default=str))
