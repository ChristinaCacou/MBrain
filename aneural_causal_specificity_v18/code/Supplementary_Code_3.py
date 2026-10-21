"""Measurement SNR adversarial control; this is a synthetic measurement model, not biology."""
from pathlib import Path
import numpy as np,pandas as pd
from numpy.random import default_rng
from scipy.stats import wilcoxon
from scipy.linalg import svdvals
OUT=Path(__file__).resolve().parent
rows=[]
for seed in range(200):
 r=default_rng(seed+89000);T=120;D=40
 t=np.arange(T)
 u=np.column_stack((np.exp(-t/25)*np.sin(t/8), np.exp(-t/30)*np.cos(t/12)))
 load=r.normal(size=(2,D)); signal=u@load
 noise=r.normal(size=(T,D));noise-=noise.mean(axis=0)
 for regime,g in [('BASE',1.),('GAIN045',.45)]:
  for nm,observed in [('NOISE_FREE',g*signal),('FIXED_NOISE',g*signal+0.2*noise),('MATCHED_SNR',g*(signal+0.2*noise))]:
   a=observed-observed.mean(axis=0);a=a/(np.linalg.norm(a)+1e-15)
   s=svdvals(a);v=s*s
   pr=(v.sum()**2)/(np.sum(v*v)+1e-15)
   # quantization-derived bitstream proxy is sensitive to noise; report only PR as primary
   rows.append((seed,regime,nm,float(pr)))
d=pd.DataFrame(rows,columns=['network_seed','regime','measurement_condition','effective_dimension'])
d.to_csv(OUT/'ExtendedData_Table5_measurement_snr.csv',index=False)
out=[]
for condition,block in d.groupby('measurement_condition'):
 piv=block.pivot(index='network_seed',columns='regime',values='effective_dimension');diff=piv['GAIN045']-piv['BASE'];
 out.append((condition,piv['BASE'].mean(),piv['GAIN045'].mean(),diff.mean(),float(wilcoxon(diff).pvalue),int((diff>0).sum()),len(diff)))
s=pd.DataFrame(out,columns=['measurement_condition','base_mean_dimension','gain045_mean_dimension','paired_difference','wilcoxon_p','number_positive','n'])
s.to_csv(OUT/'ExtendedData_Table5_summary.csv',index=False)
print(s.to_string(index=False))