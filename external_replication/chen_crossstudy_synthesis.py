import json, math
from pathlib import Path
from scipy.stats import norm, chi2

OUT=Path('analysis/chen_crossstudy_synthesis'); OUT.mkdir(parents=True,exist_ok=True)
studies=[
 {'study':'Chen et al. 2021 raw reanalysis','n':18,'rho':0.552,'status':'exploratory within-system right-behavior association'},
 {'study':'Chen et al. 2024 DANDI external test','n':8,'rho':2/3,'status':'prespecified secondary right-behavior transport association'}
]
z=[math.atanh(s['rho']) for s in studies];w=[s['n']-3 for s in studies]
zm=sum(a*b for a,b in zip(w,z))/sum(w);se=1/math.sqrt(sum(w));zstat=zm/se
Q=sum(wi*(zi-zm)**2 for wi,zi in zip(w,z));df=len(studies)-1;phet=1-chi2.cdf(Q,df);i2=max(0,(Q-df)/Q)*100 if Q>0 else 0
res={'analysis_status':'POST_HOC_DESCRIPTIVE_CROSS_STUDY_SYNTHESIS','endpoint':'right-ALM behavioral recovery versus control-trial bilateral ALM coupling','studies':studies,'fixed_effect_fisher_z':zm,'pooled_rho':math.tanh(zm),'approx_one_sided_p':1-norm.cdf(zstat),'approx_two_sided_p':2*(1-norm.cdf(abs(zstat))),'pooled_95_ci':[math.tanh(zm-1.96*se),math.tanh(zm+1.96*se)],'heterogeneity_Q':Q,'heterogeneity_p':phet,'I2_percent':i2,'interpretation':'Descriptive synthesis only. The 2024 right-side endpoint was prespecified but secondary, and its two-side Bonferroni-adjusted P value was 0.0832; therefore this pooled calculation is not confirmatory evidence and must not be used to relabel the frozen 2024 primary test as successful.'}
(OUT/'result.json').write_text(json.dumps(res,indent=2,sort_keys=True));print(json.dumps(res,indent=2,sort_keys=True))
