#!/usr/bin/env python3
"""Constructive joint interventional observational equivalence.
Synthetic Gaussian construction from v18 rho=0.6 gamma=0.25.
No real measurements. Both structural models yield identical p(M,E|Z)
although only one contains the causal edge M->E. Exclusion failure is explicit.
"""
from pathlib import Path
import numpy as np, json, pandas as pd
RHO=.6;GAMMA=.25;SEED=2026102120;N=150000
S=np.sqrt(RHO); C=(RHO-GAMMA)/RHO
V=1-(GAMMA**2+C*C*RHO+2*GAMMA*C*RHO)
B=GAMMA/S
OUT=Path(__file__).resolve().parent.parent/'results'
OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(SEED)
rows=[]
for model in ['directed_selective','null_nonselective']:
  for z in [0,1]:
    q=rng.normal(size=N); em=rng.normal(size=N); ee=rng.normal(size=N)
    m0=S*q+np.sqrt(1-RHO)*em
    if model=='directed_selective':
      m=m0+z
      e=GAMMA*m+C*S*q+np.sqrt(V)*ee
    else:
      qp=q+B*z
      m=S*qp+np.sqrt(1-RHO)*em+(1-GAMMA)*z
      e=S*qp+np.sqrt(1-RHO)*ee
    cov=np.cov(m,e,ddof=1)
    rows.append(dict(model=model,intervention_z=z,n=N,mean_M=m.mean(),mean_E=e.mean(),
                     variance_M=cov[0,0],variance_E=cov[1,1],covariance_ME=cov[0,1]))
pd.DataFrame(rows).to_csv(OUT/'joint_intervention_distribution_check_v19.csv',index=False)
assert abs(S*B-GAMMA)<1e-12
proof={'rho':RHO,'gamma':GAMMA,'q_shift_b':B,'null_direct_Z_to_M_effect':1-GAMMA,
'null_total_first_stage_M_shift':S*B+1-GAMMA,'null_downstream_E_shift':S*B,
'directed_first_stage_M_shift':1,'directed_downstream_E_shift':GAMMA,
'within_arm_covariance_both_models':[[1,RHO],[RHO,1]],
'group_Z_0_means_both_models':[0,0],'group_Z_1_means_both_models':[1,GAMMA],
'exact_distributional_equivalence':True,
'claim_scope':'Bivariate Gaussian observables M,E and randomized treatment Z; latent Q unobserved; specified intervention versions differ. Not equivalence of all counterfactual/interventional queries.',
'numerical_check_seed':SEED,'numerical_check_n_per_model_arm':N,
'data_type':'synthetic computational audit, no biological samples'}
(OUT/'joint_intervention_equivalence_proof_v19.json').write_text(json.dumps(proof,indent=2))
print(pd.DataFrame(rows).to_string(index=False,float_format=lambda f:f'{f:.6f}'))
print(json.dumps(proof,indent=2))
