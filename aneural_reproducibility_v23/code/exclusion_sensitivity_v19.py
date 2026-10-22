#!/usr/bin/env python3
"""Independent exclusion-restriction stress test for the v18 Gaussian model.
Exploratory synthetic analysis, not biological or preregistered data.

Model: Q, eps_M, eps_E independent N(0,1); M=sqrt(rho)*Q+sqrt(1-rho)*eps_M;
E=sqrt(rho)*Q+sqrt(1-rho)*eps_E under no M->E link.
Treatment Z shifts Q by b Z in the no-transfer model, so ΔE=sqrt(rho)*b.
A true selective do(M) with ΔM=1 yields ΔE=gamma in the coupled model.
Equal values of ΔE produce observational indistinguishability of the 
intervention mean alone unless off-target b is measured/bounded.

Per day: 20 observations per randomized arm, 12 independent days. Within-day
contrast noise variance=2/20. Two-sided t-test over day differences (df=11).
Under the Gaussian model, exact power follows the noncentral t distribution.
Simulation samples the analytically equivalent per-day Gaussian contrasts and
validates power with independent RNG seed. Simulation replicates are synthetic
campaigns, not actual biological replication units.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import t,nct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RHO=.6;GAMMA=.25;N_DAYS=12;N_ARM=20;ALPHA=.05;N_REP=20000;SEED=2026102119
OUT=Path(__file__).resolve().parent.parent/'results'
OUT.mkdir(exist_ok=True)
C=np.sqrt(RHO)
BSTAR=GAMMA/C
qshift_grid=[0,.05,.10,.20,BSTAR,.40,.50]
crit=t.ppf(1-ALPHA/2,N_DAYS-1)
sigma_day=np.sqrt(2/N_ARM)
rng=np.random.default_rng(SEED)
rows=[]
for b in qshift_grid:
    delta=C*b
    ncp=delta*np.sqrt(N_DAYS)/sigma_day
    exact=float(nct.sf(crit, N_DAYS-1,ncp)+nct.cdf(-crit,N_DAYS-1,ncp))
    differences=delta+sigma_day*rng.standard_normal((N_REP,N_DAYS))
    statistics=np.sqrt(N_DAYS)*differences.mean(axis=1)/differences.std(axis=1,ddof=1)
    sim=float(np.mean(np.abs(statistics)>crit))
    rows.append(dict(offtarget_Q_shift_SD=round(b,9),population_E_shift=round(delta,9),
                     analytical_rejection_rate=round(exact,6),mc_rejection_rate=round(sim,6),
                     mc_se=round(float(np.sqrt(sim*(1-sim)/N_REP)),6)))
# selective directed intervention: E noise conditioned on target intervention
# has variance c^2*rho+v in original directed model.
c=(RHO-GAMMA)/RHO
v=1-GAMMA**2-c*c*RHO-2*GAMMA*c*RHO
sigma_coupled=np.sqrt((2/N_ARM)*(c*c*RHO+v))
ncp_coupled=GAMMA*np.sqrt(N_DAYS)/sigma_coupled
coupled_power=float(nct.sf(crit,N_DAYS-1,ncp_coupled)+nct.cdf(-crit,N_DAYS-1,ncp_coupled))
rows.append(dict(offtarget_Q_shift_SD=np.nan,population_E_shift=GAMMA,
                 analytical_rejection_rate=round(coupled_power,6),mc_rejection_rate=np.nan,mc_se=np.nan))
df=pd.DataFrame(rows)
df['scenario']=['no-transfer off-target Q shift']*len(qshift_grid)+['selective directed M shift']
df.to_csv(OUT/'exclusion_restriction_sensitivity_v19.csv',index=False)

x=np.linspace(0,.5,251)
y=nct.sf(crit,N_DAYS-1,C*x*np.sqrt(N_DAYS)/sigma_day)+nct.cdf(-crit,N_DAYS-1,C*x*np.sqrt(N_DAYS)/sigma_day)
fig,ax=plt.subplots(figsize=(7.6,4.7))
ax.plot(x,y,color='black',linewidth=1.75,label='No M→E path; off-target Q shift')
ax.scatter([r['offtarget_Q_shift_SD'] for r in rows[:-1]], [r['mc_rejection_rate'] for r in rows[:-1]],
           color='black',s=22,marker='o',label='Independent Monte Carlo (20,000/case)')
ax.axvline(BSTAR,color='black',linestyle='--',alpha=.6)
ax.axhline(.05,color='black',linestyle=':',alpha=.6)
ax.annotate('Q shift that mimics a direct effect of 0.25\n(b = 0.25 / √0.60 ≈ 0.323)',
            xy=(BSTAR,float(np.interp(BSTAR,x,y))),xytext=(.14,.88),
            textcoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.7),fontsize=9)
ax.set(xlabel='Unintended treatment shift in latent Q (SD units)',ylabel='Probability of significant day-level test',
       title='Synthetic intervention specificity collapses with shared-driver leakage',ylim=(0,1.02),xlim=(0,.5))
ax.legend(loc='lower right',fontsize=8,frameon=False)
fig.tight_layout()
fig.savefig(OUT/'Figure_S3_exclusion_sensitivity.png',dpi=210)
plt.close(fig)

report={'rho':RHO,'gamma':GAMMA,'target_M_shift':1,'Q_shift_matching_target_effect':BSTAR,
'null_E_effect_per_unit_Q_shift':C,'n_days':N_DAYS,'per_arm_per_day':N_ARM,'alpha':ALPHA,
'monte_carlo_campaigns_per_condition':N_REP,'mc_seed':SEED,'df':N_DAYS-1,
'coupled_E_residual_variance':c*c*RHO+v,'coupled_analytical_day_t_power':coupled_power,
'null_analytical_day_t_power_at_matching_leakage':rows[4]['analytical_rejection_rate'],
'new_experiment':'none','data_type':'synthetic parameter sensitivity',
'validation':'analytic noncentral t computation compared against independently sampled day contrasts',
'limitations':'Assumes linear-Gaussian model, fixed exogenous effects, known coefficients, independent day contrasts and balanced randomized arms; not a biological calibration or general bound.'}
(OUT/'exclusion_sensitivity_run_manifest_v19.json').write_text(json.dumps(report,indent=2))
print(df.to_string(index=False))
print(json.dumps(report,indent=2))
