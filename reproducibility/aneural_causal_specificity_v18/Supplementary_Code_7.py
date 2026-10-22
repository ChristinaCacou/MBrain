#!/usr/bin/env python3
"""Observational-equivalence / interventional-discrimination audit.

Purely synthetic causal counterexample and design assurance. Not Stentor data,
not preregistered, and not an amendment to any registered confirmatory analysis.
Run: python Supplementary_Code_7.py --outdir /path/to/output
Python >=3.10; numpy, scipy, pandas, matplotlib.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
import matplotlib.pyplot as plt

SEED = 2026101707
RHO = .6   # identical observational covariance in both models
GAMMA = .25  # directed effect in the integrated model, absent in common cause
C = (RHO-GAMMA)/RHO
V = 1-(GAMMA**2 + C**2*RHO+2*GAMMA*C*RHO)
assert V>0

def natural_samples(n,rng,model):
    q=rng.standard_normal(n)
    e_m=rng.standard_normal(n)
    e_e=rng.standard_normal(n)
    m=np.sqrt(RHO)*q+np.sqrt(1-RHO)*e_m
    if model=='common_cause':
        e=np.sqrt(RHO)*q+np.sqrt(1-RHO)*e_e
    else:
        e=GAMMA*m+C*np.sqrt(RHO)*q+np.sqrt(V)*e_e
    return m,e

def campaign(rng, scenario, n_days=12, per_arm=20):
    # Within each recording day, intervention Z is balanced 20/20.
    # The day effect is shared and cancels in within-day contrasts.
    contrasts=[]
    for _ in range(n_days):
        z=np.repeat([0,1],per_arm)
        rng.shuffle(z)
        q=rng.standard_normal(len(z))
        ee=rng.standard_normal(len(z))
        day_effect=.30*rng.standard_normal()
        if scenario=='coupled':
            e=GAMMA*z+C*np.sqrt(RHO)*q+np.sqrt(V)*ee+day_effect
        elif scenario=='common_cause':
            e=np.sqrt(RHO)*q+np.sqrt(1-RHO)*ee+day_effect
        elif scenario=='offtarget_common_cause':
            # Off-target intervention shifts q: the exclusion restriction fails.
            # This is a common cause only, with no M->E causal arrow.
            q=q+GAMMA*z/np.sqrt(RHO)
            e=np.sqrt(RHO)*q+np.sqrt(1-RHO)*ee+day_effect
        else:raise ValueError(scenario)
        contrasts.append(e[z==1].mean()-e[z==0].mean())
    d=np.asarray(contrasts)
    t_stat=d.mean()/(d.std(ddof=1)/np.sqrt(len(d)))
    p=2*student_t.sf(abs(t_stat),len(d)-1)
    critical=student_t.ppf(.975,len(d)-1)
    interval=(d.mean()-critical*d.std(ddof=1)/np.sqrt(len(d)),d.mean()+critical*d.std(ddof=1)/np.sqrt(len(d)))
    return d.mean(),p,interval

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,required=True)
    out=ap.parse_args().outdir;out.mkdir(exist_ok=True,parents=True)
    rng=np.random.default_rng(SEED)
    # Both models are exactly bivariate normal with mean zero, var 1 and covariance rho.
    rows=[]
    for model in ('common_cause','coupled'):
        m,e=natural_samples(50000,rng,model)
        rows.append(dict(model=model,n=50000,mean_M=float(m.mean()),mean_E=float(e.mean()),
                         var_M=float(m.var(ddof=1)),var_E=float(e.var(ddof=1)),
                         covariance=float(np.cov(m,e)[0,1]),correlation=float(np.corrcoef(m,e)[0,1])))
    pd.DataFrame(rows).to_csv(out/'ExtendedData_Table9_observational_equivalence.csv',index=False)
    scenarios=['common_cause','coupled','offtarget_common_cause']
    N=4000
    summ=[]; long=[]
    for sc in scenarios:
        effects=np.empty(N);pvals=np.empty(N);lows=np.empty(N);highs=np.empty(N)
        for j in range(N):
            effects[j],pvals[j],(lows[j],highs[j])=campaign(rng,sc)
        # One planned illustration per case is indexed 0; the probabilities summarize the full ensemble.
        summ.append(dict(scenario=sc,n_campaigns=N,days=12,cells_per_arm_day=20,
                         expected_intervention_effect=0 if sc=='common_cause' else GAMMA,
                         mean_estimate=effects.mean(),sd_estimate=effects.std(ddof=1),
                         rejection_rate=(pvals<.05).mean(),
                         avg_ci_coverage=np.mean((lows <= (0 if sc=='common_cause' else GAMMA)) & ((0 if sc=='common_cause' else GAMMA)<=highs)),
                         first_campaign_effect=effects[0],first_campaign_p=pvals[0]))
        long += [dict(scenario=sc,campaign_id=j+1,effect=effects[j],pvalue=pvals[j]) for j in range(N)]
    sm=pd.DataFrame(summ);sm.to_csv(out/'ExtendedData_Table9b_intervention_assurance.csv',index=False)
    pd.DataFrame(long).to_csv(out/'ExtendedData_Table9c_intervention_trials.csv',index=False)
    fig,ax=plt.subplots(1,2,figsize=(11.5,4.1),layout='constrained')
    labels=['Common cause\nno transfer','Directed\ntransfer','Off-target\ncommon cause']
    ax[0].bar(np.arange(3),sm['mean_estimate'],color=['#6c7a89','#28666e','#a65f52'],width=.58)
    ax[0].errorbar(np.arange(3),sm['mean_estimate'],yerr=1.96*sm['sd_estimate']/np.sqrt(N),fmt='none',ecolor='black',capsize=5)
    ax[0].set_xticks(np.arange(3),labels);ax[0].set_ylabel('Mean estimated intervention contrast')
    ax[0].axhline(0,linestyle='--',color='black',lw=.8)
    ax[0].set_title('a  Causal effects distinguish only valid interventions',loc='left',fontsize=10)
    ax[1].bar(np.arange(3),sm['rejection_rate'],color=['#6c7a89','#28666e','#a65f52'],width=.58)
    ax[1].set_xticks(np.arange(3),labels);ax[1].set_ylabel('Fraction with two-sided day-level P < 0.05');ax[1].set_ylim(0,1.0)
    ax[1].axhline(.05,color='black',ls=':',lw=.8)
    ax[1].set_title('b  False identification after exclusion failure',loc='left',fontsize=10)
    fig.savefig(out/'ExtendedData_Figure7_causal_intervention.png',dpi=220,bbox_inches='tight');plt.close(fig)
    (out/'causal_run_manifest.json').write_text(json.dumps({'seed':SEED,'rho':RHO,'gamma':GAMMA,'c':C,'residual_var':V,'natural_n_per_model':50000,'n_campaigns':N,'day_count':12,'units_per_day_per_arm':20,'status':'synthetic exploratory, not biological, not preregistered'},indent=2))
    print(pd.DataFrame(rows).to_string(index=False)); print(sm.to_string(index=False))

if __name__=='__main__':main()
