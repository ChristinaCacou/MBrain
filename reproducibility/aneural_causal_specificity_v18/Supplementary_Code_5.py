#!/usr/bin/env python3
"""Independent-day assurance for an exploratory Stentor design audit.

This is a POST-REGISTRATION theoretical/synthetic analysis only. It does not
alter the publicly locked Stentor protocol, implement its biological experiment,
or supply real Stentor observations. Tested model: Gaussian cell-level paired
changes with an exchangeable random intercept for recording day.

Python >=3.10; numpy, pandas, scipy, matplotlib. Random seed 20261017.
"""
from pathlib import Path
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import nct, t, norm

OUT=Path(__file__).resolve().parent
SEED=20261017
RNG=np.random.default_rng(SEED)
ALPHA=0.05
TOTAL_SHAM=96
N_SIM=2000


def exact_sign_flip_p(day_mean: np.ndarray) -> float:
    """Two-sided exact sign-flip test, assuming independent symmetric day means.

    Only used when day_count<=12; enumeration includes observed assignment.
    """
    k=len(day_mean)
    if k>12: raise ValueError('Exact enumeration disabled for k>12')
    idx=np.arange(2**k, dtype=np.uint32)[:,None]
    shifts=np.arange(k,dtype=np.uint32)[None,:]
    signs=np.where(((idx>>shifts)&1)>0,1.0,-1.0)
    obs=abs(day_mean.sum())
    values=abs(signs@day_mean)
    return float(np.mean(values >= obs-1e-12))


def day_t_power(k,m,icc,effect):
    """Parametric model-based power; effect in cell-level paired-SD units."""
    sigma_day=np.sqrt(icc+(1-icc)/m)
    delta=effect*np.sqrt(k)/sigma_day
    crit=t.ppf(1-ALPHA/2,k-1)
    return float(nct.cdf(-crit,k-1,delta)+nct.sf(crit,k-1,delta))


def analytic_table():
    rows=[]
    for k,icc,effect in product([4,6,8,12,16,24],[0.,.02,.05,.10,.20],[0.,.4,.6]):
        m=TOTAL_SHAM//k
        rows.append(dict(
            independent_days=k,cells_per_day=m,n_sham=k*m,paired_difference_icc=icc,
            standardized_mean_shift=effect,design_effect=1+(m-1)*icc,
            effective_sample_size=TOTAL_SHAM/(1+(m-1)*icc),
            day_mean_t_power=day_t_power(k,m,icc,effect),
            exact_two_sided_signflip_min_p=2./(2**k),
            exact_signflip_can_reject_at_alpha05=(2./(2**k)<ALPHA),
            method_note='Day-level t power conditional on Gaussian independent symmetric day effects; not a biological or empirical estimate'
        ))
    return pd.DataFrame(rows)


def monte_carlo_check():
    rows=[]
    for k in [4,6,8,12]:
        m=TOTAL_SHAM//k
        for icc,effect in product([0.,.05,.10],[0.,.4]):
            effects=RNG.normal(0,np.sqrt(icc),size=(N_SIM,k,1))
            residual=RNG.normal(0,np.sqrt(1-icc),size=(N_SIM,k,m))
            x=effect+effects+residual
            d=x.mean(axis=2)
            day_t=d.mean(axis=1)/(d.std(axis=1,ddof=1)/np.sqrt(k))
            p_day=2*t.sf(abs(day_t),k-1)
            flat=x.reshape(N_SIM,-1)
            cell_t=flat.mean(axis=1)/(flat.std(axis=1,ddof=1)/np.sqrt(TOTAL_SHAM))
            p_cell=2*t.sf(abs(cell_t),TOTAL_SHAM-1)
            # Fully enumerated exact two-sided day-sign-flip test, not a Monte Carlo p-value.
            sign_idx=np.arange(2**k,dtype=np.uint16)[:,None]
            mask=np.arange(k,dtype=np.uint16)[None,:]
            signs=np.where(((sign_idx>>mask)&1)>0,1.,-1.)
            raw=abs(d@signs.T)
            observed=abs(d.sum(axis=1))[:,None]
            p_exact=np.mean(raw>=observed-1e-10,axis=1)
            for name,p in [('naive_cell_t',p_cell),('independent_day_t',p_day),('exact_day_signflip',p_exact)]:
                count=int(np.sum(p<ALPHA))
                phat=count/N_SIM
                lo,hi=norm.ppf(.975),norm.ppf(.975)
                # Wilson MC uncertainty only, NOT a biological confidence interval.
                z=norm.ppf(.975); den=1+z*z/N_SIM
                center=(phat+z*z/(2*N_SIM))/den
                half=z*np.sqrt(phat*(1-phat)/N_SIM+z*z/(4*N_SIM*N_SIM))/den
                rows.append(dict(independent_days=k,cells_per_day=m,n_sham=k*m,
                                 paired_difference_icc=icc,standardized_mean_shift=effect,
                                 inference=name,simulations=N_SIM,number_significant=count,
                                 rejection_fraction=phat,MC_wilson95_lower=center-half,
                                 MC_wilson95_upper=center+half,
                                 mathematical_signflip_min_two_sided_p=2./2**k))
    return pd.DataFrame(rows)


def make_figure(a,b):
    fig,axs=plt.subplots(1,3,figsize=(14,4.1))
    q=a[(a.standardized_mean_shift==.4)]
    for icc in [0.,.02,.05,.10,.20]:
        tdf=q[q.paired_difference_icc==icc]
        axs[0].plot(tdf.independent_days,tdf.day_mean_t_power,'o-',label=f'ICC {icc:.2f}')
    axs[0].set(xlabel='Independent recording days',ylabel='Two-sided day-level t-test power',
               title='a  Fixed 96 cells, standardized effect 0.40',ylim=(0,1.02))
    axs[0].axhline(.8,color='k',linestyle='--',linewidth=.8)
    axs[0].legend(fontsize=7,frameon=False)

    k=np.array([4,5,6,7,8,12])
    floor=2.**(1-k)
    axs[1].step(k,floor,where='mid',marker='o')
    axs[1].axhline(.05,color='k',linestyle='--',linewidth=.8,label='Nominal 0.05')
    axs[1].set(xlabel='Independent day summaries',ylabel='Smallest possible exact two-sided P',
               title='b  Combinatorial resolution of sign flips',yscale='log',ylim=(.0001,.2))
    axs[1].legend(fontsize=8,frameon=False)

    z=b[(b.standardized_mean_shift==0)&(b.paired_difference_icc==.05)]
    for name in ['naive_cell_t','independent_day_t','exact_day_signflip']:
        dd=z[z.inference==name]
        axs[2].plot(dd.independent_days,dd.rejection_fraction,'o-',label=name.replace('_',' '))
    axs[2].axhline(.05,color='k',linestyle='--',linewidth=.8)
    axs[2].set(xlabel='Independent recording days',ylabel='Empirical null rejection fraction',
               title='c  Hierarchical null, ICC=0.05',ylim=(0,.30))
    axs[2].legend(fontsize=7,frameon=False)
    fig.text(.02,.005,'Independent synthetic random-day Gaussian paired differences; no new biological Stentor data. Exact sign-flip requires day-level symmetry.',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1))
    p=OUT/'ExtendedData_Figure5_independent_day_assurance.png'
    fig.savefig(p,dpi=200,bbox_inches='tight')
    plt.close(fig)
    return p


def main():
    a=analytic_table()
    b=monte_carlo_check()
    a.to_csv(OUT/'ExtendedData_Table7_independent_day_assurance.csv',index=False)
    b.to_csv(OUT/'ExtendedData_Table7b_independent_day_montecarlo.csv',index=False)
    f=make_figure(a,b)
    print('Analytic t-test power for standardized effect 0.40 / ICC 0.05')
    print(a[(a.standardized_mean_shift==.4)&(a.paired_difference_icc==.05)][['independent_days','cells_per_day','effective_sample_size','day_mean_t_power','exact_two_sided_signflip_min_p']].to_string(index=False))
    print('\nSimulation under TRUE null, ICC=.05')
    print(b[(b.standardized_mean_shift==0)&(b.paired_difference_icc==.05)][['independent_days','inference','rejection_fraction']].to_string(index=False))
    print(f)

if __name__=='__main__':main()
