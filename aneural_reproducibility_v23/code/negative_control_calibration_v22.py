#!/usr/bin/env python3
"""Sensitivity bounds from imperfect negative-control proxy shifts.
All numbers are hypothetical. Applies to treatment-contrast data in fixed
linear-Gaussian scalar-Q structural family, not new Stentor measurements.
Note: does not identify full joint distributions involving measured W.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RHO=.60
R=np.sqrt(RHO)
M=1.0
E=.25
BASE_B=.35
OUT=Path(__file__).resolve().parent/'results_v22'
OUT.mkdir(exist_ok=True)

def theta(b,e=E,m=M):
    return (e-R*b)/(m-b/R)

def admissible_residual(t):
    return (1-RHO)*(1-t*t/RHO)

def identified_bounds(inter):
    L,U=inter
    if L>U:
        return (np.nan,np.nan)
    assert L>-R*M and U<R*M, 'Denominator sign not controlled'
    vals=[theta(L),theta(U)]
    assert min(admissible_residual(t) for t in vals)>0
    return min(vals),max(vals)

def proxy_constraint(w,measurement_error,direct_effect_bound,loading):
    """w is nominal proxy contrast; measurement_error is an ABSOLUTE contrast uncertainty budget."""
    assert loading>0 and measurement_error>=0 and direct_effect_bound>=0
    return ((w-measurement_error-direct_effect_bound)/loading,
            (w+measurement_error+direct_effect_bound)/loading)

def intersect(*ranges):
    return max(z[0] for z in ranges),min(z[1] for z in ranges)

cases=[
    {'case':'No proxy calibration','w':None,'error':None,'D':None,'loading':None},
    {'case':'Ideal unchanged proxy','w':0.,'error':.02,'D':0.,'loading':.8},
    {'case':'Bounded proxy leakage','w':0.,'error':.02,'D':.03,'loading':.8},
    {'case':'Larger allowable leakage','w':0.,'error':.02,'D':.16,'loading':.8},
    {'case':'Proxy cannot exclude zero','w':0.,'error':.02,'D':.24,'loading':.8},
    {'case':'Nonzero proxy contrast','w':.04,'error':.02,'D':.03,'loading':.8},
]
rows=[]
for c in cases:
    prior=(-BASE_B,BASE_B)
    if c['w'] is None:
        bi=prior
    else:
        bi=intersect(prior,proxy_constraint(c['w'],c['error'],c['D'],c['loading']))
    lo,hi=identified_bounds(bi)
    rows.append({**c,'b_lower':bi[0],'b_upper':bi[1],
                 'theta_lower':lo,'theta_upper':hi,
                 'zero_directed_effect_compatible':bool(lo<=0<=hi)})
results=pd.DataFrame(rows)
results.to_csv(OUT/'proxy_calibration_scenarios_v22.csv',index=False,float_format='%.10g')
# Exact threshold at which proxy contrast w=0 +/- epsilon ceases to rule out zero directed effect.
zero_b=E/R
Dcritical=.8*zero_b-.02
assert np.isclose(Dcritical,.238198889747161)
assert results.iloc[3]['theta_lower']>0
assert results.iloc[4]['theta_lower']<0

# Two proxies: overlapping bounds shrink, incompatible bounds reject a *single Q* calibrated model.
coherent=intersect((-BASE_B,BASE_B),proxy_constraint(.02,.01,.02,.8),proxy_constraint(.035,.01,.015,1.2))
incoherent=intersect((-BASE_B,BASE_B),proxy_constraint(.10,.01,.01,.8),proxy_constraint(-.10,.01,.01,1.2))
assert coherent[0] < coherent[1]
assert incoherent[0]>incoherent[1]

# Noiseless symbolic verification: θ(b) matches observed M and E for the full span
# by setting Q-shift b, Q-loading k=(rho-theta)/sqrt(rho), direct M shift (m-sqrt(rho)*b).
grid=np.linspace(-BASE_B,BASE_B,2001)
residuals=[]
for b in grid:
    t=theta(b);k=(RHO-t)/R
    assert admissible_residual(t)>0
    measured_delta_m=(M-R*b)+R*b
    measured_delta_e=t*measured_delta_m+k*b
    residuals.extend((abs(measured_delta_m-M),abs(measured_delta_e-E)))
assert max(residuals)<2e-15

fig,ax=plt.subplots(figsize=(8,4.8))
labels=[r['case'] for r in rows]
for j,r in enumerate(rows):
    y=len(rows)-1-j
    a,b=r['theta_lower'],r['theta_upper']
    ax.hlines(y,a,b,linewidth=4,color=str(.12+.11*j))
    ax.plot([a,b],[y,y],ls='none',marker='|',markersize=12,color='black')
ax.axvline(0,ls='--',color='black',linewidth=1,label='No directed effect')
ax.set_yticks(list(range(len(rows)))[::-1], labels=labels)
ax.set_xlabel('Compatible directed coefficient (unitless, hypothetical Gaussian model)')
ax.set_xlim(-.07,.42)
ax.set_title('Calibrating proxy leakage contracts the causal identification set')
ax.grid(axis='x',color='.85',alpha=.8)
ax.set_axisbelow(True)
fig.tight_layout();fig.savefig(OUT/'Figure_S7_proxy_calibration_v22.png',dpi=220);plt.close(fig)

# Sensitivity frontier: proxy leakage allowances at e=0.25, m=1, w=0 +/-0.02.
Ds=np.linspace(0,.30,601)
bounds=[]
for D in Ds:
    inter=intersect((-BASE_B,BASE_B),proxy_constraint(0.,.02,D,.8))
    bounds.append(identified_bounds(inter))
bounds=np.array(bounds)
fig,ax=plt.subplots(figsize=(7.4,4.5))
ax.fill_between(Ds,bounds[:,0],bounds[:,1],color='.7',alpha=.65)
ax.plot(Ds,bounds[:,0],color='black',lw=1.3,label='Lower compatible coefficient')
ax.plot(Ds,bounds[:,1],color='.3',ls='--',lw=1.3,label='Upper compatible coefficient')
ax.axvline(Dcritical,ls=':',color='black',lw=1.2,label='Minimum D admitting θ = 0')
ax.axhline(0,color='black',lw=1)
ax.set(xlabel='Allowed direct Z effect on proxy, D (hypothetical proxy units)',
       ylabel='Compatible directed coefficient',title='A selectivity threshold for proxy-based exclusion')
ax.legend(frameon=True,facecolor='white',edgecolor='.7',fontsize=8,loc='lower left')
fig.tight_layout();fig.savefig(OUT/'Figure_S8_proxy_leakage_frontier_v22.png',dpi=220);plt.close(fig)
receipt={'analysis':'hypothetical contrast-based negative-control proxy calibration',
         'rho':RHO,'target_m':M,'response_e':E,'prior_b_bound':BASE_B,
         'minimum_b_for_theta_zero':zero_b,'lambda':.8,
         'proxy_contrast_assumed':0.,'proxy_contrast_absolute_error':.02,
         'minimum_direct_proxy_leakage_to_admit_zero':Dcritical,
         'coherent_two_proxy_b_range':coherent,
         'incompatible_two_proxy_intersection':incoherent,
         'max_structural_identity_discrepancy':max(residuals),
         'not_empirical':True,'note':'Use proxy mean contrasts only, not the full (M,E,W) covariance.'}
(OUT/'proxy_calibration_verification_v22.json').write_text(json.dumps(receipt,indent=2))
print(results.to_string(index=False));print('\n',json.dumps(receipt,indent=2))
# With interval-valued observed contrasts, propagate the proxy-derived b constraint by eight-corner search.
def rectangle_theta_envelope(b_interval,e_interval=(.20,.30),m_interval=(.95,1.05)):
    L,U=b_interval
    assert L<=U and m_interval[0]>U/R
    v=[theta(b,e,m) for b in (L,U) for e in e_interval for m in m_interval]
    assert all(admissible_residual(t)>=0 for t in v)
    return min(v),max(v)
row2=[]
for D in [0.,.03,.16,.18,.19,.24]:
    br=intersect((-BASE_B,BASE_B),proxy_constraint(0.,.02,D,.8))
    low,high=rectangle_theta_envelope(br)
    row2.append({'direct_proxy_shift_bound_D':D,'b_lower':br[0],'b_upper':br[1],
                 'theta_lower':low,'theta_upper':high,'zero_compatible':low<=0<=high})
rect=pd.DataFrame(row2)
rect.to_csv(OUT/'proxy_plus_contrast_uncertainty_v22.csv',index=False,float_format='%.10g')
wide_critical=.8*(.20/R)-.02
assert rect.iloc[3]['theta_lower']>0 and rect.iloc[4]['theta_lower']<0
with (OUT/'proxy_calibration_verification_v22.json').open() as f: rec=json.load(f)
rec['threshold_with_hypothetical_outcome_range_0.20_to_0.30']=wide_critical
(OUT/'proxy_calibration_verification_v22.json').write_text(json.dumps(rec,indent=2))
print('\nEstimation uncertainty:\n',rect.to_string(index=False),'\nDcritical with outcome uncertainty:',wide_critical)
