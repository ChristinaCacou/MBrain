#!/usr/bin/env python3
"""Reproduce v21 illustrative robustness calculations. Assumes synthetic linear-Gaussian model.
No new biological data or biological calibration. 50k Gaussian contrast draws per full-rank design.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
RHO=0.6; R=np.sqrt(RHO); S=0.02
OUT=Path(__file__).resolve().parent.parent/'results';OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(2026102131)
def theta(e,m,b):return (e-R*b)/(m-b/R)
def variance(t):return (1-RHO)*(1-t*t/RHO)
crit=.25/R
records=[]
for B in [0,.05,.1,.15,.2,.25,.3,crit]:
 endpoints=[theta(.25,1.,z) for z in (-B,B)]
 corners=[theta(e,m,z) for e in (.2,.3) for m in (.95,1.05) for z in (-B,B)]
 assert min(variance(v) for v in corners)>0
 vals=theta(rng.uniform(.2,.3,5000),rng.uniform(.95,1.05,5000),rng.uniform(-B,B,5000))
 assert min(vals)>=min(corners)-1e-12 and max(vals)<=max(corners)+1e-12
 records.append(dict(B_Q_shift_SD=B,ideal_lower=min(endpoints),ideal_upper=max(endpoints),budget_lower=min(corners),budget_upper=max(corners),budget_contains_zero=min(corners)<=0<=max(corners),minimum_residual_variance=min(variance(v) for v in corners)))
pd.DataFrame(records).to_csv(OUT/'estimation_calibration_envelope_v21.csv',index=False,float_format='%.12g')
truth=np.array([.25,(RHO-.25)/R]);designs=[]
for name,m2,b2 in [('Separated',.32,.4),('Moderately separated',.8,.2),('Nearly parallel',.8,.081),('Exactly parallel',.8,.08)]:
 A=np.array([[1.,.1],[m2,b2]]); det=np.linalg.det(A)
 if abs(det)<1e-10:
  designs.append(dict(design=name,m1=1.,b1=.1,m2=m2,b2=b2,determinant=0.,condition_number=np.inf,SE_theta_analytic=np.nan,MC_SD_theta=np.nan,number_campaigns=0)); continue
 se=S*np.sqrt(.1**2+b2**2)/abs(det)
 E=(A@truth)[:,None]+rng.normal(0,S,size=(2,50000))
 draws=np.linalg.solve(A,E);mc=np.std(draws[0],ddof=1)
 assert abs(mc/se-1)<.025
 designs.append(dict(design=name,m1=1.,b1=.1,m2=m2,b2=b2,determinant=float(det),condition_number=float(np.linalg.cond(A)),SE_theta_analytic=float(se),MC_SD_theta=float(mc),number_campaigns=50000))
pd.DataFrame(designs).to_csv(OUT/'actuator_conditioning_v21.csv',index=False,float_format='%.12g')
plt.rcParams.update({'font.family':'serif','font.size':10})
Bs=np.linspace(0,crit,401); lower=np.array([theta(.25,1.,B) for B in Bs]);upper=np.array([theta(.25,1.,-B) for B in Bs]); l=[];u=[]
for B in Bs:
 corners=[theta(e,m,z) for e in (.2,.3) for m in (.95,1.05) for z in (-B,B)];l.append(min(corners));u.append(max(corners))
fig,ax=plt.subplots(figsize=(7.2,4.7));ax.fill_between(Bs,l,u,color='0.8',alpha=.5,label='Hypothetical estimation + calibration envelope');ax.fill_between(Bs,lower,upper,color='0.45',alpha=.25,label='Ideal fixed-contrast identification set');ax.plot(Bs,lower,'k-',lw=1.5);ax.plot(Bs,upper,'k--',lw=1.5);ax.plot(Bs,l,color='.35',ls=':',lw=1);ax.plot(Bs,u,color='.35',ls=':',lw=1);ax.axhline(0,color='black',lw=1);ax.axvline(.2/R,color='.2',lw=1,ls='-.');ax.annotate('Zero becomes compatible\\nunder assumed uncertainty',xy=(.2/R,0),xytext=(.12,-.067),arrowprops={'arrowstyle':'->','color':'black'},fontsize=8);ax.set(xlim=(0,crit),ylim=(-.11,.43),xlabel='Maximum absolute common-driver shift B (SD units)',ylabel='Directed coefficient compatible with assumptions',title='Estimand uncertainty enlarges model-conditioned identification sets');ax.legend(frameon=False,loc='upper left',fontsize=8);fig.tight_layout();fig.savefig(OUT/'Figure_S5_uncertainty_envelope_v21.png',dpi=220);plt.close(fig)
fig,ax=plt.subplots(figsize=(7.2,4.4)); v=designs[:3]; x=np.arange(3);ax.bar(x,[d['SE_theta_analytic'] for d in v],.35,facecolor='.8',edgecolor='black',label='Exact analytical SD');ax.scatter(x,[d['MC_SD_theta'] for d in v],color='black',marker='x',s=42,label='50,000-replicate Monte Carlo SD');ax.axhline(.25,color='black',ls='--',lw=.9,label='True directed coefficient = 0.25');ax.set_xticks(x,labels=[d['design'] for d in v]);ax.set_ylabel('SD of recovered directed coefficient');ax.set_title('Near-singular actuator designs amplify downstream contrast noise');ax.legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(OUT/'Figure_S6_actuator_conditioning_v21.png',dpi=220);plt.close(fig)
print(pd.DataFrame(records).to_string(index=False))
print(pd.DataFrame(designs).to_string(index=False))
