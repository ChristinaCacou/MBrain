#!/usr/bin/env python3
"""Sharp intervention partial identification, conditional on a specified Gaussian family.
Exploratory synthetic mathematical result; no Stentor molecular or physiological data.
"""
from pathlib import Path
import numpy as np, pandas as pd, json
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
RHO=.6;R=np.sqrt(RHO); DELTA_M=1.;DELTA_E=.25;B_CRIT=DELTA_E/R

def theta_of_b(b):
    return (DELTA_E-R*b)/(DELTA_M-b/R)
def residual_variance(t):
    return 1-t*t-(RHO-t)**2/RHO-2*t*(RHO-t)
def nuisance_coefficient(t):
    return (RHO-t)/R
rows=[]
for bound in [0.,.05,.10,.15,.20,.25,.30,B_CRIT]:
    lo=theta_of_b(bound); hi=theta_of_b(-bound)
    rows.append({'max_absolute_Q_shift_SD':bound,'lowest_compatible_direct_effect':lo,
                 'highest_compatible_direct_effect':hi,'zero_direct_effect_compatible':bool(lo<=0<=hi),
                 'min_residual_variance_at_endpoints':min(residual_variance(lo),residual_variance(hi))})
df=pd.DataFrame(rows);df.to_csv(OUT/'sharp_identification_bounds_v20.csv',index=False,float_format='%.12g')
assert np.all(df['min_residual_variance_at_endpoints']>0)
# Full P(Z,M,E) is EXACTLY identical for every b in [-bcrit,bcrit],
# by Gaussian conditional means and covariance. Numerically verify a dense grid.
b_grid=np.linspace(-B_CRIT,B_CRIT,1001)
max_err=0
for b in b_grid:
    t=theta_of_b(b);k=nuisance_coefficient(t);v=residual_variance(t)
    assert v>=0
    delta_m=R*b + (DELTA_M-R*b)
    delta_e=t*delta_m + k*b
    cov_me=t*1 + k*R
    var_e=t*t + k*k + 2*t*k*R +v
    max_err=max(max_err,abs(delta_m-DELTA_M),abs(delta_e-DELTA_E),abs(cov_me-RHO),abs(var_e-1))
assert max_err<1e-12,max_err
# Actuator rank test: calibrated b_i (not merely randomized Z_i) is needed
# to separate theta and the nuisance coefficient k in delta E_i = theta delta M_i + k b_i.
true_theta=.25; true_k=(RHO-true_theta)/R
D=np.array([[1.,0.],[R*.4,.4]])
responses=D@np.array([true_theta,true_k])
recovered=np.linalg.solve(D,responses)
assert np.allclose(recovered,[true_theta,true_k])
# Show exact equivalence empirically at four b values, two intervention arms.
rng=np.random.default_rng(2026102127); rec=[]
for b in [-.2,0.,.2,B_CRIT]:
 t=theta_of_b(b);k=nuisance_coefficient(t);v=residual_variance(t)
 for z in [0,1]:
  n=50000;q=rng.standard_normal(n);merr=rng.standard_normal(n);eerr=rng.standard_normal(n)
  qp=q+b*z
  m=R*qp+np.sqrt(1-RHO)*merr+(DELTA_M-R*b)*z
  e=t*m+k*qp+np.sqrt(v)*eerr
  c=np.cov(m,e)
  rec.append({'b':b,'theta':t,'z':z,'mean_M':np.mean(m),'mean_E':np.mean(e),
              'var_M':c[0,0],'var_E':c[1,1],'cov_ME':c[0,1]})
pd.DataFrame(rec).to_csv(OUT/'numerical_equivalence_grid_check_v20.csv',index=False)
# plot sharp identification interval (all black; reference dashed lines)
x=np.linspace(0,B_CRIT,400);low=np.array([theta_of_b(u) for u in x]);high=np.array([theta_of_b(-u) for u in x]);
fig,ax=plt.subplots(figsize=(7.6,4.7));ax.plot(x,low,color='black',linewidth=2,label='Lower compatible direct effect');ax.plot(x,high,color='black',linewidth=1.5,linestyle='-.',label='Upper compatible direct effect');ax.fill_between(x,low,high,color='grey',alpha=.16);ax.axhline(0,linestyle=':',color='black');ax.axhline(.25,linestyle='--',color='black',alpha=.5);ax.axvline(B_CRIT,linestyle=':',color='black');ax.set(xlabel='Upper bound B on absolute unintended Q shift (SD units)',ylabel='Compatible directed coefficient θ',title='Sharp partial identification within the stated Gaussian family',xlim=(0,B_CRIT),ylim=(-.015,.38));ax.text(.012,.02,'Synthetic identification region, not observed biology',fontsize=8,transform=ax.transAxes);ax.legend(loc='center right',frameon=False,fontsize=8);fig.tight_layout();fig.savefig(OUT/'Figure_S4_sharp_partial_identification.png',dpi=210);plt.close(fig)
report={'rho':RHO,'delta_M':DELTA_M,'delta_E':DELTA_E,'critical_Q_shift':B_CRIT,'maximum_exact_parameter_identity_error':max_err,'n_grid_checks':len(b_grid),'rank_identification_example':{'calibrated_actuator_matrix':D.tolist(),'downstream_contrasts':responses.tolist(),'determinant':float(np.linalg.det(D)),'recovered_theta':float(recovered[0]),'recovered_kappa':float(recovered[1])},'scope':'Specified scalar latent Gaussian structural model with known loadings and fixed residual model; not population-independent causal bounds, and no new Stentor observations.','proof':'For any feasible b, set theta=(e-sqrt(rho)b)/(m-b/sqrt(rho)), kappa=(rho-theta)/sqrt(rho), v=1-theta**2-kappa**2-2*theta*kappa*sqrt(rho). Set Qprime=Q+bZ; M=sqrt(rho)Qprime+sqrt(1-rho)epsM+(m-sqrt(rho)b)Z; E=theta M+kappa Qprime+sqrt(v)epsE. Conditional means (mZ,eZ) and covariance [[1,rho],[rho,1]] identical; therefore exact joint law P(Z,M,E) identical for any b with v>=0.'}
(OUT/'sharp_identification_proof_v20.json').write_text(json.dumps(report,indent=2))
print(df.to_string(index=False,float_format=lambda x:f'{x:.6f}'))
print('max exact identity err',max_err,'critical shift',B_CRIT,'rank example',report['rank_identification_example'])
