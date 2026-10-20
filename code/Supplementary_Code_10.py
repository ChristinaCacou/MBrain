#!/usr/bin/env python3
"""Post hoc missing-time-control sensitivity for previously selected Stentor panel.

Input: 14 transcript IDs x 15 sample pooled author-normalized counts extracted
and source-verified in the previous biological audit, NOT new wet-lab data.
This code does not alter the prospective Stentor validation protocol and
cannot substitute for omitted untreated recovery-time samples.
"""
from pathlib import Path
import itertools,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
counts=pd.read_csv(P/'previous_14genes.tsv',sep='\t',header=None,index_col=0)
assert counts.shape==(14,15)
SAMPLES=['DR4876_01','DR4876_02','DR4876_03','DR4876_04','DR4876_05','DR4876_06','DR4876_07','DR4876_08','DR4876_09','DR4876_10','DR4876_12','DR4876_14','DR4876_16','DR4876_18','DR4876_20']
GROUPS={'T':['DR4876_01','DR4876_04','DR4876_07'],'U':['DR4876_03','DR4876_06','DR4876_09'],'F10':['DR4876_10','DR4876_12','DR4876_14'],'F90':['DR4876_16','DR4876_18','DR4876_20']}
GLOBAL_MEAN=[7.41346,7.36912,7.42708,7.41601,7.40137,7.39913,7.44915,7.36384,7.42764,7.40139,7.40017,7.42029,7.41061,7.44015,7.39559]
counts.columns=SAMPLES
scores=np.log2(1+counts).mean(axis=0)
global_centers=pd.Series(GLOBAL_MEAN,index=SAMPLES)
rows=[]
boot={}
rng=np.random.default_rng(241018)
for flavor,vec in [('raw_log2',scores),('global_centered_log2',scores-global_centers)]:
    m={g:float(vec.loc[s].mean()) for g,s in GROUPS.items()}
    tr=m['T']-m['F90'];baseline=m['T']-m['U'];point_thresh=tr
    # sensitivity parameter d = expected untreated 90 min minus untreated immediate (unobserved).
    # DID=(F90-T)-d; negative represents a training-specific decrease.
    samples_T=vec.loc[GROUPS['T']].to_numpy();samples_F=vec.loc[GROUPS['F90']].to_numpy()
    B=80000
    draws_T=rng.choice(samples_T,size=(B,3),replace=True).mean(axis=1)
    draws_F=rng.choice(samples_F,size=(B,3),replace=True).mean(axis=1)
    potential_recovery=draws_T-draws_F
    ci=np.quantile(potential_recovery,[.025,.975])
    boot[flavor]={'observed_T_minus_F90':tr,'bootstrap_95_CI_T_minus_F90':ci.tolist(),
                  'drift_sign_flip_point':-tr,
                  'drift_magnitude_to_erase_point':tr,
                  'drift_magnitude_to_include_zero_in_resampling_interval':float(ci[0]),
                  'observed_attenuation_vs_U':tr/baseline,
                  'n_pools_per_condition':3}
    for d in [-.50,-.40,-.35,-.30,-.25,-.20,-.10,0,.10,.20,.30,.40,.50]:
        rows.append({'measurement':flavor,'counterfactual_untreated_90min_drift_log2':d,
            'observed_F90_minus_T':-tr,
            'time_control_adjusted_DID':-tr-d,
            'sampling_CI_low_given_drift':float(-ci[1]-d),
            'sampling_CI_high_given_drift':float(-ci[0]-d),
            'point_direction_training_specific_decline':bool(-tr-d<0),
            'sampling_interval_excludes_zero':bool((-ci[1]-d>0) or (-ci[0]-d<0)),
            'reference_T_minus_U':baseline})
    for BND in [0,0.10,0.20,0.25,0.30,0.35,0.40,0.50]:
        # Set identified conditional point-effect range given |untreated drift| <= BND.
        rows.append({'measurement':flavor,'counterfactual_untreated_90min_drift_log2':np.nan,
            'max_abs_untreated_drift_B':BND,
            'identified_effect_lower':float(-tr-BND),
            'identified_effect_upper':float(-tr+BND),
            'sampling_plus_sensitivity_lower':float(-ci[1]-BND),
            'sampling_plus_sensitivity_upper':float(-ci[0]+BND)})
full=pd.DataFrame(rows)
full.to_csv(P/'ExtendedData_Table11_missing_time_control_sensitivity.csv',index=False,float_format='%.9g')
(P/'ExtendedData_Table11b_sensitivity_summary.json').write_text(json.dumps(boot,indent=2))
fig,axs=plt.subplots(1,2,figsize=(10.9,3.65),dpi=170)
for ax,(flavor,name) in zip(axs,[('raw_log2','Published normalized counts'),('global_centered_log2','Transcriptome-centered counts')]):
    tr=boot[flavor]['observed_T_minus_F90'];ci=boot[flavor]['bootstrap_95_CI_T_minus_F90'];xs=np.linspace(-.55,.55,221)
    ax.axhline(0,c='black',lw=1.2)
    ax.axvline(-tr,c='black',lw=1,ls=':')
    ax.fill_between(xs,-ci[1]-xs,-ci[0]-xs,color='#dddddd',label='Pool-bootstrap 95% interval at assumed drift')
    ax.plot(xs,-tr-xs,c='black',lw=2,label='Time-control adjusted contrast')
    ax.scatter([-tr],[0],color='black',s=25,zorder=4)
    ax.set_title(name,fontsize=11,loc='left')
    ax.set_xlabel('Unobserved change in untreated pools over 90 min (log$_2$ units)',fontsize=9)
    ax.set_ylabel('Adjusted F90 − T19 effect (log$_2$ units)',fontsize=9)
    ax.grid(alpha=.17)
    ax.set_xlim(-.55,.55);ax.set_ylim(-.85,.65)
    ax.text(-.52,.53,'Untreated time drift is unmeasured',fontsize=8)
axs[0].legend(loc='lower left',fontsize=7,frameon=False)
fig.suptitle('How unobserved time-matched controls change the apparent Stentor molecular recovery signal',fontsize=12,y=1.03)
fig.tight_layout();fig.savefig(P/'ExtendedData_Figure9_time_control_tipping.png',dpi=220,bbox_inches='tight');plt.close(fig)
print(json.dumps(boot,indent=2))
