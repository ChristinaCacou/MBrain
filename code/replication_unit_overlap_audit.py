"""Post-hoc, exploratory replication-unit / shared-day positivity audit.
This script reads 19-run summary extracted from Doan et al. (2026),
GitHub blob 04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc.
Outputs are secondary analyses of published biological observations.
NO randomization or exchangeability of run/day treatment was verified.
"""
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
src=root/'Independent_Stentor_Followup2_Run_Summary.csv'
df=pd.read_csv(src)
required=['run_id','session_day','condition','test_success','test_n','late_hab_success','late_hab_n']
assert not df[required].isna().any().any()
assert len(df)==19 and df.run_id.nunique()==19
assert set(df.condition)=={'experimental','control'}
assert (df.test_success<=df.test_n).all() and (df.late_hab_success<=df.late_hab_n).all()
df['test_prop']=df.test_success/df.test_n
df['late_prop']=df.late_hab_success/df.late_hab_n
df['test_minus_late']=df.test_prop-df.late_prop
# Exploratory illustrations: treatment labels across entire runs/dates were not assigned by an established randomization.
def mean_difference(data,col):
 return data.loc[data.condition=='experimental',col].mean()-data.loc[data.condition=='control',col].mean()
def exact_reallocation(data,col):
 y=data[col].to_numpy(dtype=float)
 idx=(data.condition=='experimental').to_numpy()
 observed=y[idx].mean()-y[~idx].mean()
 nE=sum(idx);N=len(y);extreme=0;count=0
 for selected in combinations(range(N),nE):
  grp=np.zeros(N,dtype=bool);grp[list(selected)]=True
  effect=y[grp].mean()-y[~grp].mean()
  count+=1
  extreme+=abs(effect)>=abs(observed)-1e-12
 return {'effect':float(observed),'exact_two_sided_p_illustrative':float(extreme/count),'allocations':count,'extreme':int(extreme)}

summaries=[]
sessions=df.groupby(['session_day','condition'],as_index=False).agg(n_runs=('run_id','count'),test_yes=('test_success','sum'),test_n=('test_n','sum'),late_yes=('late_hab_success','sum'),late_n=('late_hab_n','sum'))
sessions['test_prop']=sessions.test_yes/sessions.test_n
sessions['late_prop']=sessions.late_yes/sessions.late_n
sessions['test_minus_late']=sessions.test_prop-sessions.late_prop
sessions.to_csv(root/'ExtendedData_Table12b_day_condition_sessions.csv',index=False)
common_dates=set(sessions.loc[sessions.condition=='experimental','session_day']) & set(sessions.loc[sessions.condition=='control','session_day'])
within=df[df.session_day.isin(common_dates)].copy()
assert len(common_dates)==1 and len(within)==4
for col in ['test_prop','test_minus_late']:
 for level,data in [('all_runs',df),('day_condition_sessions',sessions),('shared_date_runs_only',within)]:
  z=exact_reallocation(data,col)
  z.update({'comparison':col,'analysis_level':level,'n_units':len(data),'n_experimental':int(sum(data.condition=='experimental')),'n_control':int(sum(data.condition=='control')),'exchangeability_verified':False})
  summaries.append(z)
summary=pd.DataFrame(summaries)
summary.to_csv(root/'ExtendedData_Table12_cluster_and_overlap_sensitivity.csv',index=False)
# Partial identification of a DATE-STANDARDIZED descriptive contrast.
# The absent arm's binary response probability may lie anywhere in [0,1].
# These sharp algebraic bounds condition on observed date-arm proportions as fixed,
# ignore sampling uncertainty and DO NOT identify a causal treatment effect.
bydate=[]
for day,chunk in sessions.groupby('session_day',sort=True):
    e=chunk[chunk.condition=='experimental']
    c=chunk[chunk.condition=='control']
    if len(e) and len(c):
        l=u=float(e.test_prop.iloc[0]-c.test_prop.iloc[0]); support='both'
    elif len(e):
        v=float(e.test_prop.iloc[0]);l=v-1;u=v; support='paired_only'
    else:
        v=float(c.test_prop.iloc[0]);l=-v;u=1-v; support='control_only'
    bydate.append({'session_day':day,'condition_support':support,'conditional_lower_bound':l,'conditional_upper_bound':u})
bounds=pd.DataFrame(bydate)
bounds.to_csv(root/'ExtendedData_Table12c_Missing_Arm_Bounds.csv',index=False)
mean_bounds=[float(bounds.conditional_lower_bound.mean()),float(bounds.conditional_upper_bound.mean())]
assert mean_bounds[0] <=0 <= mean_bounds[1]
# Sample-weighted aggregate only (descriptive, cells not the independent replication units).
ex=df[df.condition=='experimental'];ct=df[df.condition=='control']
cell=(ex.test_success.sum()/ex.test_n.sum())-(ct.test_success.sum()/ct.test_n.sum())
day_ct=pd.crosstab(df.session_day,df.condition)
assert day_ct.shape==(8,2)
status={
 'source_repository':'https://github.com/nhithuydoan/associative_learning_Stentor',
 'source_csv':'data/dataset_followup_2.csv',
 'source_git_blob_sha':'04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc',
 'data_type':'independent published behavioral contraction observations; not RNA-seq',
 'n_cells':int(df.test_n.sum()),'n_runs':len(df),'n_dates':df.session_day.nunique(),
 'n_date_by_condition_sessions':len(sessions),'dates_both_conditions':sorted(common_dates),
 'n_dates_with_both_conditions':len(common_dates),
 'date_support_table':day_ct.to_dict(orient='index'),
 'cell_weighted_difference':float(cell),
 'date_standardized_binary_response_conditional_bounds_equal_date_weights':mean_bounds,
 'bounds_ignore_sampling_uncertainty':True,
 'not_identified':'Recovery RNA time-adjusted transcript contrasts; shared-day causal treatment effect beyond the observed date; individual consciousness',
 'interpretation':'Observational comparisons; exact label reallocations are diagnostic arithmetic and not validated randomization tests.'
}
(root/'independent_biology_overlap_status.json').write_text(json.dumps(status,indent=2))
# Published source audit for the molecular datasets: these values are metadata facts, not downloaded gene-level reanalyses.
source_audit=[
 {'dataset':'Rajan 2026 polyculture','public_source':'https://datadryad.org/dataset/doi:10.5061/dryad.z34tmpgrm','target':'untreated recovery pools at 10 and 90 minutes','confirmed':'six archived untreated-control samples named in public metadata; excluded from released 15-sample normalized matrix','independent_requantification':'NOT COMPLETE','size_issue':'rna_exp2.zip ~209.83 GB; raw reads not available to present local analysis'},
 {'dataset':'Miao 2026 transcriptomics','public_source':'https://github.com/yuemiao1017/Stentor_genome_habituation_transcriptomics','target':'21 pools across seven time points; original authors report 341 time-varying genes','confirmed':'GitHub repository documents analysis workflow; data/ currently contains README.md without processed count matrix','independent_requantification':'NOT COMPLETE','size_issue':'SRA reads and new reference would need aligned reprocessing and gene crosswalk'},
 {'dataset':'Doan 2026 common-test behavior','public_source':'https://github.com/nhithuydoan/associative_learning_Stentor','target':'independent cell contraction observations','confirmed':'19-run summary from published follow-up dataset 2','independent_requantification':'COMPLETE FOR BEHAVIOR ONLY','size_issue':'not an assay of RNA expression or recovery'}]
pd.DataFrame(source_audit).to_csv(root/'Source_Qualification_v14.csv',index=False)
# figure
fig,axs=plt.subplots(1,3,figsize=(13.6,4.0))
order=sorted(df.session_day.unique());
for i,day in enumerate(order):
 es=day_ct.loc[day,'experimental'];cs=day_ct.loc[day,'control']
 if es:axs[0].scatter(i+.5,1,s=45+35*es,marker='o',edgecolors='black',facecolors='white')
 if cs:axs[0].scatter(i+.5,0,s=45+35*cs,marker='s',edgecolors='black',facecolors='black')
axs[0].set_yticks([0,1],['Control','Paired']);axs[0].set_xticks(np.arange(len(order))+.5,[d[5:] for d in order],rotation=45,ha='right');axs[0].set_ylim(-.45,1.45);axs[0].set_title('a  Recording-day overlap')
axs[0].set_xlabel('Date (2026; MM_DD)');axs[0].grid(axis='x',alpha=.15)
for cond,marker in [('experimental','o'),('control','s')]:
 sub=df[df.condition==cond]
 x0=0 if cond=='control' else 1
 axs[1].scatter(x0+np.linspace(-.10,.10,len(sub)),sub.test_prop,marker=marker,facecolors='none' if cond=='experimental' else 'black',edgecolors='black',s=36)
 axs[1].hlines(sub.test_prop.mean(),x0-.18,x0+.18,linewidth=2,colors='black')
axs[1].set_xticks([0,1],['Control','Paired']);axs[1].set_ylim(0,.52);axs[1].set_ylabel('Contraction probability');axs[1].set_title('b  All 19 runs')
axs[1].text(.5,.47,'Across-run difference: +8.57 pp',ha='center',fontsize=9)
for cond,marker in [('experimental','o'),('control','s')]:
 sub=within[within.condition==cond];x0=0 if cond=='control' else 1
 axs[2].scatter(x0+np.array([-.065,.065]),sub.test_prop,marker=marker,facecolors='none' if cond=='experimental' else 'black',edgecolors='black',s=65)
 axs[2].hlines(sub.test_prop.mean(),x0-.18,x0+.18,lw=2,colors='black')
axs[2].set_xticks([0,1],['Control','Paired']);axs[2].set_ylim(0,.52);axs[2].set_title('c  Sole shared date: May 11')
axs[2].text(.5,.47,'Within-date difference: +2.86 pp',ha='center',fontsize=9)
axs[2].text(.5,.42,'P = 1.00 (6 hypothetical allocations)',ha='center',fontsize=8)
fig.suptitle('Independent Stentor behavior: treatment-by-date non-overlap constrains attribution',fontsize=12)
fig.text(.5,.015,'Each point in b,c is one run; treatment labels were not verified randomized. This is not molecular replication.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.05,1,.92]);fig.savefig(root/'ExtendedData_Figure11_Treatment_Date_Overlap.png',dpi=250)
print(json.dumps(status,indent=2));print(summary.to_string(index=False))
