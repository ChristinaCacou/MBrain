#!/usr/bin/env python3
"""Post-registration adversarial common-driver stress test.

A purely mechanistic common-drive null changes two locally independent sensorimotor
pathways. It may pass the three published Stentor confirmatory *statistical*
gates without constituting a model of integrated information or subjective experience.

All observations are synthetic; the unmodified locked analysis is imported only
from a matching SHA256 file. No code or threshold in that frozen pipeline is edited.

Run:
  python Supplementary_Code_6.py --locked /path/Prospective_Stentor_Validation_Analysis.py --outdir output
Dependencies: numpy pandas scipy scikit-learn statsmodels matplotlib
"""
import argparse, hashlib, importlib.util, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FROZEN_SHA = '7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151'
MASTER_SEED = 2026101706
N_DAYS = 4
CELLS_PER_DAY = 24
N_TIME = 100
N_SHAPE = 32

def load_frozen(path):
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if h != FROZEN_SHA:
        raise RuntimeError('The registered analysis file hash does not match: '+h)
    spec=importlib.util.spec_from_file_location('locked_stentor_analysis',str(path))
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def generate(condition, destination):
    """Two independent local response modules conditioned on ONE passive drive q.

    q is an exogenously specified metabolic/ionic-gain state shared by independent
    effector modules. Different locally nonlinear mode sensitivities to q deform
    trajectories even after unit-energy normalization. No m->e / e->m linkage.
    """
    dest=destination / 'synthetic_trajectories';dest.mkdir(parents=True,exist_ok=True)
    r=np.random.default_rng(MASTER_SEED)
    t=np.linspace(0,1,N_TIME)
    modes=np.stack([np.exp(-((t-.36)/.15)**2),
                    np.exp(-3*t)*np.sin(6*np.pi*t),
                    np.exp(-2*t)*np.cos(9*np.pi*t)],axis=1)
    modes[:20]*=.05
    rows=[]
    for day in range(1,N_DAYS+1):
        for j in range(CELLS_PER_DAY):
            cid=f'SYN_D{day}_C{j:02d}'
            raw_hab=0.72+r.normal(0,.08)
            qhab=float(np.clip(raw_hab,.42,.95)) if condition=='common_driver' else 1.0
            for mod in ['MECH','ELEC']:
                # Crucial: independent local loadings per modality and independent noise.
                # The two modules interact only through the exogenous nuisance q.
                W=r.normal(size=(3,N_SHAPE)); W/=np.linalg.norm(W,axis=1,keepdims=True)
                for state in ['BASE','HAB']:
                    q=1.0 if state=='BASE' else qhab
                    weights=np.array([q,1.15*q**3,1.15*q**5])
                    x=modes@(weights[:,None]*W)
                    x+=r.normal(size=x.shape)*(.08*q)
                    fn=f'{cid}_{state}_{mod}.npz'
                    np.savez_compressed(dest/fn,trajectory=x)
                    rows.append(dict(cell_id=cid,replicate_day=f'D{day}',
                        treatment_code='SHAM',state_code=state,modality_code=mod,
                        trial_id=f'{state}_{mod}',trajectory_file=fn,excluded=False,
                        exclusion_reason='none'))
    meta=destination/'SYNTHETIC_metadata.csv'
    pd.DataFrame(rows).to_csv(meta,index=False)
    return meta,dest

def analyze(locked,case,out):
    wd=out/case;wd.mkdir(parents=True,exist_ok=True)
    meta,tdir=generate(case,wd)
    result=locked.run_pipeline(meta,tdir,wd/'locked_outputs')
    # summary only, never identify this as empirical Stentor data.
    q={
      'condition':case,
      'n_sham_cells':N_DAYS*CELLS_PER_DAY,
      'primary_base_mean':result['primary']['baseline_mean'],
      'primary_hab_mean':result['primary']['habituated_mean'],
      'paired_dz':result['primary']['paired_dz'],
      'primary_p':result['primary']['paired_t_p'],
      'cross_modal_auc':result['cross_modal'].get('auc'),
      'cross_modal_perm_p':result['cross_modal'].get('permutation_p'),
      'cross_modal_boot_ci_low':result['cross_modal'].get('cluster_bootstrap_95ci',[np.nan,np.nan])[0],
      'mean_heldout_delta_log_score':result['model_score'].get('mean_delta_log_score'),
      'model_score_signflip_p':result['model_score'].get('signflip_p'),
      'model_score_boot_ci_low':result['model_score'].get('cluster_bootstrap_95ci',[np.nan,np.nan])[0],
      'pass_geometry':result['primary']['passes_primary_difference_gate'],
      'pass_crossmodal':result['cross_modal'].get('passes_cross_modal_gate',False),
      'pass_predictive':result['model_score'].get('passes_model_score_gate',False),
      'all_three':result['decision']['local_scalar_gain_null_rejected'],
    }
    return q

def figure(summary,out):
    fig,axs=plt.subplots(1,3,figsize=(11.6,3.4),layout='constrained')
    a,b=summary.iloc[0],summary.iloc[1]
    xs=np.array([0,1]);d=np.array([a.primary_base_mean,a.primary_hab_mean]);n=np.array([b.primary_base_mean,b.primary_hab_mean]);
    axs[0].plot(xs,d,'o-',label='Common physiological drive',lw=2)
    axs[0].plot(xs,n,'s--',label='No-change control',lw=1.6)
    axs[0].set_xticks(xs,['Baseline','Post training']);axs[0].set_ylabel('Electrical trajectory participation ratio')
    axs[0].legend(fontsize=7,frameon=False);axs[0].set_title('a  Normalized response geometry',loc='left',fontsize=10)
    for j,metric in enumerate(['cross_modal_auc','mean_heldout_delta_log_score']):
        yy=np.array([a[metric],b[metric]])
        axs[j+1].bar([0,1],yy,width=.52)
        axs[j+1].set_xticks([0,1],['Common drive','No change'],rotation=15,ha='right')
        axs[j+1].set_title(['b  Transfer from mechanical to electrical','c  Leave-day-out N1 versus N0'][j],loc='left',fontsize=9)
        axs[j+1].set_ylabel(['Cross-modal AUC','Mean held-out Δ log-score'][j])
        if j==0:axs[j+1].axhline(.5,color='black',lw=.7,ls=':')
        else:axs[j+1].axhline(0,color='black',lw=.7,ls=':')
    fig.savefig(out/'ExtendedData_Figure6_common_driver.png',dpi=210,bbox_inches='tight')
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--locked',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    args=ap.parse_args();out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    locked=load_frozen(args.locked)
    summary=pd.DataFrame([analyze(locked,'common_driver',out),analyze(locked,'no_change',out)])
    summary.to_csv(out/'ExtendedData_Table8_common_driver_gates.csv',index=False)
    figure(summary,out)
    (out/'RUN_RECEIPT.json').write_text(json.dumps({'status':'synthetic only','seed':MASTER_SEED,'n_days':N_DAYS,
        'cells_per_day':CELLS_PER_DAY,'registered_code_sha256':FROZEN_SHA},indent=2))
    print(summary.to_string(index=False))

if __name__=='__main__':main()
