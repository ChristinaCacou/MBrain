#!/usr/bin/env python3
"""Re-run and independently audit the numerical benchmarks in the aneural causal-specificity paper.

Synthetic-only; does not process Stentor cell data. Fails if frozen prospective code changes.
"""
import argparse, hashlib, json, platform, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

REGISTERED_SHA='7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151'

def check(flag,msg):
    if not bool(flag): raise AssertionError(msg)

def execute(cmd):
    subprocess.run(cmd,check=True)


def main():
    a=argparse.ArgumentParser()
    a.add_argument('--root',type=Path,default=Path(__file__).parent)
    a.add_argument('--out',type=Path,required=True)
    ns=a.parse_args(); root=ns.root;out=ns.out;out.mkdir(parents=True,exist_ok=True)
    locked=root/'code/Prospective_Stentor_Validation_Analysis.py'
    assert hashlib.sha256(locked.read_bytes()).hexdigest()==REGISTERED_SHA, 'Registered script changed'
    execute([sys.executable,str(root/'code/Supplementary_Code_6.py'),'--locked',str(locked),'--outdir',str(out/'common')])
    execute([sys.executable,str(root/'code/Supplementary_Code_7.py'),'--outdir',str(out/'intervention')])
    common=pd.read_csv(out/'common/ExtendedData_Table8_common_driver_gates.csv').set_index('condition')
    check(bool(common.loc['common_driver','all_three']),'common-driver counterexample did not pass')
    check(not bool(common.loc['no_change','all_three']),'matched no-change control incorrectly passed')
    check(abs(common.loc['common_driver','cross_modal_auc']-.8421223958333334)<1e-12,'AUC mismatch')
    check(abs(common.loc['common_driver','mean_heldout_delta_log_score']-.6795835991691538)<1e-12,'score mismatch')
    inter=pd.read_csv(out/'intervention/ExtendedData_Table9b_intervention_assurance.csv').set_index('scenario')
    check(len(inter)==3,'missing intervention scenario')
    for case,target in [('common_cause',.05075),('coupled',.81375),('offtarget_common_cause',.69225)]:
        check(abs(inter.loc[case,'rejection_rate']-target)<1e-12,case+' rate mismatch')
    # Independently check structural-equation covariance, without using the simulation's natural_samples().
    rho=.6;g=.25; a=(rho-g)/np.sqrt(rho)
    v=1-g*g-a*a-2*g*a*np.sqrt(rho)
    check(v>0,'residual variance negative')
    Mcoeff=np.array([np.sqrt(rho),np.sqrt(1-rho),0.])
    Ecoeff=np.array([g*np.sqrt(rho)+a,g*np.sqrt(1-rho),np.sqrt(v)])
    check(abs(np.dot(Mcoeff,Ecoeff)-rho)<1e-12,'analytic covariance wrong')
    check(abs(np.dot(Ecoeff,Ecoeff)-1)<1e-12,'analytic variance wrong')
    receipt={'audit_type':'deterministic reproduction plus independent analytic covariance check',
            'biological_data':'none','registered_sha256':REGISTERED_SHA,
            'common_driver_three_gate_pass':True,'no_change_three_gate_pass':False,
            'replicated_auc':float(common.loc['common_driver','cross_modal_auc']),
            'replicated_predictive_log_score':float(common.loc['common_driver','mean_heldout_delta_log_score']),
            'rejection_rates':{q:float(inter.loc[q,'rejection_rate']) for q in inter.index},
            'analytic_common_cause_and_directed_covariance':rho,
            'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}
    (out/'verification_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS: reproducibility checks and independent analytic identity')

if __name__=='__main__':main()
