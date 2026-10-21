#!/usr/bin/env python3
"""Reproduce the v18 synthetic analyses and post-hoc four-day audit.
Separate frozen prospective source remains byte-identical; synthetic tests are not Stentor observations.
No new confirmatory experiment was run or modified by this package.
"""
import subprocess,sys,hashlib,json,argparse,platform
from pathlib import Path
import pandas as pd
EXPECTED_LOCK='7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151'
BASE=Path(__file__).resolve().parent

def run(cmd):
 print('RUN:', ' '.join(str(c) for c in cmd),flush=True)
 subprocess.run(list(map(str,cmd)),check=True,env=None)

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--skip-causal',action='store_true');p.add_argument('--behavior-csv',type=Path)
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 lock=BASE/'frozen_Prospective_Stentor_Validation_Analysis.py'
 h=hashlib.sha256(lock.read_bytes()).hexdigest()
 if h!=EXPECTED_LOCK:raise RuntimeError('FROZEN script SHA mismatch')
 run([sys.executable,BASE/'Supplementary_Code_6.py','--locked',lock,'--outdir',a.out/'common_driver'])
 run([sys.executable,BASE/'independent_audit_v18.py','--run-dir',a.out/'common_driver','--out',a.out/'independent_audit_v18.json'])
 if not a.skip_causal:
  run([sys.executable,BASE/'Supplementary_Code_7.py','--outdir',a.out/'causal_intervention'])
 if a.behavior_csv:
  run([sys.executable,BASE/'Supplementary_Code_16_source_reproducer.py','--csv',a.behavior_csv,'--output',a.out])
 tab=pd.read_csv(a.out/'common_driver'/'ExtendedData_Table8_common_driver_gates.csv')
 d=dict(zip(tab.condition,tab.all_three))
 assert bool(d['common_driver']) and not bool(d['no_change'])
 meta={'frozen_sha256':h,'python':platform.python_version(),'status':'reproduced synthetic counterexample; NO new biological acquisition',
       'positive_synthetic_case':bool(d['common_driver']),'negative_synthetic_case':bool(d['no_change'])}
 (a.out/'REPRODUCIBILITY_RECEIPT.json').write_text(json.dumps(meta,indent=2))
 print('PASS:',json.dumps(meta))
if __name__=='__main__':main()