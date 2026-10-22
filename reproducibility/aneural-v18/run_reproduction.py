#!/usr/bin/env python3
"""Regenerate main synthetic numerical tables and run independent numerical checks.
This script does not generate biological measurements."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parent
PY=sys.executable
cmds=[
 [PY,str(ROOT/'code/adversarial_common_driver.py'),'--locked',str(ROOT/'code/frozen_protocol_unchanged.py'),'--outdir',str(ROOT/'central_reproduction')],
 [PY,str(ROOT/'code/causal_interventions.py'),'--outdir',str(ROOT/'causal_reproduction')]
]
for cmd in cmds:
 print('RUNNING:', ' '.join(cmd),flush=True)
 subprocess.run(cmd,check=True)
subprocess.run([PY,str(ROOT/'code/independent_numeric_audit.py')],check=True)
