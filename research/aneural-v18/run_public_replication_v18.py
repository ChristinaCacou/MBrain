#!/usr/bin/env python3
"""Reproduce the public v18 synthetic audit from source with pinned random seeds.
This is a post-registration exploratory computational analysis, not biology."""
from pathlib import Path
import subprocess, sys
BASE=Path(__file__).resolve().parent
SRC=BASE/'src'
def run(*args):
    print('RUN', *args, flush=True)
    subprocess.run([sys.executable, *[str(a) for a in args]], cwd=BASE, check=True)
run(SRC/'Supplementary_Code_6.py', '--locked', SRC/'Prospective_Stentor_Validation_Analysis.py', '--outdir', SRC/'frozen_audit')
run(SRC/'independent_verification.py')
run(SRC/'Supplementary_Code_7.py', '--outdir', BASE/'results'/'fresh_causal_intervention')
print('Done; results are SYNTHETIC, and four-day two-sided sign-flip P is 0.125.')
