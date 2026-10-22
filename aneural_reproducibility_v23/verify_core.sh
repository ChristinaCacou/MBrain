#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export MPLBACKEND=Agg
mkdir -p generated
python code/Supplementary_Code_6.py --locked code/Prospective_Stentor_Validation_Analysis.py --outdir generated/common_driver
python code/independent_recheck_v23.py --audit-root generated/common_driver --out generated/INDEPENDENT_ARITHMETIC_AUDIT_v23.json
python code/negative_control_calibration_v22.py
python - <<'PYTEST'
import json,pathlib,pandas as pd
r=pathlib.Path('generated/common_driver/ExtendedData_Table8_common_driver_gates.csv')
x=pd.read_csv(r).set_index('condition')
assert bool(x.loc['common_driver','all_three']) and not bool(x.loc['no_change','all_three'])
assert abs(x.loc['common_driver','cross_modal_auc']-0.842122)<1e-5
a=json.loads(pathlib.Path('generated/INDEPENDENT_ARITHMETIC_AUDIT_v23.json').read_text())
assert a['computational'][0]['exact_4day_two_sided_signflip_p']==0.125
assert a['computational'][0]['max_abs_svd_eigen_difference']<1e-10
print('ALL REPRODUCIBILITY CHECKS PASSED; synthetics only')
PYTEST
