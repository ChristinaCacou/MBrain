#!/usr/bin/env bash
set -euo pipefail
mkdir -p analysis/alm3_inventory
aws s3 ls --no-sign-request --recursive s3://crcnsarchive/alm-3/ > analysis/alm3_inventory/s3_listing.txt
python - <<'PY'
from pathlib import Path
p=Path('analysis/alm3_inventory/s3_listing.txt')
lines=p.read_text().splitlines()
print('objects',len(lines))
for ln in lines:
    low=ln.lower()
    if any(k in low for k in ['data_structure','meta_data','analysis','demo','description','filelist']):
        print(ln)
PY
