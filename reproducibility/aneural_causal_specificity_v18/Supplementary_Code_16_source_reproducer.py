#!/usr/bin/env python3
"""Recompute Extended Data Figure 12 from original Doan et al. follow-up CSV.
This program fetches/reads real recorded observations; it never creates substitute data.
The manuscript's bundled aggregate tables were extracted from the pinned GitHub blob.
"""
import argparse, csv, hashlib, io, urllib.request
from pathlib import Path
import pandas as pd
SOURCE='https://raw.githubusercontent.com/nhithuydoan/associative_learning_Stentor/main/data/dataset_followup_2.csv'
BLOB='04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc'
p=argparse.ArgumentParser();p.add_argument('--csv',type=Path);p.add_argument('--output',type=Path,default=Path('.'));a=p.parse_args()
raw=a.csv.read_bytes() if a.csv else urllib.request.urlopen(SOURCE,timeout=120).read()
git_blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
if git_blob != BLOB: raise SystemExit('REJECT: source Git blob identity changed; examine before reanalysis: '+git_blob)
df=pd.read_csv(io.BytesIO(raw));needed={'condition','folder','phase','run','stimulus','cell_number','contract'}
if not needed.issubset(df.columns):raise SystemExit('REJECT: source columns incomplete')
if len(df)!=13546:raise SystemExit('REJECT: unexpected count of source observations')
df['run_id']=df['folder'].astype(str)+'|'+df['run'].astype(str)
df['recording_date']=df['folder'].astype(str).str[:10]
test=df[(df.phase=='test') & df.contract.notna()]
if test[['run_id','cell_number']].duplicated().any():raise SystemExit('REJECT: duplicate cell test')
assert len(test)==521 and test.run_id.nunique()==19
pre=df[(df.phase=='hab') & df.contract.notna()]
assert pre.recording_date.nunique()==8
first=pre[pre.stimulus<=5];last=pre[pre.stimulus>=16]
def summary(x):return x.groupby('run_id')['contract'].agg(['sum','count'])
r=summary(first).join(summary(last),lsuffix='_early',rsuffix='_late').reset_index()
r['early_rate']=r.sum_early/r.count_early;r['late_rate']=r.sum_late/r.count_late
assert (r.late_rate<r.early_rate).all()
r.to_csv(a.output/'reproduced_run_summary.csv',index=False)
print('Verified source blob:',git_blob,'; first5 rate',first.contract.mean(),'; last5 rate',last.contract.mean())
print('All',len(r),'runs decrease; number of recording days:',pre.recording_date.nunique())