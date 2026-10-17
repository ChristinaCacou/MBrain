#!/usr/bin/env python3
# Wrapper correcting row-index handling for the 28-dpi featureCounts matrix.
import importlib.util
from pathlib import Path
import pandas as pd

p=Path(__file__).with_name('long_covid_raw_transcriptomics.py')
spec=importlib.util.spec_from_file_location('base_analysis',p)
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

def read_featurecounts_fixed(path, sep=None):
    if sep is None:
        sep = ',' if '.csv' in path.name else '\t'
    df=pd.read_csv(path,sep=sep,compression='gzip',low_memory=False)
    annotation_names={'Chr','Start','End','Strand','Length'}
    # Some featureCounts exports have an unlabeled gene-ID field; pandas then
    # promotes that extra first field to the row index. Preserve it rather
    # than incorrectly treating the first BAM column as the gene identifier.
    if 'Geneid' in df.columns:
        ids=df['Geneid'].astype(str)
        candidate=[c for c in df.columns if c not in annotation_names|{'Geneid'}]
    elif not isinstance(df.index,pd.RangeIndex):
        ids=pd.Series(df.index.astype(str),index=df.index)
        candidate=[c for c in df.columns if c not in annotation_names]
    else:
        ids=df.iloc[:,0].astype(str)
        candidate=[c for c in df.columns[1:] if c not in annotation_names]
    sample_cols=[]
    for c in candidate:
        x=pd.to_numeric(df[c],errors='coerce')
        if x.notna().mean()>.95: sample_cols.append(c)
    if not sample_cols:
        raise RuntimeError(f'No numeric sample columns found in {path}; columns={list(df.columns)[:30]}')
    tmp=df[sample_cols].copy()
    tmp.insert(0,'__gene__',ids.to_numpy())
    tmp['__gene__']=tmp['__gene__'].astype(str).str.replace(r'\.\d+$','',regex=True)
    counts=tmp.groupby('__gene__',as_index=True)[sample_cols].sum()
    return counts.apply(pd.to_numeric,errors='coerce').fillna(0)

mod.read_featurecounts=read_featurecounts_fixed
mod.main()
