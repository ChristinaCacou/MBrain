#!/usr/bin/env python3
# Final wrapper: preserves corrected GSE225517 row-index parsing and fixes
# identifier/count inference for per-sample GSE211934 raw-count tables.
import importlib.util, re, tarfile
from pathlib import Path
import pandas as pd

p=Path(__file__).with_name('long_covid_raw_transcriptomics.py')
spec=importlib.util.spec_from_file_location('base_analysis',p)
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

# Correct featureCounts parser used for GSE225517.
def read_featurecounts_fixed(path, sep=None):
    if sep is None: sep=',' if '.csv' in path.name else '\t'
    df=pd.read_csv(path,sep=sep,compression='gzip',low_memory=False)
    ann={'Chr','Start','End','Strand','Length'}
    if 'Geneid' in df.columns:
        ids=df['Geneid'].astype(str); candidate=[c for c in df.columns if c not in ann|{'Geneid'}]
    elif not isinstance(df.index,pd.RangeIndex):
        ids=pd.Series(df.index.astype(str),index=df.index); candidate=[c for c in df.columns if c not in ann]
    else:
        ids=df.iloc[:,0].astype(str); candidate=[c for c in df.columns[1:] if c not in ann]
    sample=[]
    for c in candidate:
        if pd.to_numeric(df[c],errors='coerce').notna().mean()>.95: sample.append(c)
    tmp=df[sample].copy(); tmp.insert(0,'__gene__',ids.to_numpy())
    tmp['__gene__']=tmp['__gene__'].astype(str).str.replace(r'\.\d+$','',regex=True)
    return tmp.groupby('__gene__',as_index=True)[sample].sum().apply(pd.to_numeric,errors='coerce').fillna(0)
mod.read_featurecounts=read_featurecounts_fixed

def read_g211_tar_fixed(tarpath):
    edir=mod.DATA/'g211'; edir.mkdir(exist_ok=True)
    with tarfile.open(tarpath,'r') as tf: tf.extractall(edir,filter='data')
    sample_map={'GSM6505758':(0,0),'GSM6505759':(0,0),'GSM6505760':(1,0),'GSM6505761':(1,0),'GSM6505762':(0,1),'GSM6505763':(0,1),'GSM6505764':(1,1),'GSM6505765':(1,1)}
    ss={}; diag=[]
    for f in sorted(edir.rglob('*')):
        if not f.is_file(): continue
        mm=re.search(r'(GSM\d+)',f.name)
        if not mm or mm.group(1) not in sample_map: continue
        gsm=mm.group(1)
        try: df=pd.read_csv(f,sep='\t',compression='infer',comment='#',low_memory=False)
        except Exception: df=pd.read_csv(f,sep=None,engine='python',compression='infer',comment='#',low_memory=False)
        cols=list(df.columns)
        if not diag: print('GSE211934 raw table columns:',cols); print(df.head(3).to_string())
        lower={c:str(c).strip().lower() for c in cols}
        # Prefer explicit gene identifier/symbol fields. Avoid genomic coordinates and chromosome labels.
        exact=['geneid','gene_id','gene id','gene','symbol','gene_symbol','gene symbol','gene_name','gene name','genesymbol','entrezgene','entrez_gene_id','ensembl_gene_id']
        idc=None
        for key in exact:
            for c in cols:
                if lower[c]==key: idc=c; break
            if idc is not None: break
        if idc is None:
            cand=[c for c in cols if ('gene' in lower[c] or 'symbol' in lower[c]) and not any(x in lower[c] for x in ['length','count'])]
            if cand: idc=cand[0]
        if idc is None:
            # Last resort: unique non-coordinate textual column.
            bad={'chr','chrom','chromosome','strand'}
            cand=[]
            for c in cols:
                if lower[c] in bad or lower[c] in {'start','end','length'}: continue
                numeric=pd.to_numeric(df[c],errors='coerce').notna().mean()
                if numeric<.5 and df[c].nunique()/max(len(df),1)>.5: cand.append(c)
            if cand: idc=cand[0]
        if idc is None: raise RuntimeError(f'Cannot infer gene ID column in {f}; columns={cols}')
        # Prefer explicit raw-count field.
        countc=None
        for c in cols:
            lc=lower[c]
            if lc in {'count','counts','raw_count','raw counts','read_count','read counts','numreads','reads'} or ('count' in lc and 'tpm' not in lc and 'fpkm' not in lc):
                if pd.to_numeric(df[c],errors='coerce').notna().mean()>.95: countc=c; break
        if countc is None:
            numeric=[c for c in cols if c!=idc and pd.to_numeric(df[c],errors='coerce').notna().mean()>.95]
            # Coordinates are not expression; remove them before falling back to the last numeric field.
            numeric=[c for c in numeric if lower[c] not in {'start','end','length','chr','chrom','chromosome'}]
            if not numeric: raise RuntimeError(f'Cannot infer count column in {f}; columns={cols}')
            countc=numeric[-1]
        diag.append({'gsm':gsm,'file':f.name,'id_column':str(idc),'count_column':str(countc),'rows':len(df)})
        ids=df[idc].astype(str).str.replace(r'\.\d+$','',regex=True)
        vals=pd.to_numeric(df[countc],errors='coerce').fillna(0)
        s=pd.Series(vals.to_numpy(),index=ids,name=gsm).groupby(level=0).sum(); ss[gsm]=s
    pd.DataFrame(diag).to_csv(mod.OUT/'GSE211934_parse_diagnostic.csv',index=False)
    if set(ss)!=set(sample_map): raise RuntimeError(f'GSE211934 samples found {sorted(ss)}, expected {sorted(sample_map)}')
    counts=pd.concat([ss[g] for g in sorted(ss)],axis=1).fillna(0)
    meta=pd.DataFrame([{'gsm':g,'virus':sample_map[g][0],'pff':sample_map[g][1]} for g in sorted(ss)])
    counts=counts[meta.gsm.tolist()]
    print('GSE211934 combined count matrix',counts.shape,'first IDs',list(counts.index[:10]))
    return counts,meta
mod.read_g211_tar=read_g211_tar_fixed
mod.main()
