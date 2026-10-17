#!/usr/bin/env python3
import re, tarfile, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

DATA=Path('analysis_data'); OUT=Path('analysis_outputs'); OUT.mkdir(exist_ok=True)

def read_fc(path,sep):
    df=pd.read_csv(path,sep=sep,compression='gzip',low_memory=False)
    ann={'Chr','Start','End','Strand','Length'}
    if 'Geneid' in df.columns:
        ids=df['Geneid'].astype(str); cand=[c for c in df.columns if c not in ann|{'Geneid'}]
    elif not isinstance(df.index,pd.RangeIndex):
        ids=pd.Series(df.index.astype(str),index=df.index); cand=[c for c in df.columns if c not in ann]
    else:
        ids=df.iloc[:,0].astype(str); cand=[c for c in df.columns[1:] if c not in ann]
    samples=[c for c in cand if pd.to_numeric(df[c],errors='coerce').notna().mean()>.95]
    x=df[samples].copy(); x.insert(0,'gene',ids.to_numpy()); x.gene=x.gene.astype(str).str.replace(r'\.\d+$','',regex=True)
    return x.groupby('gene')[samples].sum().astype(int)

def meta225(cols,time):
    # Mapping follows GEO sample order / library names verified in the prior raw-count run.
    if time==4:
        rows=[]
        for c in cols:
            m=re.search(r'MB-(11|15)_(Un|V1)_R(\d)',str(c));
            if not m: raise RuntimeError(f'Unrecognized 4dpi column {c}')
            rows.append({'sample':c,'line':'L'+m.group(1),'infected':'infected' if m.group(2)=='V1' else 'mock'})
        return pd.DataFrame(rows).set_index('sample')
    # 28dpi LCSB_1..16 order verified against GEO library-name metadata.
    rows=[]
    for c in cols:
        m=re.search(r'LCSB_(\d+)_',str(c));
        if not m: raise RuntimeError(f'Unrecognized 28dpi column {c}')
        n=int(m.group(1))
        if 1<=n<=4: line,inf='L11','mock'
        elif 5<=n<=8: line,inf='L15','mock'
        elif 9<=n<=12: line,inf='L11','infected'
        elif 13<=n<=16: line,inf='L15','infected'
        else: raise RuntimeError(c)
        rows.append({'sample':c,'line':line,'infected':inf})
    return pd.DataFrame(rows).set_index('sample')

def run_deseq(counts,meta,design,contrast):
    counts=counts.loc[counts.sum(axis=1)>=10].T.astype(int)
    meta=meta.loc[counts.index].copy()
    dds=DeseqDataSet(counts=counts,metadata=meta,design=design,refit_cooks=False,n_cpus=2,quiet=True)
    dds.deseq2()
    if contrast=='interaction':
        dm=dds.obsm['design_matrix']
        cols=list(dm.columns) if hasattr(dm,'columns') else list(dds.variables)
        print('Design columns:',cols)
        idx=[i for i,c in enumerate(cols) if 'virus' in str(c).lower() and 'pff' in str(c).lower()]
        if len(idx)!=1: raise RuntimeError(f'Cannot identify interaction coefficient from {cols}')
        vec=np.zeros(len(cols)); vec[idx[0]]=1.0
        ds=DeseqStats(dds,contrast=vec,cooks_filter=False,independent_filter=True,n_cpus=2,quiet=True)
    else:
        ds=DeseqStats(dds,contrast=contrast,cooks_filter=False,independent_filter=True,n_cpus=2,quiet=True)
    ds.summary()
    r=ds.results_df.copy(); r.index=r.index.astype(str)
    return r,dds

def parse211():
    tarpath=DATA/'GSE211934_RAW.tar'; edir=DATA/'g211_py'; edir.mkdir(exist_ok=True)
    with tarfile.open(tarpath,'r') as tf: tf.extractall(edir,filter='data')
    sample_map={'GSM6505758':('no','no'),'GSM6505759':('no','no'),'GSM6505760':('yes','no'),'GSM6505761':('yes','no'),'GSM6505762':('no','yes'),'GSM6505763':('no','yes'),'GSM6505764':('yes','yes'),'GSM6505765':('yes','yes')}
    ss={}
    for f in sorted(edir.rglob('*')):
        mm=re.search(r'(GSM\d+)',f.name)
        if not f.is_file() or not mm or mm.group(1) not in sample_map: continue
        gsm=mm.group(1); df=pd.read_csv(f,sep='\t',compression='infer',comment='#',low_memory=False)
        low={c:str(c).lower() for c in df.columns}
        idc=next((c for c in df.columns if low[c] in {'geneid','gene_id','gene','symbol','gene_symbol','gene name','gene_name','genesymbol'}),None)
        if idc is None: idc=next((c for c in df.columns if ('gene' in low[c] or 'symbol' in low[c]) and 'length' not in low[c] and 'count' not in low[c]),None)
        if idc is None: raise RuntimeError(f'ID col missing {df.columns}')
        countc=next((c for c in df.columns if ('count' in low[c] or low[c] in {'reads','numreads'}) and pd.to_numeric(df[c],errors='coerce').notna().mean()>.95),None)
        if countc is None:
            num=[c for c in df.columns if c!=idc and pd.to_numeric(df[c],errors='coerce').notna().mean()>.95 and low[c] not in {'start','end','length','chr','chrom','chromosome'}]
            countc=num[-1]
        s=pd.Series(pd.to_numeric(df[countc],errors='coerce').fillna(0).to_numpy(),index=df[idc].astype(str).str.replace(r'\.\d+$','',regex=True),name=gsm).groupby(level=0).sum(); ss[gsm]=s
    counts=pd.concat([ss[g] for g in sorted(ss)],axis=1).fillna(0).astype(int)
    meta=pd.DataFrame([{'sample':g,'virus':sample_map[g][0],'pff':sample_map[g][1]} for g in sorted(ss)]).set_index('sample')
    return counts,meta

def freedman_lane_module(counts,meta,genes,nperm=40320):
    # Module average log2CPM; exact residual permutation under additive null.
    idx={str(g).upper():g for g in counts.index}; found=[idx[g] for g in genes if g in idx]
    if len(found)<2: return {'n_genes':len(found),'genes_found':';'.join(map(str,found)),'interaction_beta':np.nan,'perm_p':np.nan}
    c=counts.loc[found,meta.index]
    lib=counts[meta.index].sum(axis=0).astype(float)
    y=np.log2((c.astype(float)+0.5).div(lib+1,axis=1)*1e6).mean(axis=0).to_numpy()
    v=(meta.virus=='yes').astype(float).to_numpy(); p=(meta.pff=='yes').astype(float).to_numpy()
    X0=np.column_stack([np.ones(8),v,p]); X1=np.column_stack([X0,v*p])
    b0=np.linalg.lstsq(X0,y,rcond=None)[0]; fitted=X0@b0; resid=y-fitted
    obs=np.linalg.lstsq(X1,y,rcond=None)[0][-1]
    more=0; total=0
    for perm in __import__('itertools').permutations(range(8)):
        ys=fitted+resid[list(perm)]
        bi=np.linalg.lstsq(X1,ys,rcond=None)[0][-1]
        more += abs(bi)>=abs(obs)-1e-12; total+=1
    return {'n_genes':len(found),'genes_found':';'.join(map(str,found)),'interaction_beta':float(obs),'perm_p':float((more+1)/(total+1)),'permutations':total}

def main():
    c4=read_fc(DATA/'GSE225517_4dpi.csv.gz',','); c28=read_fc(DATA/'GSE225517_28dpi.txt.gz','\t')
    m4=meta225(c4.columns,4); m28=meta225(c28.columns,28)
    r4,_=run_deseq(c4,m4,'~ line + infected',['infected','infected','mock'])
    r28,_=run_deseq(c28,m28,'~ line + infected',['infected','infected','mock'])
    r4.to_csv(OUT/'GSE225517_4dpi_PyDESeq2.csv'); r28.to_csv(OUT/'GSE225517_28dpi_PyDESeq2.csv')
    common=r4.index.intersection(r28.index); a=r4.loc[common]; b=r28.loc[common]
    sig4=set(a.index[a.padj.fillna(1)<.05]); sig28=set(b.index[b.padj.fillna(1)<.05]); ov=sorted(sig4&sig28)
    same=int((np.sign(a.loc[ov].log2FoldChange)==np.sign(b.loc[ov].log2FoldChange)).sum()) if ov else 0
    rho,p_rho=stats.spearmanr(a.log2FoldChange,b.log2FoldChange,nan_policy='omit')
    s225={'common_genes':len(common),'FDR_4dpi':len(sig4),'FDR_28dpi':len(sig28),'FDR_both':len(ov),'same_direction_both':same,'same_direction_fraction_both':same/len(ov) if ov else None,'overlap_hypergeom_p':float(stats.hypergeom.sf(len(ov)-1,len(common),len(sig4),len(sig28))) if ov else 1.0,'spearman_all_effects':float(rho),'spearman_all_effects_p':float(p_rho)}
    pd.DataFrame({'lfc4':a.log2FoldChange,'q4':a.padj,'lfc28':b.log2FoldChange,'q28':b.padj}).sort_values('q28').head(500).to_csv(OUT/'GSE225517_PyDESeq2_persistence_top500.csv')

    c211,m211=parse211(); r211,dds211=run_deseq(c211,m211,'~ virus + pff + virus:pff','interaction'); r211.to_csv(OUT/'GSE211934_interaction_PyDESeq2.csv')
    s211={'filtered_genes':len(r211),'interaction_p_lt_0_05':int((r211.pvalue.fillna(1)<.05).sum()),'interaction_FDR_lt_0_05':int((r211.padj.fillna(1)<.05).sum()),'interaction_abs_log2FC_ge_1':int((r211.log2FoldChange.abs()>=1).sum())}
    modules={'source_apoptosis_qpcr':['FADD','BAX','BID','CASP9'],'source_autophagy_qpcr':['WDR45','TECPR2','MAP1LC3A','SQSTM1'],'dopaminergic_presynaptic_hypothesis':['TH','DDC','SLC18A2','SLC6A3','SNAP25','SYN1','SYT1','VAMP2']}
    mod=[]
    for name,genes in modules.items():
        z=freedman_lane_module(c211,m211,genes); z['module']=name
        idx={str(g).upper():g for g in r211.index}; found=[idx[g] for g in genes if g in idx]
        if found:
            z['mean_PyDESeq2_interaction_log2FC']=float(r211.loc[found,'log2FoldChange'].mean())
            z['individual_interaction_pvalues']=';'.join(f'{g}:{r211.loc[g,"pvalue"]:.4g}' for g in found)
        mod.append(z)
    pd.DataFrame(mod).to_csv(OUT/'GSE211934_module_exact_permutation.csv',index=False)
    summary={'method':'PyDESeq2 0.5.x negative-binomial reanalysis plus exact Freedman-Lane residual permutation for prespecified modules','GSE225517':s225,'GSE211934':s211,'modules':mod}
    (OUT/'transcriptomic_pydeseq2_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
