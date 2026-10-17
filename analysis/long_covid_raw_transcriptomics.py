#!/usr/bin/env python3
import os, re, json, tarfile, gzip, urllib.request, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path('analysis_outputs')
OUT.mkdir(exist_ok=True)
DATA = Path('analysis_data')
DATA.mkdir(exist_ok=True)

URLS = {
 'g225_4': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225517/suppl/GSE225517_feature_counts_covid19_data.csv.gz',
 'g225_28': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225517/suppl/GSE225517_feature_counts_covid19_midbrain_2024.txt.gz',
 'g211_tar': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211934/suppl/GSE211934_RAW.tar',
}

def download(url, path):
    if not path.exists():
        print('Downloading', url)
        urllib.request.urlretrieve(url, path)
    print(path, path.stat().st_size, 'bytes')

def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / np.arange(1,n+1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked,1)
    return q

def read_featurecounts(path, sep=None):
    if sep is None:
        sep = ',' if '.csv' in path.name else '\t'
    df = pd.read_csv(path, sep=sep, compression='gzip', low_memory=False)
    # Prefer Geneid; otherwise first column.
    idcol = 'Geneid' if 'Geneid' in df.columns else df.columns[0]
    ann = {idcol,'Chr','Start','End','Strand','Length'}
    sample_cols=[]
    for c in df.columns:
        if c in ann: continue
        x=pd.to_numeric(df[c],errors='coerce')
        if x.notna().mean() > .95:
            sample_cols.append(c)
    if not sample_cols:
        raise RuntimeError(f'No numeric sample columns found in {path}; columns={list(df.columns)[:30]}')
    counts=df[[idcol]+sample_cols].copy()
    counts[idcol]=counts[idcol].astype(str).str.replace(r'\.\d+$','',regex=True)
    counts=counts.groupby(idcol,as_index=True)[sample_cols].sum()
    counts=counts.apply(pd.to_numeric,errors='coerce').fillna(0)
    return counts

def logcpm(counts):
    lib=counts.sum(axis=0).astype(float)
    return np.log2((counts.astype(float)+0.5).div(lib+1,axis=1)*1e6)

def ols_same_design(Y, X, coef_idx):
    # Y genes x samples, X samples x p
    XtX_inv=np.linalg.inv(X.T@X)
    B=XtX_inv@X.T@Y.T
    fitted=X@B
    resid=Y.T-fitted
    df=X.shape[0]-X.shape[1]
    s2=(resid**2).sum(axis=0)/df
    se=np.sqrt(np.maximum(s2,0)*XtX_inv[coef_idx,coef_idx])
    beta=B[coef_idx]
    t=np.divide(beta,se,out=np.zeros_like(beta),where=se>0)
    p=2*stats.t.sf(np.abs(t),df)
    return beta,se,t,p,df

def geo_sample_text(gsm):
    url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
    with urllib.request.urlopen(url,timeout=60) as r:
        return r.read().decode('utf-8','replace')

def map_225_columns(count_cols, timepoint):
    meta=[]
    if timepoint==4:
        gsms=[f'GSM70497{i}' for i in range(53,69)]
    else:
        gsms=[f'GSM91393{i}' for i in range(62,78)]
    for gsm in gsms:
        txt=geo_sample_text(gsm)
        title=re.search(r'!Sample_title\s*=\s*(.+)',txt)
        title=title.group(1).strip() if title else ''
        line=re.search(r'!Sample_characteristics_ch1\s*=\s*cell line:\s*(\d+)',txt,re.I)
        if not line:
            line=re.search(r'Line\s+(11|15)',title,re.I)
        line=int(line.group(1)) if line else None
        treatment='infected' if ('SARS-CoV-2 infection' in title or 'V1' in title) else 'mock'
        desc=re.findall(r'!Sample_description\s*=\s*(.+)',txt)
        tokens=[gsm]
        m=re.search(r'\[([^\]]+)\]',title)
        if m: tokens.append(m.group(1))
        for d in desc:
            lm=re.search(r'Library name:\s*([^\s;]+)',d)
            if lm: tokens.append(lm.group(1))
        meta.append({'gsm':gsm,'title':title,'line':line,'treatment':treatment,'tokens':tokens})
    mapping={}
    for c in count_cols:
        cs=str(c)
        hits=[]
        for m in meta:
            for tok in m['tokens']:
                if tok and (tok in cs or tok in os.path.basename(cs) or cs in tok):
                    hits.append(m); break
        if len(hits)==1:
            mapping[c]=hits[0]
    if len(mapping) != len(count_cols):
        diag={'timepoint':timepoint,'columns':list(map(str,count_cols)),'mapped':{str(k):v['gsm'] for k,v in mapping.items()},'metadata':meta}
        (OUT/f'gse225517_mapping_diagnostic_{timepoint}dpi.json').write_text(json.dumps(diag,indent=2))
        raise RuntimeError(f'Could map only {len(mapping)}/{len(count_cols)} count columns at {timepoint} dpi. Diagnostic written.')
    return pd.DataFrame([{'column':c,'gsm':mapping[c]['gsm'],'line':mapping[c]['line'],'infected':1 if mapping[c]['treatment']=='infected' else 0} for c in count_cols])

def analyze_225(path4,path28):
    c4=read_featurecounts(path4,',')
    c28=read_featurecounts(path28,'\t')
    print('GSE225517 4dpi shape',c4.shape,'cols',list(c4.columns))
    print('GSE225517 28dpi shape',c28.shape,'cols',list(c28.columns))
    m4=map_225_columns(list(c4.columns),4)
    m28=map_225_columns(list(c28.columns),28)
    # Reorder by metadata rows; filter low-expression within each timepoint as source did.
    c4=c4[m4['column'].tolist()]; c28=c28[m28['column'].tolist()]
    c4=c4.loc[c4.sum(axis=1)>=10]; c28=c28.loc[c28.sum(axis=1)>=10]
    y4=logcpm(c4); y28=logcpm(c28)
    def fit(y,m):
        line15=(m['line'].to_numpy()==15).astype(float)
        inf=m['infected'].to_numpy().astype(float)
        X=np.column_stack([np.ones(len(m)),line15,inf])
        b,se,t,p,df=ols_same_design(y.to_numpy(),X,2)
        q=bh(p)
        return pd.DataFrame({'gene':y.index,'log2cpm_infection_effect':b,'se':se,'t':t,'p':p,'q':q,'df':df}).set_index('gene')
    r4=fit(y4,m4); r28=fit(y28,m28)
    common=r4.index.intersection(r28.index)
    a=r4.loc[common]; b=r28.loc[common]
    sig4=set(a.index[a.q<.05]); sig28=set(b.index[b.q<.05]); ov=sorted(sig4&sig28)
    M=len(common); n4=len(sig4); n28=len(sig28); k=len(ov)
    hyper=float(stats.hypergeom.sf(k-1,M,n4,n28)) if k else 1.0
    same=0
    if k:
        same=int((np.sign(a.loc[ov,'log2cpm_infection_effect'])==np.sign(b.loc[ov,'log2cpm_infection_effect'])).sum())
    sign_p=float(stats.binomtest(same,k,.5).pvalue) if k else 1.0
    rho_all,p_rho_all=stats.spearmanr(a.log2cpm_infection_effect,b.log2cpm_infection_effect)
    if k>=3:
        rho_ov,p_rho_ov=stats.spearmanr(a.loc[ov,'log2cpm_infection_effect'],b.loc[ov,'log2cpm_infection_effect'])
    else: rho_ov=p_rho_ov=np.nan
    union=len(sig4|sig28); jacc=k/union if union else np.nan
    # Threshold-free directional concordance of the strongest decile at either time point.
    rank_strength=np.maximum(np.abs(a.log2cpm_infection_effect.to_numpy()),np.abs(b.log2cpm_infection_effect.to_numpy()))
    cutoff=np.quantile(rank_strength,.90)
    top=common[rank_strength>=cutoff]
    top_same=int((np.sign(a.loc[top,'log2cpm_infection_effect'])==np.sign(b.loc[top,'log2cpm_infection_effect'])).sum())
    top_p=float(stats.binomtest(top_same,len(top),.5).pvalue)
    summary={
      'dataset':'GSE225517','common_filtered_genes':M,'sig_FDR_4dpi':n4,'sig_FDR_28dpi':n28,
      'sig_both':k,'overlap_hypergeom_p':hyper,'jaccard_sig':jacc,'same_direction_both':same,
      'same_direction_fraction_both':same/k if k else np.nan,'sign_concordance_p_both':sign_p,
      'spearman_effect_all':float(rho_all),'spearman_effect_all_p':float(p_rho_all),
      'spearman_effect_overlap':float(rho_ov) if np.isfinite(rho_ov) else np.nan,'spearman_effect_overlap_p':float(p_rho_ov) if np.isfinite(p_rho_ov) else np.nan,
      'top_decile_n':len(top),'top_decile_same_direction':top_same,'top_decile_same_direction_fraction':top_same/len(top),'top_decile_sign_p':top_p,
      'analysis':'raw counts; within-time log2CPM OLS adjusted for iPSC line; cross-time persistence uses within-time effects only'
    }
    both=pd.DataFrame({'effect_4dpi':a.log2cpm_infection_effect,'q_4dpi':a.q,'effect_28dpi':b.log2cpm_infection_effect,'q_28dpi':b.q})
    both['persistent_fdr']=(both.q_4dpi<.05)&(both.q_28dpi<.05)
    both['same_direction']=np.sign(both.effect_4dpi)==np.sign(both.effect_28dpi)
    both['min_abs_effect']=np.minimum(np.abs(both.effect_4dpi),np.abs(both.effect_28dpi))
    both.sort_values(['persistent_fdr','min_abs_effect'],ascending=[False,False]).head(250).to_csv(OUT/'GSE225517_persistence_top250.csv')
    r4.sort_values('q').head(250).to_csv(OUT/'GSE225517_4dpi_top250.csv')
    r28.sort_values('q').head(250).to_csv(OUT/'GSE225517_28dpi_top250.csv')
    m4.to_csv(OUT/'GSE225517_4dpi_sample_mapping.csv',index=False); m28.to_csv(OUT/'GSE225517_28dpi_sample_mapping.csv',index=False)
    return summary,both

def read_g211_tar(tarpath):
    edir=DATA/'g211'; edir.mkdir(exist_ok=True)
    with tarfile.open(tarpath,'r') as tf:
        tf.extractall(edir,filter='data')
    sample_map={
      'GSM6505758':(0,0),'GSM6505759':(0,0),
      'GSM6505760':(1,0),'GSM6505761':(1,0),
      'GSM6505762':(0,1),'GSM6505763':(0,1),
      'GSM6505764':(1,1),'GSM6505765':(1,1),
    }
    ss={}
    for f in sorted(edir.rglob('*')):
        if not f.is_file(): continue
        m=re.search(r'(GSM\d+)',f.name)
        if not m or m.group(1) not in sample_map: continue
        gsm=m.group(1)
        try:
            df=pd.read_csv(f,sep='\t',compression='infer',comment='#',low_memory=False)
        except Exception:
            df=pd.read_csv(f,sep=None,engine='python',compression='infer',comment='#',low_memory=False)
        num=[c for c in df.columns if pd.to_numeric(df[c],errors='coerce').notna().mean()>.95]
        if not num: raise RuntimeError(f'No count column in {f}; cols={list(df.columns)}')
        idc=df.columns[0]; vc=num[-1]
        s=pd.Series(pd.to_numeric(df[vc],errors='coerce').fillna(0).to_numpy(),index=df[idc].astype(str).str.replace(r'\.\d+$','',regex=True),name=gsm)
        s=s.groupby(level=0).sum(); ss[gsm]=s
    if set(ss)!=set(sample_map): raise RuntimeError(f'GSE211934 samples found {sorted(ss)}, expected {sorted(sample_map)}')
    counts=pd.concat([ss[g] for g in sorted(ss)],axis=1).fillna(0)
    meta=pd.DataFrame([{'gsm':g,'virus':sample_map[g][0],'pff':sample_map[g][1]} for g in sorted(ss)])
    counts=counts[meta.gsm.tolist()]
    return counts,meta

def analyze_211(tarpath):
    counts,meta=read_g211_tar(tarpath)
    counts=counts.loc[counts.sum(axis=1)>=10]
    y=logcpm(counts)
    v=meta.virus.to_numpy(float); pff=meta.pff.to_numpy(float)
    X=np.column_stack([np.ones(len(meta)),v,pff,v*pff])
    b,se,t,p,df=ols_same_design(y.to_numpy(),X,3); q=bh(p)
    r=pd.DataFrame({'gene':y.index,'interaction_log2cpm':b,'se':se,'t':t,'p':p,'q':q,'df':df}).set_index('gene')
    # Cell means and raw factorial contrasts for transparency.
    groups={
      'control':(v==0)&(pff==0),'virus':(v==1)&(pff==0),'pff':(v==0)&(pff==1),'pff_virus':(v==1)&(pff==1)
    }
    for name,mask in groups.items(): r[f'mean_{name}']=y.loc[:,mask].mean(axis=1)
    r['pff_effect_no_virus']=r.mean_pff-r.mean_control
    r['pff_effect_with_virus']=r.mean_pff_virus-r.mean_virus
    r['virus_effect_no_pff']=r.mean_virus-r.mean_control
    r['virus_effect_with_pff']=r.mean_pff_virus-r.mean_pff
    # Pre-specified source-qPCR modules + hypothesis module.
    modules={
      'source_apoptosis_qpcr':['FADD','BAX','BID','CASP9'],
      'source_autophagy_qpcr':['WDR45','TECPR2','MAP1LC3A','SQSTM1'],
      'dopaminergic_presynaptic_hypothesis':['TH','DDC','SLC18A2','SLC6A3','SNAP25','SYN1','SYT1','VAMP2']
    }
    idx_upper={str(g).upper():g for g in y.index}
    modrows=[]
    for mn,genes in modules.items():
        found=[idx_upper[g] for g in genes if g in idx_upper]
        if len(found)<2:
            modrows.append({'module':mn,'n_genes':len(found),'genes_found':';'.join(map(str,found)),'interaction_beta':np.nan,'p':np.nan})
            continue
        score=y.loc[found].mean(axis=0).to_numpy()[None,:]
        mb,mse,mt,mp,mdf=ols_same_design(score,X,3)
        modrows.append({'module':mn,'n_genes':len(found),'genes_found':';'.join(map(str,found)),'interaction_beta':float(mb[0]),'se':float(mse[0]),'t':float(mt[0]),'p':float(mp[0]),'df':mdf})
    mod=pd.DataFrame(modrows)
    mod.to_csv(OUT/'GSE211934_module_interactions.csv',index=False)
    r.sort_values('p').head(500).to_csv(OUT/'GSE211934_interaction_top500.csv')
    meta.to_csv(OUT/'GSE211934_sample_mapping.csv',index=False)
    # Raw effect-size amplification diagnostics, without claiming source DESeq2 significance.
    pff0=r.pff_effect_no_virus
    v_in_pff=r.virus_effect_with_pff
    amp_up=((pff0>0)&(v_in_pff>0))
    amp_down=((pff0<0)&(v_in_pff<0))
    strong_up=amp_up&(pff0.abs()>=.5)&(v_in_pff.abs()>=.5)
    strong_down=amp_down&(pff0.abs()>=.5)&(v_in_pff.abs()>=.5)
    summary={
      'dataset':'GSE211934','filtered_genes':len(r),'gene_interaction_q_lt_0_05':int((r.q<.05).sum()),
      'gene_interaction_p_lt_0_05':int((r.p<.05).sum()),'interaction_abs_ge_1':int((r.interaction_log2cpm.abs()>=1).sum()),
      'directional_amplification_up_all':int(amp_up.sum()),'directional_amplification_down_all':int(amp_down.sum()),
      'directional_amplification_up_both_effects_ge_0_5':int(strong_up.sum()),'directional_amplification_down_both_effects_ge_0_5':int(strong_down.sum()),
      'analysis':'raw counts; log2CPM factorial OLS virus + PFF + virus:PFF; n=2/cell, so gene-wise p/q are exploratory and module/effect-size interpretation is primary'
    }
    return summary,r,mod

def main():
    p4=DATA/'GSE225517_4dpi.csv.gz'; p28=DATA/'GSE225517_28dpi.txt.gz'; ptar=DATA/'GSE211934_RAW.tar'
    download(URLS['g225_4'],p4); download(URLS['g225_28'],p28); download(URLS['g211_tar'],ptar)
    s225,_=analyze_225(p4,p28)
    s211,_,mod=analyze_211(ptar)
    summary={'generated_utc':pd.Timestamp.utcnow().isoformat(),'GSE225517':s225,'GSE211934':s211,'GSE211934_modules':mod.replace({np.nan:None}).to_dict(orient='records')}
    (OUT/'transcriptomic_raw_analysis_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    rows=[]
    for ds,d in [('GSE225517',s225),('GSE211934',s211)]:
        for k,v in d.items(): rows.append({'dataset':ds,'metric':k,'value':v})
    for _,rr in mod.iterrows(): rows.append({'dataset':'GSE211934','metric':'module_'+str(rr['module'])+'_interaction_beta','value':rr.get('interaction_beta')}); rows.append({'dataset':'GSE211934','metric':'module_'+str(rr['module'])+'_p','value':rr.get('p')})
    pd.DataFrame(rows).to_csv(OUT/'transcriptomic_raw_analysis_summary.csv',index=False)
    print(json.dumps(summary,indent=2,default=str))

if __name__=='__main__': main()
