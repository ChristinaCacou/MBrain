#!/usr/bin/env python3
"""iPSC-line robustness for GSE225517 raw-count infection effects."""
import urllib.request, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import importlib.util

# Import tested PyDESeq2 helpers without invoking its main().
p=Path(__file__).with_name('long_covid_pydeseq2.py')
spec=importlib.util.spec_from_file_location('de','analysis/long_covid_pydeseq2.py')
de=importlib.util.module_from_spec(spec); spec.loader.exec_module(de)
DATA=Path('analysis_data'); OUT=Path('analysis_outputs'); DATA.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
URLS={
'4':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225517/suppl/GSE225517_feature_counts_covid19_data.csv.gz',
'28':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225517/suppl/GSE225517_feature_counts_covid19_midbrain_2024.txt.gz'}
for t,url in URLS.items():
    path=DATA/('GSE225517_4dpi.csv.gz' if t=='4' else 'GSE225517_28dpi.txt.gz')
    if not path.exists(): urllib.request.urlretrieve(url,path)

def corr_metrics(r1,r2):
    common=r1.index.intersection(r2.index)
    a=r1.loc[common]; b=r2.loc[common]
    ok=a.log2FoldChange.notna()&b.log2FoldChange.notna()
    rho,p=stats.spearmanr(a.loc[ok,'log2FoldChange'],b.loc[ok,'log2FoldChange'])
    sig1=set(a.index[a.padj.fillna(1)<.05]); sig2=set(b.index[b.padj.fillna(1)<.05]); ov=sorted(sig1&sig2)
    same=sum(np.sign(a.loc[g,'log2FoldChange'])==np.sign(b.loc[g,'log2FoldChange']) for g in ov)
    return {'common':int(ok.sum()),'spearman':float(rho),'spearman_p':float(p),'fdr_a':len(sig1),'fdr_b':len(sig2),'fdr_overlap':len(ov),'same_direction_overlap':int(same),'same_direction_fraction_overlap':same/len(ov) if ov else None,'overlap_hypergeom_p':float(stats.hypergeom.sf(len(ov)-1,len(common),len(sig1),len(sig2))) if ov else 1.0}

def interaction_fit(counts,meta):
    # numeric contrast auto-identification from design matrix, mirroring main interaction analysis.
    counts=counts.loc[counts.sum(axis=1)>=10].T.astype(int); meta=meta.loc[counts.index].copy()
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    dds=DeseqDataSet(counts=counts,metadata=meta,design='~ line + infected + line:infected',refit_cooks=False,n_cpus=2,quiet=True); dds.deseq2()
    cols=list(dds.obsm['design_matrix'].columns)
    idx=[i for i,c in enumerate(cols) if 'line' in str(c).lower() and 'infected' in str(c).lower()]
    if len(idx)!=1: raise RuntimeError(f'Interaction column ambiguous: {cols}')
    vec=np.zeros(len(cols)); vec[idx[0]]=1
    ds=DeseqStats(dds,contrast=vec,cooks_filter=False,independent_filter=True,n_cpus=2,quiet=True); ds.summary(); r=ds.results_df.copy()
    return {'interaction_fdr_lt_0_05':int((r.padj.fillna(1)<.05).sum()),'interaction_p_lt_0_05':int((r.pvalue.fillna(1)<.05).sum()),'median_abs_interaction_lfc':float(r.log2FoldChange.abs().median()),'design_columns':';'.join(map(str,cols))}

def analyze_time(t,path,sep):
    c=de.read_fc(path,sep); m=de.meta225(c.columns,int(t))
    line_results={}
    for line in ['L11','L15']:
        mm=m[m.line==line]
        rr,_=de.run_deseq(c[mm.index],mm,'~ infected',['infected','infected','mock'])
        rr.to_csv(OUT/f'GSE225517_{t}dpi_{line}_PyDESeq2.csv'); line_results[line]=rr
    x=corr_metrics(line_results['L11'],line_results['L15'])
    x.update(interaction_fit(c,m)); return x,line_results

r4,lr4=analyze_time('4',DATA/'GSE225517_4dpi.csv.gz',',')
r28,lr28=analyze_time('28',DATA/'GSE225517_28dpi.txt.gz','\t')
# Within-line cross-time concordance.
cross={line:corr_metrics(lr4[line],lr28[line]) for line in ['L11','L15']}
summary={'within_time_line_concordance':{'4dpi':r4,'28dpi':r28},'within_line_cross_time_concordance':cross,'interpretation_note':'Line-specific fits test biological-background robustness within each sequencing wave. Cross-time correlations remain cross-wave comparisons and cannot separate temporal biology from batch/acquisition differences.'}
(OUT/'GSE225517_line_sensitivity_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
