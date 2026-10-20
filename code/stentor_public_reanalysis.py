#!/usr/bin/env python3
"""Independent secondary reanalysis of publicly processed Stentor learning RNA-seq.

Inputs from https://github.com/aralbright/2025_stentor_learning , pinned by Git blob SHA.
This is NOT processing the raw FASTQ reads and does not execute the prospective
Stentor protocol. Analyses here are explicitly exploratory and small-sample.

Run: python stentor_public_reanalysis.py --out reanalysis_outputs
Use --input-dir DIRECTORY for manually downloaded source CSV files.
"""
import argparse, hashlib, io, json, pathlib, urllib.request
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

BASE='https://raw.githubusercontent.com/aralbright/2025_stentor_learning/main/'
FILES={
 'mono4':('exp1/exp1_sleuth_4_wald_v2.csv','1413206cb156839e6450898258b352ee73121f94'),
 'mono19':('exp1/exp1_sleuth_19_wald_v2.csv','db78b7206d9e77adea197cef61a43bb8817077a1'),
 'poly19':('exp2/exp2_rna_t19ut19_sleuth_wald_v2.csv','1a741431cc96b2c9ab5ca23686610759ceb71681'),
 'recovery':('exp2/exp2_rna_t19f90_sleuth_wald.csv','9200d6749ae0bee8b810536a3b0511915f082610'),
 'norm':('exp2/exp2_rna_sleuth_obsnorm.csv','48e76227459806a280b72675c00685e19700a86e'),
}
GROUPS={
 'T':['DR4876_01','DR4876_04','DR4876_07'],
 'U':['DR4876_03','DR4876_06','DR4876_09'],
 'F10':['DR4876_10','DR4876_12','DR4876_14'],
 'F90':['DR4876_16','DR4876_18','DR4876_20'],
 'T0':['DR4876_02','DR4876_05','DR4876_08'],
}

def load_pinned(label,input_dir):
    path,sha=FILES[label]
    local=pathlib.Path(input_dir)/pathlib.Path(path).name if input_dir else None
    if local and local.exists(): raw=local.read_bytes()
    else:
        print('Retrieving:',BASE+path,flush=True)
        with urllib.request.urlopen(BASE+path,timeout=120) as h: raw=h.read()
    digest=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
    if digest!=sha: raise ValueError(f'{label}: source SHA mismatch: {digest} != {sha}')
    return pd.read_csv(io.BytesIO(raw),index_col=0)

def compare(a,b,which_a,which_b):
    ids=a.index.intersection(b.index)
    x=a.loc[ids,'b'].astype(float).to_numpy(); y=b.loc[ids,'b'].astype(float).to_numpy()
    ap=a.loc[ids,'qval'].astype(float).to_numpy()<0.05
    bp=b.loc[ids,'qval'].astype(float).to_numpy()<0.05
    overlap=ap&bp
    return dict(contrasts=[which_a,which_b],n_common=len(ids),
                pearson_effect_correlation=float(np.corrcoef(x,y)[0,1]),
                significant_a=int(ap.sum()),significant_b=int(bp.sum()),
                significant_overlap=int(overlap.sum()),
                overlap_same_direction=int(np.sum(overlap&(x*y>0))),
                overlap_expected=float(ap.sum()*bp.sum()/len(ids)),
                exploratory_hypergeom_p=float(hypergeom.sf(overlap.sum()-1,len(ids),int(ap.sum()),int(bp.sum()))))

def fold_select_counts(log_expr,k,groups=GROUPS,flips=(0,0,0)):
    """Within 3 matched sample blocks, train on 2 pairs, test 1 pair;
    all gene selection is repeated in each training fold. Assumes filtered input.
    """
    margins=[]
    for held in range(3):
        tr=[i for i in range(3) if i!=held]
        diff=np.mean(np.column_stack([
           log_expr[groups['T'][i]].to_numpy()-log_expr[groups['U'][i]].to_numpy()
           if flips[i]==0 else
           log_expr[groups['U'][i]].to_numpy()-log_expr[groups['T'][i]].to_numpy()
           for i in tr]),axis=1)
        pick=np.argsort(-np.abs(diff),kind='mergesort')[:k]
        delta=(log_expr[groups['T'][held]].to_numpy()-log_expr[groups['U'][held]].to_numpy())
        if flips[held]:delta=-delta
        margins.append(float(np.mean(np.sign(diff[pick])*delta[pick])))
    return margins

def run(input_dir,outdir):
    d={key:load_pinned(key,input_dir) for key in FILES}
    for key in ('mono4','mono19','poly19','recovery'):
        t=d[key]
        if 'target_id' not in t.columns:raise ValueError(f'{key}: target_id column missing')
        d[key]=t.set_index('target_id')[['pval','qval','b']].apply(pd.to_numeric,errors='raise')
    compares=[compare(d['mono4'],d['mono19'],'mono4','mono19'),
              compare(d['mono19'],d['poly19'],'mono19','poly19'),
              compare(d['poly19'],d['recovery'],'poly19','recovery')]
    ids=d['mono4'].index.intersection(d['mono19'].index)
    subset=d['mono4'].loc[ids]
    m=d['mono19'].loc[ids]
    chosen=ids[(subset.qval<.05)&(m.qval<.05)&(subset.b*m.b>0)].tolist()
    if not chosen:raise ValueError('No selected markers; check metadata')
    norm=d['norm'].apply(pd.to_numeric,errors='raise')
    filt=((norm>=10).sum(axis=1)>=6)
    expr=np.log2(1+norm.loc[filt]).astype('float64')
    if not set(chosen).issubset(set(expr.index)):raise ValueError('Panel genes not all represented in expression matrix')
    gene_table=[]
    for gene in chosen:
        row=dict(gene_id=gene,mono4_q=float(d['mono4'].loc[gene,'qval']),mono19_q=float(d['mono19'].loc[gene,'qval']),
                 poly19_q=float(d['poly19'].loc[gene,'qval']),mono4_b=float(d['mono4'].loc[gene,'b']),
                 mono19_b=float(d['mono19'].loc[gene,'b']),poly19_b=float(d['poly19'].loc[gene,'b']))
        for group,cols in GROUPS.items():row['mean_'+group]=float(expr.loc[gene,cols].mean())
        row['moves_toward_U_at_F90']=abs(row['mean_F90']-row['mean_U'])<abs(row['mean_T']-row['mean_U'])
        gene_table.append(row)
    ptab=pd.DataFrame(gene_table)
    group_scores={k:expr.loc[chosen,v].mean(axis=0).tolist() for k,v in GROUPS.items()}
    contrasts={k+'_minus_U':[float(a-b) for a,b in zip(v,group_scores['U'])]
               for k,v in group_scores.items() if k!='U'}
    cv={str(k):fold_select_counts(expr,k) for k in (10,25,50,100,250,500)}
    permutations=[]
    for mask in range(8):
        flips=tuple((mask>>i)&1 for i in range(3))
        margins=fold_select_counts(expr,50,flips=flips)
        permutations.append(dict(mask=mask,positive_folds=sum(x>0 for x in margins),margins=margins))
    p_exact=sum(x['positive_folds']>=permutations[0]['positive_folds'] for x in permutations)/8
    transfer=[]
    ext=d['mono19'].sort_values('pval')
    for k in (10,25,50,100,250,500,1000):
        sub=ext.index.intersection(expr.index)[:0]  # reserved, sorting completed below
        sub=ext.loc[ext.index.intersection(expr.index)].sort_values('pval').head(k)
        sg=np.sign(sub.b.to_numpy())
        mar=[float(np.mean(sg*(expr.loc[sub.index,GROUPS['T'][j]].to_numpy()-
                                 expr.loc[sub.index,GROUPS['U'][j]].to_numpy()))) for j in range(3)]
        transfer.append(dict(k=k,pair_margins=mar,positive_pairs=sum(x>0 for x in mar)))
    summary=dict(input_type='author-released processed normalized RNA-seq and Wald tables, not original FASTQ',
                 sample_unit='biological pool; not individual cell/donor',samples_per_group=3,
                 genes_normalized=int(norm.shape[0]),genes_expressed_filter=int(expr.shape[0]),
                 published_significant={k:int((d[k].qval<.05).sum()) for k in ('mono4','mono19','poly19','recovery')},
                 overlap_tests=compares,chosen_n=len(chosen),
                 chosen_genes=chosen,chosen_gene_same_direction_in_poly19=int((ptab.poly19_b>0).sum()),
                 gene_panel_group_mean={k:float(np.mean(v)) for k,v in group_scores.items()},
                 gene_panel_scores_per_sample=group_scores,panel_paired_contrasts=contrasts,
                 genes_toward_untrained_at_90_min=int(ptab.moves_toward_U_at_F90.sum()),
                 cv_within_poly= cv,cv_50_fold_label_permutation_exact_p=p_exact,
                 cv_50_permutations=permutations,mono19_to_poly_transfer=transfer,
                 proteomic_note='Independent proteomics sample-level matrix not contained in GitHub; no proteomics reprocessing represented as completed.')
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'real_biological_results.json').write_text(json.dumps(summary,indent=2))
    ptab.to_csv(outdir/'real_14gene_panel.csv',index=False)
    pd.DataFrame(transfer).to_csv(outdir/'cross_culture_transfer.csv',index=False)
    pd.DataFrame(permutations).to_csv(outdir/'within_cohort_exact_permutations.csv',index=False)
    print('REAL DATA SUMMARY:',json.dumps({k:summary[k] for k in ('genes_expressed_filter','published_significant','chosen_n','cv_50_fold_label_permutation_exact_p')},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input-dir',default=None);ap.add_argument('--out',default='stentor_reanalysis_outputs')
    args=ap.parse_args();run(args.input_dir,pathlib.Path(args.out))
