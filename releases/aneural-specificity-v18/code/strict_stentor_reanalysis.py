#!/usr/bin/env python3
"""Strict, provenance-preserving secondary reanalysis of two Stentor RNA-seq cohorts.

Never estimates missing biological values. RNA abundance quantification must already
have been executed on *all* relevant libraries using one validated reference.
"""
import argparse,csv,hashlib,json,math,re,sys,zipfile,subprocess,shutil
from pathlib import Path
import numpy as np
import pandas as pd

R_IDS=[f'DR4876_{k:02d}' for k in range(1,22)]
M_IDS=[f'SRR{k}' for k in range(39453085,39453106)]
CONTROL_IDS=['DR4876_11','DR4876_13','DR4876_15','DR4876_17','DR4876_19','DR4876_21']
G={
'T':['DR4876_01','DR4876_04','DR4876_07'],
'U':['DR4876_03','DR4876_06','DR4876_09'],
'F10':['DR4876_10','DR4876_12','DR4876_14'],
'U10':['DR4876_11','DR4876_13','DR4876_15'],
'F90':['DR4876_16','DR4876_18','DR4876_20'],
'U90':['DR4876_17','DR4876_19','DR4876_21']}


def write_json(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')

def audit(args):
    obs=pd.read_csv(args.matrix,nrows=1)
    names=set(obs.columns)
    found=sorted(set(R_IDS)&names)
    missing=sorted(set(R_IDS)-names)
    d={'dataset':'Rajan et al. 2026 polyculture RNA-seq','file':str(args.matrix),
       'observed_columns':found,'missing_columns':missing,
       'archived_control_sample_ids':CONTROL_IDS,
       'has_required_six_controls':all(x in names for x in CONTROL_IDS),
       'can_estimate_did':all(x in names for x in R_IDS) and len(pd.read_csv(args.matrix,nrows=1).columns)>20,
       'status':'READY' if set(R_IDS)<=names else 'STOP: missing measurements; no estimates or imputation performed'}
    write_json(args.out,d);print(json.dumps(d,indent=2))
    return 0 if not missing else 2

def zip_audit(args):
    with zipfile.ZipFile(args.zip) as z:
        members=z.infolist(); matches={}
        for sample in R_IDS:
            v=[]
            for q in members:
                nm=q.filename.replace('\\','/'); base=nm.split('/')[-1]
                if re.search(r'(?<!\d)'+re.escape(sample)+r'(?!\d)',nm,re.I) and re.search(r'R[12]',base,re.I) and '.fastq' in base.lower():
                    v.append({'member':q.filename,'compressed_bytes':q.compress_size,'bytes':q.file_size})
            matches[sample]=v
    audit={'archive':str(args.zip),'size_bytes':Path(args.zip).stat().st_size,
           'members':len(members),'fastq_members':matches,
           'missing_r1_r2':[s for s,v in matches.items() if not any(re.search('R1',x['member']) for x in v) or not any(re.search('R2',x['member']) for x in v)]}
    write_json(args.out,audit);print(json.dumps({'members':len(members),'samples_detected':sum(bool(x) for x in matches.values()),'missing':audit['missing_r1_r2']},indent=2))
    return 0 if not audit['missing_r1_r2'] else 2

def hash_file(p):
    s=hashlib.sha256()
    with open(p,'rb') as f:
        for q in iter(lambda:f.read(4*1024*1024),b''):s.update(q)
    return s.hexdigest()

def quantify(args):
    """Requires already acquired 209.8 GB archive and verified original Kallisto index.
    Does not run unless explicitly given --execute."""
    if not args.execute:
        print('DRY-RUN ONLY. --execute required; no download or biological analysis performed.');return 2
    if not Path(args.index).is_file():raise SystemExit('STOP: exact/reference-verified Kallisto index missing')
    if shutil.which('kallisto') is None:raise SystemExit('STOP: kallisto not found on PATH')
    dest=Path(args.outdir);dest.mkdir(parents=True,exist_ok=True)
    selected=CONTROL_IDS if args.controls_only else R_IDS
    with zipfile.ZipFile(args.zip) as z:
        inf=z.infolist()
        for sid in selected:
            candidate={}
            for mate in (1,2):
                v=[x for x in inf if sid.lower() in x.filename.lower() and re.search(r'R'+str(mate)+r'(?:_001)?\.(?:fastq|fq)\.gz$',x.filename,re.I)]
                if len(v)!=1:raise SystemExit(f'STOP: {sid} read {mate} has {len(v)} matching members; inspect archive names manually')
                candidate[mate]=v[0]
            out=dest/sid
            if (out/'abundance.tsv').exists():
                print(sid,'already quantified, skipping');continue
            free=shutil.disk_usage(dest).free
            need=sum(x.file_size for x in candidate.values())+2*1024**3
            if free<need:raise SystemExit(f'STOP: {sid} needs {need/1024**3:.1f} GB temporary free space; have {free/1024**3:.1f}')
            tmp=dest/'_temporary_fastq';tmp.mkdir(exist_ok=True)
            reads=[]
            for mate in (1,2):
                fn=tmp/f'{sid}_R{mate}.fastq.gz';reads.append(fn)
                with z.open(candidate[mate]) as src,open(fn,'wb') as target:
                    shutil.copyfileobj(src,target,length=4*1024*1024)
            cmd=['kallisto','quant','-i',str(Path(args.index).resolve()),'-o',str(out),'-b','30','-t',str(args.threads)]+[str(x) for x in reads]
            rc=subprocess.run(cmd,check=False)
            for f in reads:f.unlink(missing_ok=True)
            if rc.returncode or not (out/'abundance.tsv').exists():raise SystemExit(f'STOP: kallisto failed on {sid}')
            write_json(out/'secondary_quantification_provenance.json',{
              'sample_id':sid,'archive_sha256':hash_file(args.zip) if args.hash_entire_archive else 'not computed',
              'index_sha256':hash_file(args.index), 'original_archive_members':[candidate[1].filename,candidate[2].filename],
              'method':'kallisto quant bootstrap=30','study':'Rajan 2026','note':'Secondary processing; reference must be validated against original study.'})
    print('Completed requested quantifications. Do not infer a recovery effect until all 21 sample estimates are jointly processed.')
    return 0

def matrix(args):
    """Combine all sample-level kallisto estimated counts; log-normalize jointly.
    Requires all 21. Inference is descriptive due n=3 pools/condition.
    """
    ix=[]; arr=[]; xs={}
    for sample in R_IDS:
        p=Path(args.quantdir)/sample/'abundance.tsv'
        if not p.exists():raise SystemExit(f'STOP: missing measured abundance.tsv: {sample}. Do not impute.')
        x=pd.read_csv(p,sep='\t')
        if not {'target_id','est_counts'}<=set(x.columns):raise SystemExit(f'STOP: {p} lacks target_id/est_counts')
        if x.target_id.duplicated().any():raise SystemExit(f'STOP: duplicate transcript target IDs in {sample}')
        xs[sample]=x.set_index('target_id')['est_counts']
    raw=pd.concat(xs,axis=1)
    if raw.isna().any().any() or not np.isfinite(raw.to_numpy()).all() or (raw.to_numpy()<0).any():raise SystemExit('STOP: incompatible references/invalid counts')
    if raw.shape[0]<1000:raise SystemExit('STOP: full transcriptome needed, not a 14-gene extract')
    # Median-ratio library size factors based on positive counts present in all samples.
    vals=raw.to_numpy(dtype=float); complete=(vals>0).all(axis=1)
    if complete.sum()<200:raise SystemExit('STOP: insufficient universally positive genes for stable median-ratio size factors')
    logs=np.log(vals[complete]);gm=np.exp(logs.mean(axis=1));sf=np.median(vals[complete]/gm[:,None],axis=0)
    if not np.isfinite(sf).all() or (sf<=0).any():raise SystemExit('STOP: invalid size factors')
    normalized=raw/sf
    # The 14-gene signature is from separate previously analyzed discovery data.
    panel=pd.read_csv(args.panel)
    if 'transcript_id' not in panel.columns:raise SystemExit('STOP: panel CSV must contain transcript_id')
    genes=list(panel.transcript_id.astype(str))
    missing=sorted(set(genes)-set(raw.index))
    if missing:raise SystemExit('STOP: gene IDs absent / changed reference: '+','.join(missing))
    log=np.log2(normalized+1)
    base=log.loc[genes].mean(axis=0)
    centered=base-log.mean(axis=0)
    lines=[]
    for metric,score in [('raw_log2',base),('centered_log2',centered)]:
        means={k:float(score[v].mean()) for k,v in G.items()}
        for minute in (10,90):
            before=means[f'F{minute}']-means['T']
            drift=means[f'U{minute}']-means['U']
            lines.append({'score':metric,'minutes':minute,'trained_recovery_delta':before,
                'untreated_time_drift':drift,'difference_in_differences':before-drift,
                'n_pools_per_condition':3,'inference':'DESCRIPTIVE; no pairing or causal assumptions verified'})
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    raw.to_csv(out/'all21_estimated_counts.csv',index_label='transcript_id')
    normalized.to_csv(out/'all21_jointly_normalized_counts.csv',index_label='transcript_id')
    pd.DataFrame({'sample':R_IDS,'size_factor':sf}).to_csv(out/'size_factors.csv',index=False)
    pd.DataFrame({'sample':R_IDS,'raw_panel':base.to_numpy(),'centered_panel':centered.to_numpy()}).to_csv(out/'panel_sample_scores.csv',index=False)
    pd.DataFrame(lines).to_csv(out/'recovery_contrasts_descriptive.csv',index=False)
    write_json(out/'analysis_status.json',{'status':'DESCRIPTIVE_REANALYSIS_COMPLETE','n_samples':21,'n_genes':raw.shape[0],
      'controls_observed':CONTROL_IDS,'sensitivity':'uncertainty across only three pooled replicates not reliably estimated',
      'limitations':['not the original sleuth normalization','source index identity and all biological metadata require audit',
      'genome reference can differ from independent Miao cohort','no claim about consciousness or individual-cell memory']})
    print(pd.DataFrame(lines).to_string(index=False));return 0

def independent(args):
    """Exploratory orthology-checked transfer of a fixed Rajan signature to Miao time series.

    Requires *measured* all-21 sample count matrix, original author sample metadata,
    and high-confidence one-to-one sequence-based orthology mapping. No assumptions
    are made about gene-name equivalence across independently assembled genomes.
    """
    if not args.miao_counts or not args.crosswalk or not args.metadata or not args.panel:
        raise SystemExit('STOP: measured Miao matrix, validated orthology crosswalk, seven-stage sample metadata and frozen Rajan panel are all required.')
    m=pd.read_csv(args.miao_counts);x=pd.read_csv(args.crosswalk);meta=pd.read_csv(args.metadata);panel=pd.read_csv(args.panel)
    if not {'Miao_gene_id','Rajan_transcript_id','mapping_confidence'}<=set(x.columns):raise SystemExit('STOP: crosswalk schema incomplete')
    if not {'run_accession','timepoint_order'}<=set(meta.columns):raise SystemExit('STOP: metadata must include run_accession and timepoint_order')
    if 'transcript_id' not in panel.columns:raise SystemExit('STOP: Rajan frozen panel must include transcript_id')
    if x.Miao_gene_id.duplicated().any() or x.Rajan_transcript_id.duplicated().any():raise SystemExit('STOP: many-to-many orthology mapping requires manual resolution')
    if not all(x.mapping_confidence.astype(str).str.lower()=='high'):raise SystemExit('STOP: only manually reviewed high-confidence one-to-one mapping accepted')
    if set(meta.run_accession)!=set(M_IDS) or meta.run_accession.duplicated().any():raise SystemExit('STOP: incomplete or duplicated accession metadata')
    if set(m.columns)&set(M_IDS)!=set(M_IDS):raise SystemExit(f'STOP: expected all 21 measured Miao libraries, found {len(set(m.columns)&set(M_IDS))}')
    if 'gene_id' not in m.columns:raise SystemExit('STOP: gene_id column absent')
    if len(m)<1000:raise SystemExit('STOP: complete transcriptome is needed for stable size factors')
    group_sizes=meta.groupby('timepoint_order').run_accession.nunique()
    if len(group_sizes)!=7 or not all(group_sizes==3):raise SystemExit('STOP: 7 stages x 3 samples not confirmed')
    frozen=set(panel.transcript_id.astype(str));mappings=x[x.Rajan_transcript_id.astype(str).isin(frozen)]
    if set(mappings.Rajan_transcript_id.astype(str))!=frozen:raise SystemExit('STOP: frozen panel lacks complete one-to-one mapping across reference genomes')
    if m.gene_id.duplicated().any():raise SystemExit('STOP: duplicate gene IDs')
    cnt=m.set_index('gene_id')[M_IDS]
    if not np.isfinite(cnt.to_numpy(dtype=float)).all() or (cnt.to_numpy(dtype=float)<0).any():raise SystemExit('STOP: non-finite or negative measured counts')
    vals=cnt.to_numpy(dtype=float);mask=(vals>0).all(axis=1)
    if mask.sum()<200:raise SystemExit('STOP: insufficient positive genes for size factors')
    gm=np.exp(np.log(vals[mask]).mean(axis=1));sf=np.median(vals[mask]/gm[:,None],axis=0)
    if not np.isfinite(sf).all() or (sf<=0).any():raise SystemExit('STOP: invalid size factors')
    norm=np.log2(cnt/sf+1)
    gene_ids=mappings.Miao_gene_id.astype(str).tolist()
    missing=set(gene_ids)-set(norm.index)
    if missing:raise SystemExit('STOP: mapped Miao genes absent from full measured matrix: '+','.join(sorted(missing)))
    score=norm.loc[gene_ids].mean(axis=0)
    dat=meta.copy().set_index('run_accession').loc[M_IDS]
    dat['mapped_panel_mean_log2']=score.values
    stages=dat.groupby('timepoint_order').mapped_panel_mean_log2.mean().sort_index()
    from scipy.stats import spearmanr
    obs=float(spearmanr(list(range(7)),stages.values).statistic)
    import itertools
    vals=stages.values
    ranks=np.argsort(np.argsort(vals)).astype(float)
    target=np.arange(7).astype(float)
    # Under sharp exchangeability of the seven *stage labels* only; exploratory.
    den=np.sqrt(((target-target.mean())**2).sum()*((ranks-ranks.mean())**2).sum())
    p=0;n=0
    for perm in itertools.permutations(ranks):
        pr=np.asarray(perm,dtype=float)
        rho=float(((target-target.mean())*(pr-pr.mean())).sum()/den)
        p+=(abs(rho)>=abs(obs)-1e-12);n+=1
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    dat.reset_index().to_csv(out/'measured_miao_panel_scores.csv',index=False)
    write_json(out/'independent_miao_exploratory_result.json',{
      'n_measured_samples':21,'n_mapped_genes':len(gene_ids),'n_timepoints':7,
      'stage_means':{str(i):float(v) for i,v in stages.items()},
      'spearman_stage_mean_trend':obs,'exploratory_stage_label_permutation_p':p/n,
      'limitations':['Cross-study orthology only after independently verified sequence mapping',
      'Miao time course is not the Rajan untreated-recovery experiment',
      'Temporal correlation is not a molecular memory mechanism, causal validation, or proof of consciousness',
      'Stage-level exchangeability not guaranteed; permutation result sensitivity only']})
    print(f'EXPLORATORY measured independent analysis complete: 21 samples, {len(gene_ids)} mapped genes, rho={obs:.3f}, stage-label p={p/n:.4f}')
    return 0

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='cmd',required=True)
    a=sub.add_parser('audit-processed');a.add_argument('--matrix',required=True);a.add_argument('--out',required=True);a.set_defaults(fn=audit)
    a=sub.add_parser('audit-archive');a.add_argument('--zip',required=True);a.add_argument('--out',required=True);a.set_defaults(fn=zip_audit)
    a=sub.add_parser('quantify-rajan');a.add_argument('--zip',required=True);a.add_argument('--index',required=True);a.add_argument('--outdir',required=True);a.add_argument('--threads',type=int,default=8);a.add_argument('--controls-only',action='store_true');a.add_argument('--execute',action='store_true');a.add_argument('--hash-entire-archive',action='store_true');a.set_defaults(fn=quantify)
    a=sub.add_parser('assemble-rajan');a.add_argument('--quantdir',required=True);a.add_argument('--panel',required=True);a.add_argument('--outdir',required=True);a.set_defaults(fn=matrix)
    a=sub.add_parser('independent-miao');a.add_argument('--miao-counts');a.add_argument('--crosswalk');a.add_argument('--metadata');a.add_argument('--panel');a.add_argument('--outdir',default='miao_exploratory_outputs');a.set_defaults(fn=independent)
    args=p.parse_args();sys.exit(args.fn(args))
if __name__=='__main__':main()
