#!/usr/bin/env python3
"""Exploratory HLA-DRB1/3/4/5 allele-panel exact-kmer analysis of original public SRA reads.
Partial prefix of two FASTQ mates is NOT full allele-aware remapping / personalized genotype inference.
No human reads are retained in output; only aggregate diagnostics."""
import csv,gzip,hashlib,io,json,sys,time,urllib.request
from collections import Counter,defaultdict
from pathlib import Path
out=Path("hla_read_pilot");out.mkdir(exist_ok=True)
K=31
PREFIX=int(sys.argv[1]) if len(sys.argv)>1 else 268435456
wanted=['HN22','HBA22','HN35','HBA35']
manifest=list(csv.DictReader(open('sra_run_manifest.tsv'),delimiter='\t'))
byname={r['sample_title']:r for r in manifest}
assert all(x in byname for x in wanted)
db='https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_nuc.fasta'
request=urllib.request.Request(db,headers={'User-Agent':'Scientific-HLA-reanalysis/1.0'})
with urllib.request.urlopen(request,timeout=90) as f: fasta=f.read().decode()
fa_sha=hashlib.sha256(fasta.encode()).hexdigest()
seqs=[]
name=None;parts=[]
for line in fasta.splitlines():
 if line.startswith('>'):
  if name is not None:seqs.append((name,''.join(parts)))
  name=line[1:];parts=[]
 else:parts.append(line.strip().upper())
if name is not None:seqs.append((name,''.join(parts)))
gcounts=Counter(); bit={'DRB1':1,'DRB3':2,'DRB4':4,'DRB5':8}; index={}
for hdr,seq in seqs:
 gene=next((g for g in bit if g+'*' in hdr),None)
 if not gene:continue
 gcounts[gene]+=1
 for i in range(len(seq)-K+1):
  k=seq[i:i+K]
  if set(k).issubset(set('ACGT')):index[k]=index.get(k,0)|bit[gene]
assert all(gcounts[g]>0 for g in bit),gcounts
translate=str.maketrans('ACGTN','TGCAN')
def classify(seq):
 seq=seq.upper()
 seen=0
 for st in [0,20,40,60,80,100,120]:
  if st+K > len(seq):continue
  k=seq[st:st+K]
  seen |= index.get(k,0)
  seen |= index.get(k.translate(translate)[::-1],0)
 return seen
def fetchprefix(url,path):
 req=urllib.request.Request(url,headers={'Range':f'bytes=0-{PREFIX-1}','User-Agent':'Scientific-HLA-reanalysis/1.0'})
 status=None; n=0
 with urllib.request.urlopen(req,timeout=120) as resp,open(path,'wb') as out:
  status=resp.status
  while n<PREFIX:
   block=resp.read(min(1024*1024,PREFIX-n))
   if not block:break
   out.write(block);n+=len(block)
 return status,n
results=[]
for sam in wanted:
 r=byname[sam]; urls=r['fastq_ftp'].split(';')
 assert len(urls)==2
 counts=Counter(); nbytes=0; statuses=[]
 for mate,ftp in enumerate(urls,1):
  dest=out/f'temp_{sam}_R{mate}.fastq.gz.partial'
  status,n=fetchprefix('https://'+ftp,dest);statuses.append(status);nbytes+=n
  f=gzip.open(dest,'rt',errors='replace')
  try:
   while True:
    try:
     h=f.readline();q=f.readline();plus=f.readline();qual=f.readline()
    except (EOFError,OSError,UnicodeError):break
    if not h or not q or not plus or not qual:break
    if not h.startswith('@') or not plus.startswith('+'):continue
    counts['reads_evaluated']+=1
    mask=classify(q.strip())
    if mask:
     counts['DRB_family_kmer_supported_reads']+=1
     if mask in bit.values():
      gene=next(g for g,b in bit.items() if b==mask)
      counts['exclusive_'+gene]+=1
     else:counts['multi_locus_ambiguous']+=1
  finally:f.close()
  dest.unlink()
 summary={'sample':sam,'run':r['run_accession'],'mates':2,'http_statuses':statuses,
  'prefix_bytes_per_mate':PREFIX,'bytes_downloaded':nbytes,'read_counts':dict(counts)}
 results.append(summary)
 print(json.dumps(summary),flush=True)
report={'classification':'exploratory incomplete FASTQ-prefix allele-panel exact-31mer read-support audit, NOT donor HLA typing, complete read mapping, or validated HLA quantification',
'K':K,'positions_tested_per_read':[0,20,40,60,80,100,120],
'source':'GSE174409 SRP319708 ENA public FASTQ runs','allele_database_url':db,'allele_fasta_sha256':fa_sha,
'allele_count_by_locus':dict(gcounts),'distinct_index_kmers':len(index),
'limitations':['only first compressed prefix of each mate; not random read subsampling','exact matches miss sequencing errors and reads crossing noncanonical splice junctions','no donor allele calls; no read quality/model-based alignment or whole-transcriptome calibration','multi-locus reads cannot be uniquely assigned','cannot validate original HLA-DRB1/DRB5 per-donor read mapping or associated differential-expression FDR'], 'results':results}
(out/'HLA_FASTQ_PILOT.json').write_text(json.dumps(report,indent=2)+'\n')
print('REPORT '+str(out/'HLA_FASTQ_PILOT.json'))
