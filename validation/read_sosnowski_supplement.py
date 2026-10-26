#!/usr/bin/env python3
from pathlib import Path
import csv,json,re,sys,hashlib
p=Path(sys.argv[1]);b=p.read_bytes()
if b'<html' in b[:500].lower() or b'captcha' in b[:2000].lower():raise ValueError('Not a real data table; HTML challenge')
lines=p.read_text(encoding='utf-8',errors='replace').splitlines()
if len(lines)<100:raise ValueError('Too few lines for published full differential table')
hits=[]
for ix,line in enumerate(lines):
  cells=next(csv.reader([line],delimiter='\t')) if '\t' in line else next(csv.reader([line]))
  if any(re.fullmatch('EGR1',c.strip().strip('"'),flags=re.I) for c in cells):
    hits.append({'line_number':ix+1,'cells':cells[:25]})
header=lines[0][:1000]
summary={'source':'Sosnowski et al 2022 DLPFC opioid overdose published Supplementary Table S1 (mmc5.txt)','url':sys.argv[2],'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b),'lines':len(lines),'header':header,'EGR1_matching_rows':hits,'limitation':'Extracted published summary results, not a new raw-read or adjusted-count reanalysis; already included in Carter et al 2026 pooled meta-analysis.'}
Path('sosnowski_output').mkdir(exist_ok=True)
Path('sosnowski_output/EGR1_from_published_Sosnowski_TableS1.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2)[:7000])
if not hits:print('NO EGR1 row matched; report missing, do not invent')
