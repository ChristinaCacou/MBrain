#!/usr/bin/env python3
"""Audit EGR1 rows in Wei 2023 independently recruited 95-donor ventral midbrain public supplement.
Published gene-level statistics are not donor-level raw counts and do not constitute our rerun."""
from pathlib import Path
import json,re,hashlib,openpyxl
out=Path("wei_2023_egr1_audit");out.mkdir(exist_ok=True)
matches=[]
summary=[]
for f in sorted(Path(".").glob("Wei2023_*.xlsx")):
  if not f.exists():continue
  wb=openpyxl.load_workbook(f,read_only=True,data_only=True)
  sha=hashlib.sha256(f.read_bytes()).hexdigest()
  meta={"file":f.name,"sha256":sha,"bytes":f.stat().st_size,"sheets":[]}
  for sh in wb.worksheets:
    count=0; hit=0; header=None
    for i,row in enumerate(sh.iter_rows(values_only=True),1):
      r=list(row)
      if i==1:header=[str(x)[:130] if x is not None else "" for x in r[:35]]
      count+=1
      cols=[j for j,v in enumerate(r) if isinstance(v,str) and (v.strip().upper()=="EGR1" or re.fullmatch(r"ENSG00000120738(?:\\.[0-9]+)?",v.strip()))]
      if cols:
        hit+=1
        matches.append({"file":f.name,"sheet":sh.title,"row":i,"matched_columns":cols,
        "header":header,"cells":[str(x)[:220] if x is not None else "" for x in r[:45]]})
    meta["sheets"].append({"name":sh.title,"rows_scanned":count,"EGR1_rows":hit,"header":header})
  summary.append(meta)
result={"study":"Wei et al Nat Commun 2023 n95 human VM, GSE240457","scope":"direct extraction of published supplementary result tables, NOT a novel pseudobulk refit","source_files":summary,"EGR1_rows":matches}
(out/"Published_Wei2023_EGR1_Extraction.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"files":len(summary),"EGR1_matches":len(matches),"rows":matches[:30],"file_metadata":summary},indent=2)[:28000])
