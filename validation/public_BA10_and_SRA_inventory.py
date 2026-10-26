#!/usr/bin/env python3
"""Public-only BA10 collection and independently public GSE174409 SRA run manifest.
Does not attempt to fetch restricted human data."""
import json,urllib.parse,urllib.request,csv,io,hashlib
from pathlib import Path
out=Path("public_data_access_audit");out.mkdir(exist_ok=True)
def get(url,accept="application/json"):
 req=urllib.request.Request(url,headers={"User-Agent":"OUD-Research-Public-Data-Audit/1.0","Accept":accept})
 try:
  with urllib.request.urlopen(req,timeout=45) as r:
   body=r.read(3_000_000)
   return {"status":r.status,"final_url":r.url,"length":len(body),"sha256":hashlib.sha256(body).hexdigest(),"text":body.decode("utf-8",errors="replace")}
 except Exception as e:return {"error":repr(e),"url":url}
urls={
"nemo_col_plain":"https://assets.nemoarchive.org/api/collection/col-c9h2ph2",
"nemo_col_prefixed":"https://assets.nemoarchive.org/api/collection/nemo:col-c9h2ph2",
"nemo_col_website":"https://assets.nemoarchive.org/col-c9h2ph2",
"scorch_pub_portal":"https://scorch.igs.umaryland.edu",
"scorch_brainome":"https://brainome.ucsd.edu/SCORCH",
"ena_sra_manifest":"https://www.ebi.ac.uk/ena/portal/api/filereport?"+urllib.parse.urlencode({"accession":"SRP319708","result":"read_run","fields":"run_accession,experiment_accession,sample_accession,sample_title,fastq_ftp,fastq_bytes,library_layout","format":"tsv"}),
}
results={}
for key,url in urls.items():
 x=get(url,"text/plain" if key=="ena_sra_manifest" else "application/json")
 (out/(key+".txt")).write_text(x.get("text",x.get("error",""))[:3000000])
 results[key]={k:v for k,v in x.items() if k!="text"}
 if "text" in x:
  s=x["text"]
  if key=="ena_sra_manifest":
   rows=list(csv.DictReader(io.StringIO(s),delimiter="\t"))
   results[key]["rows"]=len(rows)
   results[key]["sample_titles_first_10"]=[{k:r.get(k) for k in ("run_accession","experiment_accession","sample_title","fastq_bytes")} for r in rows[:10]]
   (out/"sra_run_manifest.tsv").write_text(s)
  if key.startswith("nemo_col") and x.get("status")==200:
   try:
    o=json.loads(s)
    results[key]["json_keys"]=list(o) if isinstance(o,dict) else []
    results[key]["public_download_candidate_strings"]=list({v for v in __import__('re').findall(r'https?://[^"\\\s]{10,200}',s) if any(ext in v.lower() for ext in ('.h5ad','.csv','.tsv','.rds','.h5','.gz','.zip'))})[:30]
   except Exception:pass
print(json.dumps(results,indent=2))
(out/"PUBLIC_ACCESS_AUDIT.json").write_text(json.dumps(results,indent=2)+"\n")
