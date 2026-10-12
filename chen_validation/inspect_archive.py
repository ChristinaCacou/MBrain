from pathlib import Path
import json, re, zipfile, urllib.request

URL='https://zenodo.org/records/6713616/files/analyses_scripts.zip?download=1'
out=Path('analysis/chen_archive_inspection'); out.mkdir(parents=True,exist_ok=True)
zip_path=Path('chen_validation/analyses_scripts.zip')
if not zip_path.exists():
    urllib.request.urlretrieve(URL, zip_path)
root=Path('chen_validation/analyses_scripts')
root.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(zip_path) as z:
    z.extractall(root)
files=[p for p in root.rglob('*') if p.is_file()]
(out/'file_list.txt').write_text('\n'.join(str(p.relative_to(root)) for p in files))
patterns=[r'recovery',r'interhemispheric',r'Figure\s*5',r'fig\s*5',r'deltaR',r'dR2',r'behavior',r'modular']
hits=[]
for p in files:
    if p.suffix.lower() not in {'.m','.txt','.md','.py','.csv','.json'}: continue
    try: txt=p.read_text(errors='ignore')
    except Exception: continue
    lines=txt.splitlines()
    for i,line in enumerate(lines):
        if any(re.search(pat,line,re.I) for pat in patterns):
            lo=max(0,i-2); hi=min(len(lines),i+3)
            hits.append({'file':str(p.relative_to(root)),'line':i+1,'context':'\n'.join(lines[lo:hi])})
(out/'keyword_hits.json').write_text(json.dumps(hits[:1000],indent=2))
# summarize binary/Mat files that may be processed intermediates
matlike=[str(p.relative_to(root)) for p in files if p.suffix.lower() in {'.mat','.h5','.hdf5','.csv'}]
(out/'processed_candidates.txt').write_text('\n'.join(matlike))
print('files',len(files),'hits',len(hits),'processed_candidates',len(matlike))
for h in hits[:80]:
    print('\n---',h['file'],h['line'],'---\n',h['context'])
