from pathlib import Path
import urllib.request, zipfile, re, json, os
ROOT=Path('chen_validation/raw_source'); ROOT.mkdir(parents=True,exist_ok=True)
zip_path=ROOT/'analyses_scripts.zip'
if not zip_path.exists():
    urllib.request.urlretrieve('https://zenodo.org/records/6713616/files/analyses_scripts.zip?download=1',zip_path)
ex=ROOT/'scripts'; ex.mkdir(exist_ok=True)
with zipfile.ZipFile(zip_path) as z:
    z.extractall(ex)
# also fetch compact metadata helper directly
try:
    urllib.request.urlretrieve('https://zenodo.org/records/6713616/files/meta_animal_info.m?download=1',ROOT/'meta_animal_info.m')
except Exception as e:
    print('meta download failed',e)
# collect text sources and lines that mention standard task/session lists or BAYLOR identifiers
hits=[]
for p in list(ex.rglob('*.m'))+[ROOT/'meta_animal_info.m']:
    if not p.exists(): continue
    try: txt=p.read_text(errors='ignore')
    except: continue
    for i,line in enumerate(txt.splitlines(),1):
        lo=line.lower()
        if ('standard' in lo or 'baylorgc' in lo or 'filelist' in lo or 'exp_list' in lo or 'session' in lo and ('animal' in lo or 'task' in lo)):
            hits.append({'file':str(p.relative_to(ROOT)),'line':i,'text':line[:500]})
# list source files likely responsible for Fig2/3/4/5/6 and metadata
files=[str(p.relative_to(ROOT)) for p in ex.rglob('*.m') if any(k in p.name.lower() for k in ['figure2','figure3','figure4','figure5','figure6','meta','load_data','modular','rank'])]
(ROOT/'source_hits.json').write_text(json.dumps(hits,indent=2))
(ROOT/'candidate_files.json').write_text(json.dumps(files,indent=2))
print('hits',len(hits),'files',len(files))
for h in hits[:500]: print(f"{h['file']}:{h['line']}: {h['text']}")
