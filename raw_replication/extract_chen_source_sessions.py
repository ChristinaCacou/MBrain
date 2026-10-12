from pathlib import Path
import re, zipfile, urllib.request, json
url='https://zenodo.org/records/6713616/files/analyses_scripts.zip?download=1'
z=Path('analysis_scripts.zip')
urllib.request.urlretrieve(url,z)
out=Path('analysis/chen_source_session_inspection');out.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(z) as zz: zz.extractall('analysis_scripts')
pat=re.compile(r'(?:data_structure_)?(BAYLORGC\d+)[_/\\-]?(\d{4})[_-]?(\d{2})[_-]?(\d{2})',re.I)
records=[]; contexts=[]
for p in Path('analysis_scripts').rglob('*.m'):
    txt=p.read_text(errors='ignore'); lines=txt.splitlines()
    for i,line in enumerate(lines):
        if 'BAYLORGC' in line or 'standard' in line.lower() or 'session' in line.lower() and ('task' in line.lower() or 'mouse' in line.lower()):
            if 'BAYLORGC' in line:
                contexts.append({'file':str(p),'line':i+1,'text':line[:1000]})
            for m in pat.finditer(line):
                records.append({'animal':m.group(1).upper(),'date':f'{m.group(2)}_{m.group(3)}_{m.group(4)}','file':str(p),'line':i+1})
uniq=[];seen=set()
for r in records:
    k=(r['animal'],r['date'])
    if k not in seen:seen.add(k);uniq.append(r)
(out/'session_ids.json').write_text(json.dumps(uniq,indent=2))
(out/'baylor_context.json').write_text(json.dumps(contexts,indent=2))
# Also save targeted relevant-script excerpts around session/task selection tokens.
with open(out/'relevant_excerpts.txt','w') as f:
    for p in Path('analysis_scripts').rglob('*.m'):
        if p.name in ['analysis_compile_all_data_2ALM_CDprojection.m','analysis_compile_all_data_2ALM_selectivity.m']:
            lines=p.read_text(errors='ignore').splitlines();f.write(f'### {p}\n')
            for i,line in enumerate(lines):
                low=line.lower()
                if any(k in low for k in ['baylor','standard','session','task','mouse_name','animal_name','data_structure']):
                    a=max(0,i-2);b=min(len(lines),i+3)
                    f.write(f'-- around line {i+1} --\n')
                    for j in range(a,b):f.write(f'{j+1:05d}: {lines[j]}\n')
print(json.dumps({'n_unique_ids':len(uniq),'n_context_lines':len(contexts)},indent=2))
