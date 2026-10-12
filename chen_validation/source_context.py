from pathlib import Path
import zipfile, urllib.request
URL='https://zenodo.org/records/6713616/files/analyses_scripts.zip?download=1'
z=Path('chen_validation/analyses_scripts.zip')
if not z.exists(): urllib.request.urlretrieve(URL,z)
r=Path('chen_validation/analyses_scripts'); r.mkdir(parents=True,exist_ok=True)
if not any(r.iterdir()):
    with zipfile.ZipFile(z) as zz: zz.extractall(r)
for fn, ranges in {
 'analysis_compile_all_data_2ALM_CDprojection.m':[(610,720),(920,990)],
 'analysis_compile_all_data_2ALM_selectivity.m':[(670,780),(1100,1210)]
}.items():
    p=next(r.rglob(fn)); lines=p.read_text(errors='ignore').splitlines()
    print('\n###',fn)
    for a,b in ranges:
        print(f'\n--- lines {a}-{b} ---')
        for i in range(a-1,min(b,len(lines))): print(f'{i+1:04d}: {lines[i]}')
