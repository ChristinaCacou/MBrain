import ast, json
from itertools import permutations
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path('analysis/pui_blind_fixedC')
OUT=ROOT/'diagnostics'; OUT.mkdir(parents=True,exist_ok=True)
COMMIT='a17677ef2ab38c198761ad66f8a8548e4b857901'
URL=f'https://raw.githubusercontent.com/elliottower/neural-geometry-reliability/{COMMIT}/experiments/shared_bundle.py'


def effects():
    txt=urlopen(URL).read().decode(); tree=ast.parse(txt)
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SILENCING_EFFECTS' for t in node.targets):
            return {str(k):float(v) for k,v in ast.literal_eval(node.value).items()}
    raise RuntimeError('SILENCING_EFFECTS not found')


def rho(x,y):
    x=rankdata(np.asarray(x,float)); y=rankdata(np.asarray(y,float)); x-=x.mean(); y-=y.mean()
    d=np.sqrt((x@x)*(y@y)); return float(x@y/d) if d else np.nan


def exact_p_positive(x,y):
    obs=rho(x,y); yy=tuple(np.asarray(y,float)); n=0; ge=0
    for p in permutations(yy):
        n+=1; ge += rho(x,p) >= obs-1e-15
    return obs, ge/n


def main():
    f=pd.read_csv(ROOT/'blind_features_region.csv')
    o=effects(); d=pd.DataFrame({'region':list(o),'silencing_effect':list(o.values())}).merge(f,on='region')
    rows=[]
    for threshold in [1,3,5,7]:
        s=d[d.n_sessions>=threshold].copy().sort_values('region')
        if len(s)<=8:
            r,p=exact_p_positive(s.pui_mean,s.silencing_effect)
            method='exact all label permutations'
        else:
            r=rho(s.pui_mean,s.silencing_effect); p=np.nan; method='descriptive only (>8 regions)'
        rows.append({'min_sessions':threshold,'n_regions':len(s),'regions':','.join(s.region),'spearman_rho':r,'one_sided_permutation_p':p,'method':method})
    pd.DataFrame(rows).to_csv(OUT/'coverage_sensitivity.csv',index=False)
    (OUT/'coverage_sensitivity.json').write_text(json.dumps({'analysis_type':'post-reveal diagnostic; primary analysis unchanged','rows':rows},indent=2))
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__': main()
