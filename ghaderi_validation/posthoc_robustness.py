import itertools, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

# Values are the audited outputs of the locked Ghaderi workflow (run 38192392565),
# itself extracted from the two hashed source files. This is a post-reveal robustness analysis.
rows=[
('A1','audio',61.794043,38.874605),('A1','delay',53.364070,7.785660),('A1','whisker',55.663611,12.758641),
('ALM','audio',55.710011,7.800500),('ALM','delay',58.713999,32.092592),('ALM','whisker',68.447492,52.341964),
('wM2','audio',57.666080,33.049962),('wM2','delay',61.741684,44.528599),('wM2','whisker',69.161736,48.247632),
('wS1','audio',54.306252,19.451337),('wS1','delay',51.942839,15.976635),('wS1','whisker',57.635138,37.640712),
('wS2','audio',54.921825,25.997587),('wS2','delay',51.662188,8.182996),('wS2','whisker',58.683634,35.487251)]
df=pd.DataFrame(rows,columns=['region','epoch','decoder_accuracy','causal_impairment_pp'])
REG=['A1','wS1','wS2','wM2','ALM']; EPS=['audio','delay','whisker']

def stat(tab,regs,eps):
    xs=[]; ys=[]
    for ep in eps:
        q=tab[tab.epoch==ep].set_index('region').loc[regs]
        rx=rankdata(q.decoder_accuracy.to_numpy(float)); ry=rankdata(q.causal_impairment_pp.to_numpy(float)); rx-=rx.mean(); ry-=ry.mean(); xs.append(rx); ys.append(ry)
    x=np.concatenate(xs); y=np.concatenate(ys); return float((x@y)/np.sqrt((x@x)*(y@y)))

def exact(tab,regs,eps):
    obs=stat(tab,regs,eps); vals=[]
    for perm in itertools.permutations(regs):
        mp=dict(zip(regs,perm)); ca=tab[['region','epoch','causal_impairment_pp']].copy(); ca['region']=ca.region.map(mp)
        t=tab.drop(columns=['causal_impairment_pp']).merge(ca,on=['region','epoch'])
        vals.append(stat(t,regs,eps))
    vals=np.asarray(vals); return obs,float(np.mean(vals>=obs-1e-12)),len(vals)

out={'analysis_status':'post-reveal influence/robustness analysis; not part of the locked Ghaderi primary test','leave_one_region_out':{},'leave_one_epoch_out':{}}
for r in REG:
    regs=[x for x in REG if x!=r]; q=df[df.region.isin(regs)]; s,p,n=exact(q,regs,EPS); out['leave_one_region_out'][r]={'statistic':s,'exact_one_sided_p':p,'n_permutations':n}
for ep in EPS:
    eps=[x for x in EPS if x!=ep]; q=df[df.epoch.isin(eps)]; s,p,n=exact(q,REG,eps); out['leave_one_epoch_out'][ep]={'statistic':s,'exact_one_sided_p':p,'n_permutations':n}
out['source_note']='Input values are rounded to six decimals from the audited locked-workflow output; results were checked against full-precision calculations before reporting.'
Path('analysis/ghaderi_matched_validation').mkdir(parents=True,exist_ok=True)
Path('analysis/ghaderi_matched_validation/posthoc_robustness.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
