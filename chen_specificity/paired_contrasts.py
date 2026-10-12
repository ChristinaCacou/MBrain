import itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
IN=Path('analysis_inputs/chen_mouse_rows.csv'); OUT=Path('analysis/chen_specificity'); OUT.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(IN); x=df.coupling.to_numpy(float); rng=np.random.default_rng(20261012)
def rho(y): return float(spearmanr(x,y).statistic)
def exact_swap(a,b):
 a=np.asarray(a,float); b=np.asarray(b,float); obs=rho(a)-rho(b); n=len(a); ge=tw=0; vals=[]
 # exact 2^n within-mouse outcome-label swaps
 for mask in range(1<<n):
  bits=((mask>>np.arange(n))&1).astype(bool); aa=np.where(bits,b,a); bb=np.where(bits,a,b); d=rho(aa)-rho(bb); vals.append(d); ge += d>=obs-1e-15; tw += abs(d)>=abs(obs)-1e-15
 vals=np.asarray(vals)
 # paired mouse bootstrap
 boots=[]
 for _ in range(50000):
  idx=rng.integers(0,n,n); xx=x[idx]; aa=a[idx]; bb=b[idx]
  boots.append(float(spearmanr(xx,aa).statistic-spearmanr(xx,bb).statistic))
 return {'delta_rho':obs,'exact_one_sided_p':ge/(1<<n),'exact_two_sided_p':tw/(1<<n),'bootstrap_95_ci':[float(np.nanpercentile(boots,2.5)),float(np.nanpercentile(boots,97.5))],'n_mice':n}
res={
 'right_behavior_vs_left_behavior':exact_swap(df.right_behavior_recovery,df.left_behavior_recovery),
 'right_behavior_vs_right_neuronal':exact_swap(df.right_behavior_recovery,df.right_neuronal_recovery),
 'right_neuronal_vs_left_neuronal':exact_swap(df.right_neuronal_recovery,df.left_neuronal_recovery),
 'interpretation':'Post hoc paired specificity contrasts. Significant association in one endpoint and nonsignificant association in another does not itself establish a significant difference; these exact within-mouse swap tests directly test the correlation differences.'
}
(OUT/'paired_contrasts.json').write_text(json.dumps(res,indent=2)); print(json.dumps(res,indent=2))
