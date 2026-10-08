import io, json, math, zipfile
from pathlib import Path
from collections import defaultdict
import requests, numpy as np, pandas as pd, networkx as nx
from scipy.io import loadmat
from scipy.stats import wilcoxon, spearmanr
from bs4 import BeautifulSoup
SEED=20261008; GAIN=.85; M0=5.; TAU=.05; DT=.01; FREQS=np.arange(0,100.0001,.5)
out=Path('v9_results'); out.mkdir(exist_ok=True); (out/'graphs').mkdir(exist_ok=True); (out/'ccep_subject_matrices').mkdir(exist_ok=True)
def dl(u,p):
 r=requests.get(u,timeout=180,allow_redirects=True,headers={'User-Agent':'Mozilla/5.0'}); r.raise_for_status(); Path(p).write_bytes(r.content)
def rp(G): return sum(1 for u,v in G.edges if u<v and G.has_edge(v,u))
def rf(G): return 2*rp(G)/G.number_of_edges() if G.number_of_edges() else np.nan
def opt(G,maxi=True,nprop=6000,seed=0,strong=False):
 rg=np.random.default_rng(seed); H=G.copy(); cur=rp(H); ed=list(H.edges())
 for _ in range(nprop):
  if len(ed)<2: break
  i,j=rg.choice(len(ed),2,replace=False); a,b=ed[i]; c,d=ed[j]
  if len({a,b,c,d})<4 or a==d or c==b or H.has_edge(a,d) or H.has_edge(c,b): continue
  H.remove_edge(a,b); H.remove_edge(c,d); H.add_edge(a,d); H.add_edge(c,b)
  if strong and not nx.is_strongly_connected(H):
   H.remove_edge(a,d); H.remove_edge(c,b); H.add_edge(a,b); H.add_edge(c,d); continue
  nr=rp(H); ok=nr>=cur if maxi else nr<=cur
  if ok: cur=nr; ed=list(H.edges())
  else: H.remove_edge(a,d); H.remove_edge(c,b); H.add_edge(a,b); H.add_edge(c,d)
 return H
def assign(G,pool,seed):
 rg=np.random.default_rng(seed); H=G.copy(); ed=sorted(H.edges()); order=rg.permutation(len(pool))
 for k,(u,v) in enumerate(ed): H[u][v]['reliability'],H[u][v]['latency_ms']=map(float,pool[order[k%len(pool)]])
 return H
def mat(G):
 ns=list(G.nodes()); ix={n:i for i,n in enumerate(ns)}; N=len(ns); W=np.zeros((N,N)); D=np.zeros((N,N))
 for u,v,d in G.edges(data=True): W[ix[v],ix[u]]=max(float(d.get('reliability',1)),1e-9); D[ix[v],ix[u]]=float(d.get('latency_ms',0))/1000
 for i in range(N):
  s=W[i].sum()
  if s>0: W[i]*=GAIN/s
 return W,D
def bp(W,D):
 N=W.shape[0]; I=np.eye(N); one=np.ones(N); Q=np.eye(N)*M0; vals=[]
 for f in FREQS:
  h=1/(1+1j*2*np.pi*f*TAU); K=W*h*np.exp(-1j*2*np.pi*f*D); H=np.linalg.inv(I-K); S=H@Q@H.conj().T; vals.append(float(np.real(one@S@one)/(N*N)))
 vals=np.array(vals); z={}
 for name,a,b in [('slow',0,2),('theta',4,8),('gamma',30,100)]:
  m=(FREQS>=a)&(FREQS<=b); z[name]=float(np.trapezoid(vals[m],FREQS[m]))
 return z
def sim(W,D,seed,T=60):
 rg=np.random.default_rng(seed); N=len(W); steps=int(T/DT); burn=int(10/DT); db=np.clip(np.rint(D/DT).astype(int),1,25); mx=int(db.max()); ring=np.zeros((mx+1,N)); state=np.zeros(N); decay=np.exp(-DT/TAU); jump=(1-decay)/DT; target=np.full(N,M0); mu=np.clip(target-W@target,1e-8,None); pop=np.zeros(steps); tot=0
 masks=[(k,np.where(db==k,W,0)) for k in range(1,mx+1) if np.any(db==k)]
 for t in range(steps):
  slot=t%(mx+1); arr=ring[slot].copy(); ring[slot]=0; state=state*decay+jump*arr; ev=rg.poisson(np.clip(mu+state,0,1000)*DT); tot+=ev.sum(); pop[t]=ev.sum()/(N*DT)
  if ev.any():
   for k,M in masks: ring[(t+k)%(mx+1)]+=M@ev
 x=pop[burn:]; lag=10; ac=float(np.corrcoef(x[:-lag],x[lag:])[0,1]) if x.std()>0 else np.nan
 return tot/(N*T),ac
def save(G,p): pd.DataFrame([{'source':u,'target':v,'reliability':d.get('reliability',''),'latency_ms':d.get('latency_ms','')} for u,v,d in G.edges(data=True)]).to_csv(p,index=False)
def pval(x):
 x=np.asarray(x,float); x=x[np.isfinite(x)]
 try: return float(wilcoxon(x).pvalue) if len(x)>1 and not np.allclose(x,0) else np.nan
 except: return np.nan
def ci(x,B=5000,seed=SEED):
 x=np.asarray(x,float); x=x[np.isfinite(x)]; rg=np.random.default_rng(seed); n=len(x); b=np.median(x[rg.integers(0,n,(B,n))],axis=1); return [float(np.quantile(b,.025)),float(np.quantile(b,.975))]
# data
cc=out/'ccepData_V1.mat'; dl('https://data.nemar.org/on004080/v1.0.0/derivatives/av_ccep/ccepData_V1.mat',cc)
h=requests.get('https://www.nature.com/articles/s41467-025-66988-y',headers={'User-Agent':'Mozilla/5.0'},timeout=90); h.raise_for_status(); soup=BeautifulSoup(h.text,'html.parser'); src=None
for a in soup.find_all('a',href=True):
 if 'Source data file' in ' '.join(a.stripped_strings) and 'zip' in ' '.join(a.stripped_strings).lower(): src=a['href']; break
rz=requests.get(src,headers={'User-Agent':'Mozilla/5.0'},timeout=120); rz.raise_for_status(); zipfile.ZipFile(io.BytesIO(rz.content)).extractall(out/'glioma_source')
# phenotype axis
f=pd.read_excel(out/'glioma_source'/'SourceData'/'Fig. 5b'/'Fig. 5b.xlsx',header=None)
def after(label):
 rows=[]; on=False
 for _,r in f.iterrows():
  x=str(r.iloc[0]).strip() if pd.notna(r.iloc[0]) else ''
  if x.lower()==label.lower(): on=True; continue
  if on:
   v=pd.to_numeric(r.iloc[:2],errors='coerce').to_numpy(float)
   if np.all(np.isfinite(v)): rows.append(v)
   elif rows: break
 return np.array(rows)
sp,hj=after('Spontaneous'),after('Hijacked'); ev=np.log(hj.mean(0)/sp.mean(0)); en=float(np.linalg.norm(ev)); eu=ev/en
axis={'spontaneous_replicates':sp.tolist(),'hijacked_replicates':hj.tolist(),'delta_log_theta':float(ev[0]),'delta_log_gamma':float(ev[1]),'definition':'alpha=dot(model delta log power, experimental delta log power)/||experimental delta log power||^2. This is a dimensionless cross-platform spectral-direction benchmark, not a migration regression or CMC prediction.'}; (out/'phenotype_axis.json').write_text(json.dumps(axis,indent=2))
def proj(pl,ph):
 d=np.array([np.log(ph['theta']/pl['theta']),np.log(ph['gamma']/pl['gamma'])]); return float(d@ev/(ev@ev)),float(d@eu)
# reproducible synthetic ensemble
primary=[]
for k in range(24):
 rg=np.random.default_rng(SEED+1000+k)
 while True:
  U=nx.fast_gnp_random_graph(40,.15,seed=int(rg.integers(1,2**31-1)))
  if 90<=U.number_of_edges()<=155 and sum(U.degree(i)==0 for i in U)<=3: break
 order=rg.permutation(40); rank=np.empty(40,int); rank[order]=np.arange(40); lo=nx.DiGraph(); lo.add_nodes_from(range(40))
 for a,b in U.edges: lo.add_edge(a,b) if rank[a]<rank[b] else lo.add_edge(b,a)
 hi=opt(lo,True,3200,SEED+9000+k)
 for G in [lo,hi]:
  for u,v in G.edges: G[u][v]['reliability']=1.; G[u][v]['latency_ms']=0.
 save(lo,out/'graphs'/f'primary_{k+1:02d}_low.csv'); save(hi,out/'graphs'/f'primary_{k+1:02d}_high.csv'); Wl,Dl=mat(lo); Wh,Dh=mat(hi); pl,ph=bp(Wl,Dl),bp(Wh,Dh); ss=[(sim(Wl,Dl,SEED+100000+k*10+r),sim(Wh,Dh,SEED+200000+k*10+r)) for r in range(4)]; al,_=proj(pl,ph)
 primary.append({'pair':k+1,'recip_low':rp(lo),'recip_high':rp(hi),'rate_low':np.mean([q[0][0] for q in ss]),'rate_high':np.mean([q[1][0] for q in ss]),'ac100_low':np.mean([q[0][1] for q in ss]),'ac100_high':np.mean([q[1][1] for q in ss]),'slow_low':pl['slow'],'slow_high':ph['slow'],'theta_low':pl['theta'],'theta_high':ph['theta'],'gamma_low':pl['gamma'],'gamma_high':ph['gamma'],'axis_alpha':al})
pdf=pd.DataFrame(primary); pdf.to_csv(out/'primary_synthetic_results.csv',index=False)
# fixed-gain falsification
fixed=[]
for k in range(16):
 rg=np.random.default_rng(SEED+30000+k)
 while True:
  G=nx.gnm_random_graph(40,int(rg.integers(110,166)),seed=int(rg.integers(1,2**31-1)),directed=True); G.remove_edges_from(nx.selfloop_edges(G))
  if nx.is_strongly_connected(G) and min(dict(G.in_degree()).values())>0: break
 lo=opt(G,False,5000,SEED+31000+k,True); hi=opt(G,True,5000,SEED+32000+k,True)
 for H in [lo,hi]:
  for u,v in H.edges: H[u][v]['reliability']=1.; H[u][v]['latency_ms']=0.
 save(lo,out/'graphs'/f'fixed_{k+1:02d}_low.csv'); save(hi,out/'graphs'/f'fixed_{k+1:02d}_high.csv'); Wl,Dl=mat(lo); Wh,Dh=mat(hi); pl,ph=bp(Wl,Dl),bp(Wh,Dh)
 fixed.append({'pair':k+1,'recip_low':rp(lo),'recip_high':rp(hi),'rho_low':float(max(abs(np.linalg.eigvals(Wl)))),'rho_high':float(max(abs(np.linalg.eigvals(Wh)))),'slow_low':pl['slow'],'slow_high':ph['slow']})
fx=pd.DataFrame(fixed); fx.to_csv(out/'fixed_gain_falsification.csv',index=False)
# subject regional matrices
cd=loadmat(cc,simplify_cells=True)['ccepData']; cd=[cd] if isinstance(cd,dict) else list(cd); subs=[]; edges=[]; cres=[]
def roi(x):
 try:
  v=int(float(str(x).strip())); return v if v>0 else None
 except: return None
def hmap(s):
 e=s.get('electrodes',{}); m={}
 if isinstance(e,dict):
  names=np.atleast_1d(e.get('name',[])).astype(str); hv=np.atleast_1d(e.get('jsonHemi',e.get('hemisphere',[]))).astype(str)
  if len(hv)!=len(names): hv=np.array(['']*len(names))
  for n,h in zip(names,hv): m[n]='L' if 'L' in h.upper() else ('R' if 'R' in h.upper() else '')
 return m
for si,s in enumerate(cd):
 sid=str(s['id']); age=float(s['age']); hm=hmap(s); tested=defaultdict(int); det=defaultdict(int); lats=defaultdict(list); runs=[s['run']] if isinstance(s['run'],dict) else list(s['run'])
 for r in runs:
  ch=np.atleast_1d(r['channel_names']).astype(str); pa=np.atleast_1d(r['stimpair_names']).astype(str); cr=np.atleast_1d(r['channel_DestrieuxNr']).astype(str); pr=np.asarray(r['stimpair_DestrieuxNr']).astype(str); n1=np.asarray(r['n1_peak_sample'],float); tt=np.asarray(r['tt'],float).reshape(-1); good=np.atleast_1d(r['good_channels']).astype(int).reshape(-1)-1; good=good[(good>=0)&(good<len(ch))]; gs=set(ch[good])
  for j,pair in enumerate(pa):
   pe=pair.split('-')
   if len(pe)!=2 or pe[0] not in gs or pe[1] not in gs: continue
   rr=np.atleast_1d(pr[j]); a,b=roi(rr[0]),roi(rr[1])
   if a is None or b is None or a!=b: continue
   h1,h2=hm.get(pe[0],''),hm.get(pe[1],'')
   if h1 and h2 and h1!=h2: continue
   src=f'{h1 or h2 or "?"}:{a}'
   for i in good:
    tr=roi(cr[i]); th=hm.get(ch[i],'')
    if tr is None: continue
    tgt=f'{th or "?"}:{tr}'
    if src==tgt: continue
    key=(src,tgt); tested[key]+=1; v=n1[i,j]
    if np.isfinite(v):
     q=int(round(v))-1
     if 0<=q<len(tt):
      lat=float(tt[q])*1000
      if 8<=lat<=100: det[key]+=1; lats[key].append(lat)
 nodes=sorted(set(sum(([u,v] for u,v in tested),[]))); ix={n:i for i,n in enumerate(nodes)}; n=len(nodes); T=np.zeros((n,n),int); R=np.zeros((n,n),int); P=np.full((n,n),np.nan); L=np.full((n,n),np.nan); G=nx.DiGraph(); G.add_nodes_from(nodes)
 for (u,v),nt in tested.items():
  nd=det[(u,v)]; p=nd/nt; T[ix[v],ix[u]]=nt; R[ix[v],ix[u]]=nd; P[ix[v],ix[u]]=p
  if nd and lats[(u,v)]:
   lat=float(np.median(lats[(u,v)])); L[ix[v],ix[u]]=lat; G.add_edge(u,v,reliability=p,latency_ms=lat); edges.append({'subject':sid,'age':age,'source_region':u,'target_region':v,'n_tested':nt,'n_detected':nd,'n1_detection_probability':p,'median_latency_ms':lat})
 np.savez_compressed(out/'ccep_subject_matrices'/f'{sid}.npz',nodes=np.array(nodes),n_tested=T,n_detected=R,detection_probability=P,latency_ms=L,subject=sid,age=age); ntp=sum(v>0 for v in tested.values()); subs.append({'subject':sid,'age':age,'n_regions':n,'n_tested_region_pairs':ntp,'n_responsive_region_edges':G.number_of_edges(),'responsive_fraction':G.number_of_edges()/max(ntp,1),'reciprocal_pairs_observed':rp(G),'reciprocity_fraction_observed':rf(G)})
 if n<6 or G.number_of_edges()<10: continue
 Wo,Do=mat(G); po=bp(Wo,Do); lo=opt(G,False,7000,SEED+40000+si); hi=opt(G,True,7000,SEED+50000+si); pool=[(d['reliability'],d['latency_ms']) for _,_,d in G.edges(data=True)]; lo=assign(lo,pool,SEED+60000+si); hi=assign(hi,pool,SEED+60000+si); save(G,out/'graphs'/f'ccep_{sid}_observed.csv'); save(lo,out/'graphs'/f'ccep_{sid}_low.csv'); save(hi,out/'graphs'/f'ccep_{sid}_high.csv'); Wl,Dl=mat(lo); Wh,Dh=mat(hi); pl,ph=bp(Wl,Dl),bp(Wh,Dh); ss=[(sim(Wl,Dl,SEED+700000+si*10+r),sim(Wh,Dh,SEED+800000+si*10+r)) for r in range(2)]; al,sg=proj(pl,ph)
 cres.append({'subject':sid,'age':age,'n_regions':n,'edges':G.number_of_edges(),'tested_region_pairs':ntp,'observed_reciprocal_pairs':rp(G),'observed_reciprocity_fraction':rf(G),'observed_slow':po['slow'],'observed_theta':po['theta'],'observed_gamma':po['gamma'],'low_reciprocal_pairs':rp(lo),'high_reciprocal_pairs':rp(hi),'rate_low':np.mean([q[0][0] for q in ss]),'rate_high':np.mean([q[1][0] for q in ss]),'ac100_low':np.mean([q[0][1] for q in ss]),'ac100_high':np.mean([q[1][1] for q in ss]),'slow_low':pl['slow'],'slow_high':ph['slow'],'theta_low':pl['theta'],'theta_high':ph['theta'],'gamma_low':pl['gamma'],'gamma_high':ph['gamma'],'axis_alpha':al,'axis_signed':sg})
pd.DataFrame(subs).to_csv(out/'ccep_subject_summary.csv',index=False); pd.DataFrame(edges).to_csv(out/'ccep_regional_edges.csv',index=False); cdf=pd.DataFrame(cres); cdf.to_csv(out/'ccep_subject_level_hawkes.csv',index=False)
with zipfile.ZipFile(out/'ccep_subject_matrices.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in (out/'ccep_subject_matrices').glob('*.npz'): z.write(p,p.name)
def summ(df,label):
 r={'label':label,'n_subjects':len(df)}
 for m,l,h in [('rate','rate_low','rate_high'),('ac100','ac100_low','ac100_high'),('slow','slow_low','slow_high'),('theta','theta_low','theta_high'),('gamma','gamma_low','gamma_high')]:
  d=(df[h]-df[l]).to_numpy(float); r[m]={'median_low':float(df[l].median()),'median_high':float(df[h].median()),'median_difference':float(np.median(d)),'bootstrap_95CI':ci(d,seed=SEED+len(df)),'positive':int(np.sum(d>0)),'wilcoxon_p':pval(d)}
 a=df.axis_alpha.to_numpy(float); r['axis_alpha']={'median':float(np.median(a)),'bootstrap_95CI':ci(a,seed=SEED+99+len(df)),'positive':int(np.sum(a>0)),'wilcoxon_p':pval(a)}; return r
assoc={}
for m in ['observed_slow','observed_theta','observed_gamma']:
 ro,p=spearmanr(cdf.observed_reciprocity_fraction,cdf[m]); assoc[m]={'rho':float(ro),'p':float(p)}
summary={'seed':SEED,'eligible_subjects':len(cdf),'adult_eligible_subjects':int((cdf.age>=18).sum()),'all_subjects':summ(cdf,'all'),'adult_age_ge_18':summ(cdf[cdf.age>=18],'adult'),'observed_reciprocity_associations':assoc,'phenotype_axis':axis,'primary_synthetic':{'n_pairs':24,'recip_diff_median':float(np.median(pdf.recip_high-pdf.recip_low)),'rate_diff_median':float(np.median(pdf.rate_high-pdf.rate_low)),'rate_p':pval(pdf.rate_high-pdf.rate_low),'ac100_diff_median':float(np.median(pdf.ac100_high-pdf.ac100_low)),'ac100_p':pval(pdf.ac100_high-pdf.ac100_low),'slow_diff_median':float(np.median(pdf.slow_high-pdf.slow_low)),'slow_p':pval(pdf.slow_high-pdf.slow_low),'axis_alpha_median':float(np.median(pdf.axis_alpha)),'axis_alpha_p':pval(pdf.axis_alpha)},'fixed_gain_control':{'n_pairs':16,'recip_diff_median':float(np.median(fx.recip_high-fx.recip_low)),'rho_low_median':float(fx.rho_low.median()),'rho_high_median':float(fx.rho_high.median()),'slow_diff_median':float(np.median(fx.slow_high-fx.slow_low)),'slow_p':pval(fx.slow_high-fx.slow_low)}}
(out/'v9_summary.json').write_text(json.dumps(summary,indent=2)); (out/'requirements.txt').write_text('numpy\npandas\nscipy\nnetworkx\nrequests\nbeautifulsoup4\nopenpyxl\n'); (out/'PROVENANCE.txt').write_text('CCEP: NEMAR on004080 ccepData_V1.mat. Regional edges retain tested nonresponses. Coupling is N1 detection probability, not voltage amplitude. Glioma benchmark is a theta/gamma phenotype-axis projection, not a fitted migration receiver. Seed 20261008.\n'); (out/'README_V9.txt').write_text(f'Frozen reproducibility package. Eligible subject-level CCEP networks: {len(cdf)}/74; adult eligible: {(cdf.age>=18).sum()}. Includes exact graph edge lists, subject matrices, analysis outputs, requirements, provenance, and reproduce_v9.py.\n'); (out/'reproduce_v9.py').write_text(Path(__file__).read_text()); print(json.dumps(summary,indent=2))
