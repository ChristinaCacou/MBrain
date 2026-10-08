import io, json, math, os, zipfile, hashlib
from pathlib import Path
from collections import defaultdict, Counter
import requests
import numpy as np
import pandas as pd
import networkx as nx
from scipy.io import loadmat
from scipy.stats import wilcoxon, spearmanr

SEED=20261008
rng=np.random.default_rng(SEED)
out=Path('empirical_results'); out.mkdir(exist_ok=True)

# ----------------------- data retrieval -----------------------
base='https://data.nemar.org/on004080/v1.0.0/derivatives/av_ccep/'
def dl(url,path):
    r=requests.get(url,timeout=180,allow_redirects=True,headers={'User-Agent':'Mozilla/5.0'})
    r.raise_for_status(); Path(path).write_bytes(r.content); return len(r.content)
cc_path=out/'ccepData_V1.mat'
dl(base+'ccepData_V1.mat',cc_path)
article='https://www.nature.com/articles/s41467-025-66988-y'
h=requests.get(article,timeout=90,headers={'User-Agent':'Mozilla/5.0'}); h.raise_for_status()
from bs4 import BeautifulSoup
soup=BeautifulSoup(h.text,'html.parser')
src=None
for a in soup.find_all('a',href=True):
    txt=' '.join(a.stripped_strings)
    if 'Source data file' in txt and 'zip' in txt.lower(): src=a['href']; break
if not src: raise RuntimeError('Source data ZIP link not found')
if src.startswith('/'): src='https://www.nature.com'+src
rz=requests.get(src,timeout=120,headers={'User-Agent':'Mozilla/5.0'}); rz.raise_for_status()
gzip=out/'glioma_source_data.zip'; gzip.write_bytes(rz.content)
with zipfile.ZipFile(io.BytesIO(rz.content)) as z: z.extractall(out/'glioma_source')

# ----------------------- CCEP subject matrices -----------------------
cd=loadmat(cc_path,simplify_cells=True)['ccepData']
if isinstance(cd,dict): cd=[cd]
subjects=[]; edge_rows=[]
matrix_dir=out/'ccep_subject_matrices'; matrix_dir.mkdir(exist_ok=True)

for s in cd:
    sid=str(s['id']); age=float(s['age'])
    runs=s['run']; runs=[runs] if isinstance(runs,dict) else list(runs)
    all_pairs=sorted({str(x) for r in runs for x in np.atleast_1d(r['stimpair_names'])})
    all_chans=sorted({str(x) for r in runs for x in np.atleast_1d(r['channel_names'])})
    pi={x:i for i,x in enumerate(all_pairs)}; ci={x:i for i,x in enumerate(all_chans)}
    tested=np.zeros((len(all_chans),len(all_pairs)),dtype=np.int16)
    detected=np.zeros_like(tested)
    latency_lists=[[[] for _ in all_pairs] for __ in all_chans]
    for r in runs:
        chans=np.atleast_1d(r['channel_names']).astype(str)
        pairs=np.atleast_1d(r['stimpair_names']).astype(str)
        n1=np.asarray(r['n1_peak_sample'],float)
        tt=np.asarray(r['tt'],float).reshape(-1)
        good=np.atleast_1d(r['good_channels']).astype(int).reshape(-1)-1
        good=good[(good>=0)&(good<len(chans))]
        good_names=set(chans[good].tolist())
        for pj,pair in enumerate(pairs):
            pe=pair.split('-')
            stim_good=(len(pe)==2 and pe[0] in good_names and pe[1] in good_names)
            if not stim_good: continue
            for ri in good:
                ch=chans[ri]; ii=ci[ch]; jj=pi[pair]
                tested[ii,jj]+=1
                v=n1[ri,pj]
                if np.isfinite(v):
                    idx=int(round(v))-1
                    if 0<=idx<len(tt):
                        lat=float(tt[idx])*1000.0
                        if 8.0 <= lat <= 100.0:
                            detected[ii,jj]+=1; latency_lists[ii][jj].append(lat)
    strength=np.divide(detected,tested,out=np.full(tested.shape,np.nan,float),where=tested>0)
    latency=np.full(tested.shape,np.nan,float)
    for i in range(len(all_chans)):
        for j in range(len(all_pairs)):
            if latency_lists[i][j]: latency[i,j]=float(np.median(latency_lists[i][j]))
    np.savez_compressed(matrix_dir/f'{sid}.npz',
        latency_ms=latency,response_strength=strength,n_tested=tested,n_detected=detected,
        response_channels=np.array(all_chans),stim_pairs=np.array(all_pairs),subject=sid,age=age)
    pos=np.isfinite(latency) & np.isfinite(strength) & (strength>0)
    ntested=int((tested>0).sum()); nresp=int(pos.sum())
    vals_lat=latency[pos]; vals_str=strength[pos]
    subjects.append(dict(subject=sid,age=age,n_channels=len(all_chans),n_stim_pairs=len(all_pairs),
                         n_tested_pairs=ntested,n_responsive_pairs=nresp,response_fraction=nresp/max(ntested,1),
                         median_latency_ms=float(np.median(vals_lat)) if nresp else np.nan,
                         q1_latency_ms=float(np.quantile(vals_lat,.25)) if nresp else np.nan,
                         q3_latency_ms=float(np.quantile(vals_lat,.75)) if nresp else np.nan,
                         median_response_strength=float(np.median(vals_str)) if nresp else np.nan))
    for i,j in zip(*np.where(pos)):
        edge_rows.append(dict(subject=sid,age=age,stim_pair=all_pairs[j],response_channel=all_chans[i],
                              latency_ms=float(latency[i,j]),response_strength=float(strength[i,j]),
                              n_tested=int(tested[i,j]),n_detected=int(detected[i,j])))

summ=pd.DataFrame(subjects); edges=pd.DataFrame(edge_rows)
summ.to_csv(out/'ccep_subject_summary.csv',index=False)
edges.to_csv(out/'ccep_responsive_edges.csv',index=False)
with zipfile.ZipFile(out/'ccep_subject_matrices.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(matrix_dir.glob('*.npz')): z.write(p,p.name)

eligible={sid:g[['latency_ms','response_strength']].to_numpy(float) for sid,g in edges.groupby('subject') if len(g)>=20}
if len(eligible)<20: raise RuntimeError(f'Only {len(eligible)} subjects with >=20 responsive edges')
elig_ids=sorted(eligible)

# ----------------------- phenotype-linked spectral receiver -----------------------
sroot=out/'glioma_source'/'SourceData'
f5b=pd.read_excel(sroot/'Fig. 5b'/'Fig. 5b.xlsx',header=None)
def after_label(df,label):
    rows=[]; active=False
    for _,r in df.iterrows():
        first=str(r.iloc[0]).strip() if pd.notna(r.iloc[0]) else ''
        if first.lower()==label.lower(): active=True; continue
        if active:
            vals=pd.to_numeric(r.iloc[:2],errors='coerce').to_numpy(float)
            if np.all(np.isfinite(vals)): rows.append(vals)
            elif rows: break
    return np.array(rows,float)
sp_spec=after_label(f5b,'Spontaneous'); hj_spec=after_label(f5b,'Hijacked')
if sp_spec.shape!=(3,2) or hj_spec.shape!=(3,2): raise RuntimeError((sp_spec.shape,hj_spec.shape,f5b.to_string()))
sp_mean=sp_spec.mean(axis=0); hj_mean=hj_spec.mean(axis=0)
delta_emp=np.log(hj_mean/sp_mean)

f4=pd.read_excel(sroot/'Fig. 4l'/'Fig. 4l.xlsx',header=None)
f5h=pd.read_excel(sroot/'Fig. 5h'/'Fig. 5h.xlsx',header=None)
rows4=pd.to_numeric(f4.iloc[:,0],errors='coerce').notna(); d4=f4.loc[rows4].copy()
rows5=pd.to_numeric(f5h.iloc[:,0],errors='coerce').notna(); d5=f5h.loc[rows5].copy()
days=pd.to_numeric(d4.iloc[:,0],errors='coerce').to_numpy(float)
sp_cmc=pd.to_numeric(d4.iloc[:,3],errors='coerce').to_numpy(float)
sp_sd=pd.to_numeric(d4.iloc[:,4],errors='coerce').to_numpy(float)
hj_cmc=pd.to_numeric(d5.iloc[:,5],errors='coerce').to_numpy(float)
hj_sd=pd.to_numeric(d5.iloc[:,6],errors='coerce').to_numpy(float)
auc_sp=float(np.trapezoid(sp_cmc,days)); auc_hj=float(np.trapezoid(hj_cmc,days)); delta_auc=auc_hj-auc_sp
beta=delta_auc*delta_emp/float(delta_emp@delta_emp)
rep_delta=np.log(hj_spec/sp_spec)
rep_pred=np.array([float(beta@x) for x in rep_delta])
receiver={
 'definition':'Delta CMC AUC (days 1-7) = beta_theta * Delta log(theta power) + beta_gamma * Delta log(gamma power)',
 'theta_spontaneous_mean':float(sp_mean[0]),'gamma_spontaneous_mean':float(sp_mean[1]),
 'theta_hijacked_mean':float(hj_mean[0]),'gamma_hijacked_mean':float(hj_mean[1]),
 'empirical_delta_log_theta':float(delta_emp[0]),'empirical_delta_log_gamma':float(delta_emp[1]),
 'spontaneous_cmc_auc':auc_sp,'hijacked_cmc_auc':auc_hj,'empirical_delta_cmc_auc':delta_auc,
 'beta_theta':float(beta[0]),'beta_gamma':float(beta[1]),
 'paired_spectral_replicate_implied_delta_auc':rep_pred.tolist(),
 'calibration_scope':'two-condition aggregate calibration (spontaneous vs hijacked); not an independently validated cell-level predictor'
}
(out/'tumor_receiver_fit.json').write_text(json.dumps(receiver,indent=2))

# ----------------------- empirical-CCEP Hawkes validation -----------------------
N=40; GAIN=.85; M0=5.0; TAU=.050; DT=.010

def recips(G): return sum(1 for u,v in G.edges if G.has_edge(v,u))//2

def make_pair(seed):
    rg=np.random.default_rng(seed)
    while True:
        ug=nx.fast_gnp_random_graph(N,0.15,seed=int(rg.integers(1,2**31-1)),directed=False)
        if 90<=ug.number_of_edges()<=155 and sum(1 for x in ug if ug.degree(x)==0)<=3: break
    order=rg.permutation(N); rank=np.empty(N,int); rank[order]=np.arange(N)
    G=nx.DiGraph(); G.add_nodes_from(range(N))
    for a,b in ug.edges:
        G.add_edge(a,b) if rank[a]<rank[b] else G.add_edge(b,a)
    H=G.copy(); cur=recips(H); ed=list(H.edges())
    for _ in range(3200):
        if len(ed)<2: break
        i,j=rg.choice(len(ed),2,replace=False); a,b=ed[i]; c,d=ed[j]
        if len({a,b,c,d})<4: continue
        ne1=(a,d); ne2=(c,b)
        if a==d or c==b or H.has_edge(*ne1) or H.has_edge(*ne2): continue
        H.remove_edge(a,b); H.remove_edge(c,d); H.add_edge(*ne1); H.add_edge(*ne2)
        nr=recips(H)
        if nr>=cur:
            cur=nr; ed=list(H.edges())
        else:
            H.remove_edge(*ne1); H.remove_edge(*ne2); H.add_edge(a,b); H.add_edge(c,d)
    return G,H

def field_for_pair(seed):
    rg=np.random.default_rng(seed)
    d=np.zeros((N,N),float); strength=np.zeros((N,N),float); contrib=[]
    for src in range(N):
        for tgt in range(N):
            if src==tgt: continue
            sid=elig_ids[int(rg.integers(len(elig_ids)))]; pool=eligible[sid]
            row=pool[int(rg.integers(len(pool)))]
            d[tgt,src]=row[0]/1000.0; strength[tgt,src]=max(row[1],1e-6); contrib.append(sid)
    return d,strength,contrib

def weighted(G,S):
    W=np.zeros((N,N),float)
    for src,tgt in G.edges: W[tgt,src]=S[tgt,src]
    for i in range(N):
        sm=W[i].sum()
        if sm>0: W[i]*=GAIN/sm
    return W

def bandpowers(W,D):
    freqs=np.arange(0,100.0001,0.5); one=np.ones(N); vals=np.empty_like(freqs)
    I=np.eye(N); diag=np.eye(N)*M0
    for k,f in enumerate(freqs):
        h=1/(1+1j*2*np.pi*f*TAU)
        K=W*h*np.exp(-1j*2*np.pi*f*D)
        H=np.linalg.inv(I-K); SS=H@diag@H.conj().T
        vals[k]=float(np.real(one@SS@one)/(N*N))
    ret={}
    for name,a,b in [('slow',0,2),('theta',4,8),('gamma',30,100)]:
        mask=(freqs>=a)&(freqs<=b); ret[name]=float(np.trapezoid(vals[mask],freqs[mask]))
    return ret

def simulate(W,D,seed,T=60.0):
    rg=np.random.default_rng(seed); steps=int(T/DT); burn=int(10/DT)
    decay=math.exp(-DT/TAU); delay_bins=np.clip(np.rint(D/DT).astype(int),1,20)
    maxd=int(delay_bins.max()); ring=np.zeros((maxd+1,N),float); state=np.zeros(N,float)
    mu=np.full(N,M0)-W@np.full(N,M0); mu=np.clip(mu,1e-6,None)
    pop=np.zeros(steps,float); node_sum=np.zeros(N,float)
    masks=[]
    for db in range(1,maxd+1):
        M=np.where(delay_bins==db,W,0.0)
        if np.any(M): masks.append((db,M))
    for t in range(steps):
        slot=t%(maxd+1); arrivals=ring[slot].copy(); ring[slot]=0
        state=state*decay+arrivals/TAU
        lam=np.clip(mu+state,0,500)
        ev=rg.poisson(lam*DT)
        node_sum+=ev; pop[t]=ev.sum()/(N*DT)
        if ev.any():
            for db,M in masks:
                target=(t+db)%(maxd+1); ring[target]+=M@ev
    x=pop[burn:]
    x=np.convolve(x,np.ones(2)/2,mode='valid')
    lag=int(round(.1/DT)); a=x[:-lag]; b=x[lag:]
    ac=float(np.corrcoef(a,b)[0,1]) if a.std()>0 and b.std()>0 else np.nan
    return float(node_sum.sum()/(N*T)),ac

rows=[]; contributors=Counter()
for k in range(24):
    Glo,Ghi=make_pair(SEED+1000+k)
    D,S,cc=field_for_pair(SEED+5000+k); contributors.update(cc)
    Wlo=weighted(Glo,S); Whi=weighted(Ghi,S)
    plo=bandpowers(Wlo,D); phi=bandpowers(Whi,D)
    sim=[]
    for rep in range(3):
        sim.append((simulate(Wlo,D,SEED+100000+k*20+rep),simulate(Whi,D,SEED+200000+k*20+rep)))
    rate_lo=np.mean([z[0][0] for z in sim]); rate_hi=np.mean([z[1][0] for z in sim])
    ac_lo=np.mean([z[0][1] for z in sim]); ac_hi=np.mean([z[1][1] for z in sim])
    dlog=np.array([math.log(phi['theta']/plo['theta']),math.log(phi['gamma']/plo['gamma'])])
    pred=float(beta@dlog)
    rows.append(dict(pair=k+1,low_reciprocal_pairs=recips(Glo),high_reciprocal_pairs=recips(Ghi),
                     mean_rate_low=rate_lo,mean_rate_high=rate_hi,autocorr100_low=ac_lo,autocorr100_high=ac_hi,
                     slow_low=plo['slow'],slow_high=phi['slow'],theta_low=plo['theta'],theta_high=phi['theta'],
                     gamma_low=plo['gamma'],gamma_high=phi['gamma'],
                     delta_log_theta=float(dlog[0]),delta_log_gamma=float(dlog[1]),predicted_delta_cmc_auc=pred))
res=pd.DataFrame(rows); res.to_csv(out/'empirical_ccep_hawkes_results.csv',index=False)

def paired_stats(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); d=b-a
    try: w=wilcoxon(d,alternative='two-sided').pvalue
    except Exception: w=np.nan
    return {'low_median':float(np.median(a)),'high_median':float(np.median(b)),'median_difference':float(np.median(d)),
            'mean_difference':float(np.mean(d)),'positive_pairs':int((d>0).sum()),'n':int(len(d)),'wilcoxon_p':float(w)}
stats={
 'ccep':{
   'subjects_total':int(len(summ)),'subjects_with_responsive_edges':int(edges.subject.nunique()),
   'responsive_edge_entries':int(len(edges)),'pooled_latency_median_ms':float(edges.latency_ms.median()),
   'pooled_latency_q1_ms':float(edges.latency_ms.quantile(.25)),'pooled_latency_q3_ms':float(edges.latency_ms.quantile(.75)),
   'subject_median_latency_median_ms':float(summ.median_latency_ms.median()),
   'subject_response_fraction_median':float(summ.response_fraction.median()),
   'subjects_contributing_to_hawkes_fields':int(sum(v>0 for v in contributors.values()))
 },
 'receiver':receiver,
 'hawkes':{
   'mean_rate':paired_stats(res.mean_rate_low,res.mean_rate_high),
   'autocorr100':paired_stats(res.autocorr100_low,res.autocorr100_high),
   'slow_power':paired_stats(res.slow_low,res.slow_high),
   'theta_power':paired_stats(res.theta_low,res.theta_high),
   'gamma_power':paired_stats(res.gamma_low,res.gamma_high),
   'predicted_delta_cmc_auc':{
      'median':float(res.predicted_delta_cmc_auc.median()),'mean':float(res.predicted_delta_cmc_auc.mean()),
      'positive_pairs':int((res.predicted_delta_cmc_auc>0).sum()),'n':24,
      'wilcoxon_vs_zero_p':float(wilcoxon(res.predicted_delta_cmc_auc).pvalue)
   }
 },
 'analysis_note':'CCEP latency is measured N1 timing. Response strength is the within-subject fraction of eligible runs in which the official N1 detector identified a response for each stimulation-pair/recording-channel pair. The glioma receiver is a two-condition aggregate spectral calibration, not independent patient-level validation. Hawkes graph pairs are an independent validation ensemble regenerated from the manuscript-described graph constraints because the original edge-list reproducibility package was not supplied.'
}
(out/'empirical_results_summary.json').write_text(json.dumps(stats,indent=2))
text=[]
text.append(f"Actual CCEP matrices were reconstructed for {stats['ccep']['subjects_total']} subjects; {stats['ccep']['responsive_edge_entries']:,} responsive stimulation-response entries were retained. Pooled N1 latency was {stats['ccep']['pooled_latency_median_ms']:.1f} ms (IQR {stats['ccep']['pooled_latency_q1_ms']:.1f}-{stats['ccep']['pooled_latency_q3_ms']:.1f} ms).")
for key,label in [('autocorr100','100-ms autocorrelation'),('slow_power','0-2 Hz power'),('theta_power','theta power'),('gamma_power','gamma power')]:
    q=stats['hawkes'][key]; text.append(f"{label}: median high-minus-low difference {q['median_difference']:.5g}; {q['positive_pairs']}/{q['n']} positive; Wilcoxon P={q['wilcoxon_p']:.3g}.")
q=stats['hawkes']['predicted_delta_cmc_auc']; text.append(f"Phenotype-linked receiver: empirical spontaneous-to-hijacked CMC AUC difference {receiver['empirical_delta_cmc_auc']:.1f}; beta_theta={receiver['beta_theta']:.1f}, beta_gamma={receiver['beta_gamma']:.1f}. Applied to Hawkes theta/gamma shifts, median predicted high-minus-low CMC AUC change was {q['median']:.1f} ({q['positive_pairs']}/24 positive; P={q['wilcoxon_vs_zero_p']:.3g}).")
(out/'manuscript_ready_results.txt').write_text('\n'.join(text))

import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(6.4,4.2)); ax.hist(edges.latency_ms,bins=np.arange(8,102,2)); ax.set_xlabel('Measured N1 latency (ms)'); ax.set_ylabel('Responsive CCEP entries'); ax.set_title('Subject-level CCEP N1 latency distribution'); fig.tight_layout(); fig.savefig(out/'fig_ccep_latency.png',dpi=220); plt.close(fig)
fig,ax=plt.subplots(figsize=(6.4,4.2)); ax.scatter(res.slow_high/res.slow_low, res.predicted_delta_cmc_auc); ax.axhline(0,lw=1); ax.set_xlabel('High / low 0-2 Hz power'); ax.set_ylabel('Predicted Delta CMC AUC'); ax.set_title('Empirically calibrated tumor receiver'); fig.tight_layout(); fig.savefig(out/'fig_receiver_validation.png',dpi=220); plt.close(fig)

print(json.dumps(stats,indent=2))
