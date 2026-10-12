import json, os, re, urllib.request
from pathlib import Path
import numpy as np, h5py
from scipy.stats import spearmanr

SEED=20261012
ASSET=os.environ['ASSET_ID']
OUT=Path(os.environ.get('OUTDIR','analysis/chen2024_sessions')); OUT.mkdir(parents=True,exist_ok=True)
FP=Path('chen2024_session.nwb')

def dec(x):
    if isinstance(x,bytes):return x.decode(errors='replace')
    if isinstance(x,np.generic):return x.item()
    return x

def arrstr(ds): return np.array([str(dec(x)) for x in np.asarray(ds[()]).ravel()],dtype=object)
def parse_json_string(x):
    try:return json.loads(str(dec(x)))
    except:return None

def side_from_opto_location(raw):
    z=parse_json_string(raw)
    if not isinstance(z,list):return 'unknown'
    sides=set()
    for q in z:
        try:ml=float(q.get('ml_location'))
        except:continue
        if ml<0:sides.add('left')
        elif ml>0:sides.add('right')
    if sides=={'left'}:return 'left'
    if sides=={'right'}:return 'right'
    if sides=={'left','right'}:return 'bilateral'
    return 'unknown'

def get_subject(f):
    for q in ['general/subject/description','general/subject/subject_id']:
        if q in f:
            s=str(dec(f[q][()]))
            if s:return s
    return 'unknown'

def map_stim_trials(f,start,stop):
    n=len(start); lab=np.array(['none']*n,dtype=object); detail=[]
    base='acquisition/BehavioralEvents/photostim_start_times'
    if base+'/timestamps' not in f:return lab,detail
    ts=np.asarray(f[base+'/timestamps'][()],float).ravel()
    ctl=np.asarray(f[base+'/control'][()]).ravel() if base+'/control' in f else np.full(len(ts),-1)
    code_side={}
    if 'general/optogenetics' in f:
        for name,g in f['general/optogenetics'].items():
            m=re.search(r'_(\d+)$',name)
            if not m or 'location' not in g:continue
            code_side[int(m.group(1))]=side_from_opto_location(g['location'][()])
    for t,c in zip(ts,ctl):
        i=int(np.searchsorted(start,t,side='right')-1)
        if i<0 or i>=n or t>stop[i]:continue
        s=code_side.get(int(c),'unknown')
        if lab[i]=='none':lab[i]=s
        elif lab[i]!=s:lab[i]='mixed'
        detail.append({'trial_index':i,'timestamp':float(t),'control_code':int(c),'side':s})
    return lab,detail

def unit_sides_and_spikes(f):
    q=arrstr(f['units/unit_quality']) if 'units/unit_quality' in f else np.array(['']*len(f['units/id']),object)
    cl=arrstr(f['units/classification']) if 'units/classification' in f else np.array(['']*len(q),object)
    good=np.array([(a.lower()=='good' or b.lower()=='good') for a,b in zip(q,cl)])
    eloc=arrstr(f['general/extracellular_ephys/electrodes/location'])
    flat=np.asarray(f['units/electrodes'][()],int).ravel(); eind=np.asarray(f['units/electrodes_index'][()],int).ravel() if 'units/electrodes_index' in f else np.arange(1,len(flat)+1)
    st=np.asarray(f['units/spike_times'][()],float).ravel(); sind=np.asarray(f['units/spike_times_index'][()],int).ravel()
    sides=[]; spikes=[]; a=0; ea=0
    for i,b in enumerate(sind):
        eb=eind[i] if i<len(eind) else ea+1
        eis=flat[ea:eb]; ea=eb
        side='other'
        for ei in eis:
            if 0<=ei<len(eloc):
                z=parse_json_string(eloc[ei]); br=str((z or {}).get('brain_regions','')).lower() if isinstance(z,dict) else ''
                if 'left alm' in br:side='left';break
                if 'right alm' in br:side='right';break
        spikes.append(st[a:b]);a=b
        sides.append(side if good[i] else 'bad')
    return np.array(sides,object),spikes

def count_matrix(spikes,idx,go):
    X=np.zeros((len(idx),len(spikes)),float)
    for j,s in enumerate(spikes):
        lo=go[idx]-0.4; hi=go[idx]
        X[:,j]=np.searchsorted(s,hi,side='left')-np.searchsorted(s,lo,side='left')
    return X

def cd(X,y):
    if len(X)==0 or np.sum(y==1)<1 or np.sum(y==0)<1:return None
    w=X[y==1].mean(0)-X[y==0].mean(0); z=np.linalg.norm(w)
    return None if not np.isfinite(z) or z<=0 else w/z

def coupling(XL,XR,y):
    rr=np.where(y==1)[0];ll=np.where(y==0)[0]
    if len(rr)<10 or len(ll)<10:return None,None,0
    rng=np.random.default_rng(SEED); vals=[]; within=[]
    for _ in range(30):
        r=rr.copy();l=ll.copy();rng.shuffle(r);rng.shuffle(l);nr=len(r)//2;nl=len(l)//2
        tr=np.r_[r[:nr],l[:nl]];te=np.r_[r[nr:],l[nl:]]
        wL=cd(XL[tr],y[tr]);wR=cd(XR[tr],y[tr])
        if wL is None or wR is None or len(te)<6:continue
        pL=XL[te]@wL;pR=XR[te]@wR;rho=spearmanr(pL,pR).statistic
        if np.isfinite(rho):vals.append(float(rho))
        vv=[]
        for z in [0,1]:
            k=y[te]==z
            if k.sum()>=3:
                q=spearmanr(pL[k],pR[k]).statistic
                if np.isfinite(q):vv.append(float(q))
        if vv:within.append(float(np.mean(vv)))
    return (float(np.mean(vals)) if vals else None),(float(np.mean(within)) if within else None),len(vals)

def bacc(mask,instr,outcome):
    v=[]
    for s in ['left','right']:
        k=mask&(instr==s)
        if k.sum()<5:return None
        v.append(float(np.mean(outcome[k]=='hit')))
    return float(np.mean(v))
def sep(sc,mask,instr):
    a=mask&(instr=='right');b=mask&(instr=='left')
    if a.sum()<5 or b.sum()<5:return None
    return float(np.mean(sc[a])-np.mean(sc[b]))

def main():
    urllib.request.urlretrieve(f'https://api.dandiarchive.org/api/assets/{ASSET}/download/',FP)
    try:
        with h5py.File(FP,'r') as f:
            subject=get_subject(f); identifier=str(dec(f['identifier'][()])) if 'identifier' in f else ASSET
            tr=f['intervals/trials'];start=np.asarray(tr['start_time'][()],float);stop=np.asarray(tr['stop_time'][()],float);n=len(start)
            instr=arrstr(tr['trial_instruction']);outcome=arrstr(tr['outcome']);early=arrstr(tr['early_lick']) if 'early_lick' in tr else np.array(['no early']*n,object)
            auto=np.asarray(tr['auto_water'][()],int) if 'auto_water' in tr else np.zeros(n,int);free=np.asarray(tr['free_water'][()],int) if 'free_water' in tr else np.zeros(n,int)
            go=np.asarray(f['acquisition/BehavioralEvents/go_start_times/timestamps'][()],float).ravel()[:n]
            stim,stimdetail=map_stim_trials(f,start,stop)
            valid=(np.isin(outcome,['hit','miss']))&(early=='no early')&(auto==0)&(free==0)&np.isin(instr,['left','right'])
            ctrl=valid&(stim=='none');correct=ctrl&(outcome=='hit')
            sides,spikes=unit_sides_and_spikes(f); li=np.where(sides=='left')[0];ri=np.where(sides=='right')[0]
            idx=np.arange(n);XL=count_matrix([spikes[i] for i in li],idx,go) if len(li) else np.empty((n,0));XR=count_matrix([spikes[i] for i in ri],idx,go) if len(ri) else np.empty((n,0))
            cc=np.where(correct)[0];y=(instr[cc]=='right').astype(int)
            coup=within=None;nrep=0
            if len(li)>5 and len(ri)>5:coup,within,nrep=coupling(XL[cc],XR[cc],y)
            ba0=bacc(ctrl,instr,outcome); rec={}
            for s in ['left','right']:
                bs=bacc(valid&(stim==s),instr,outcome);rec[s]=None if ba0 is None or bs is None or ba0<=0 else float(bs/ba0)
            rec['side_average']=float(np.mean([rec['left'],rec['right']])) if rec['left'] is not None and rec['right'] is not None else None
            neuronal={}
            for s,X in [('left',XL),('right',XR)]:
                w=cd(X[cc],y) if X.shape[1]>5 else None
                if w is None:neuronal[s]=None;continue
                sc=X@w;q0=sep(sc,ctrl,instr);qs=sep(sc,valid&(stim==s),instr)
                neuronal[s]=None if q0 is None or qs is None or abs(q0)<1e-12 else float(qs/q0)
            res={'asset_id':ASSET,'identifier':identifier,'subject':subject,'n_trials':n,'n_good_units_left_ALM':int(len(li)),'n_good_units_right_ALM':int(len(ri)),'n_correct_control_left':int(np.sum(correct&(instr=='left'))),'n_correct_control_right':int(np.sum(correct&(instr=='right'))),'stim_trial_counts':{s:int(np.sum(valid&(stim==s))) for s in ['left','right','bilateral','unknown','mixed']},'control_trial_count':int(ctrl.sum()),'control_bacc':ba0,'control_bilateral_coupling':coup,'control_bilateral_within_choice_coupling':within,'n_split_reps':nrep,'behavior_recovery':rec,'neuronal_recovery':neuronal,'stim_code_side_examples':stimdetail[:20],'lock':'external_replication/CHEN2024_EXTERNAL_LOCK.md'}
    finally: FP.unlink(missing_ok=True)
    (OUT/f'{ASSET}.json').write_text(json.dumps(res,indent=2,default=str));print(json.dumps(res,indent=2,default=str))
if __name__=='__main__':main()
