import os,re,json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.stats import rankdata
import mat73

SEED=20261012
OUT=Path('analysis/alm3_sessions'); OUT.mkdir(parents=True,exist_ok=True)
PATH=Path(os.environ.get('LOCAL_FILE','session.mat'))
KEY=os.environ.get('S3_KEY','')

def fld(x,n,d=None): return x.get(n,d) if isinstance(x,dict) else getattr(x,n,d)
def s1(x):
    if x is None:return ''
    if isinstance(x,bytes):return x.decode(errors='replace')
    if isinstance(x,str):return x
    a=np.atleast_1d(x).ravel()
    return '' if len(a)==0 else s1(a[0]) if a[0] is not x else str(x)
def alist(x):
    if x is None:return []
    if isinstance(x,list):return x
    if isinstance(x,tuple):return list(x)
    if isinstance(x,np.ndarray):return list(x.ravel())
    return [x]
def vec(x,n=None,dtype=float):
    a=np.atleast_1d(np.asarray(x).squeeze()).ravel()
    try:a=a.astype(dtype)
    except Exception:a=np.array([dtype(np.asarray(z).squeeze()) for z in a],dtype=dtype)
    if n is not None and len(a)<n:
        b=np.full(n,np.nan if dtype is float else 0,dtype=dtype); b[:len(a)]=a; a=b
    return a
def hval(h,key):
    if isinstance(h,list):
        for z in h:
            if s1(fld(z,'keyNames','')).lower()==key.lower():return fld(z,'value')
        raise KeyError(key)
    ks=[s1(z) for z in alist(fld(h,'keyNames',[]))]; vs=alist(fld(h,'value',[]))
    for i,k in enumerate(ks):
        if k.lower()==key.lower():return vs[i]
    raise KeyError(key)
def hfirst(h,names,n,dtype=float):
    for k in names:
        try:return vec(hval(h,k),n,dtype)
        except Exception:pass
    raise KeyError(names)
def masks(obj,n):
    names=[s1(z) for z in alist(fld(obj,'trialTypeStr',[]))]; M=np.asarray(fld(obj,'trialTypeMat'))
    if M.ndim==1:M=M.reshape(1,-1)
    if M.shape[0]!=len(names) and M.ndim==2 and M.shape[1]==len(names):M=M.T
    def r(k):
        for i,nm in enumerate(names):
            if nm.lower()==k.lower():
                x=np.asarray(M[i]).squeeze().astype(float); return np.pad(x,(0,max(0,n-len(x))))[:n]>0
        return np.zeros(n,bool)
    return {k:r(k) for k in ['HitR','HitL','ErrR','ErrL','LickEarly']}
def units(obj):
    H=fld(obj,'eventSeriesHash')
    U=[fld(z,'value') for z in H if fld(z,'value') is not None] if isinstance(H,list) else alist(fld(H,'value',[]))
    out=[]
    for u in U:
        c=s1(fld(u,'cellType','')).lower()
        if 'pyramidal' not in c or 'fs' in c:continue
        et=fld(u,'eventTrial',None)
        if et is None:et=fld(u,'eventTrials',None)
        if et is not None and len(np.atleast_1d(et)):out.append(u)
    return out
def matrix(U,idx,start,cue):
    idx=np.asarray(idx,int); X=np.zeros((len(idx),len(U)),float)
    for j,u in enumerate(U):
        ev=np.asarray(fld(u,'eventTimes',[]),float).ravel(); et=fld(u,'eventTrial',None)
        if et is None:et=fld(u,'eventTrials',[])
        et=np.asarray(et,int).ravel();
        for i,t0 in enumerate(idx):
            tid=t0+1; lo=start[t0]+cue[t0]-0.4; hi=start[t0]+cue[t0]
            sp=ev[et==tid]; X[i,j]=np.sum((sp>=lo)&(sp<hi))
    return X
def fitcd(X,y):
    d=X[y==1].mean(0)-X[y==0].mean(0); z=np.linalg.norm(d)
    return None if not np.isfinite(z) or z==0 else d/z
def auc(sc,y):
    a=sc[y==1]; b=sc[y==0]
    if len(a)<1 or len(b)<1:return np.nan
    z=np.r_[a,b]; r=rankdata(z); n1=len(a); n0=len(b); u=r[:n1].sum()-n1*(n1+1)/2
    q=float(u/(n1*n0)); return max(q,1-q)
def cvauc(X,y):
    n1=(y==1).sum(); n0=(y==0).sum(); k=int(min(5,n1,n0))
    if k<2:return np.nan
    rng=np.random.default_rng(SEED); split=[]
    for lab in [0,1]:
        q=np.where(y==lab)[0].copy(); rng.shuffle(q); split.append(np.array_split(q,k))
    sc=np.full(len(y),np.nan)
    for f in range(k):
        te=np.r_[split[0][f],split[1][f]]; tr=np.setdiff1d(np.arange(len(y)),te)
        cd=fitcd(X[tr],y[tr]);
        if cd is not None:sc[te]=X[te]@cd
    ok=np.isfinite(sc); return auc(sc[ok],y[ok])
def sep(sc,R,L):
    return np.nan if R.sum()<5 or L.sum()<5 else float(sc[R].mean()-sc[L].mean())
def bacc(m,hitR,hitL,R,L):
    rr=m&R; ll=m&L
    return np.nan if rr.sum()<5 or ll.sum()<5 else float((hitR[rr].mean()+hitL[ll].mean())/2)

def load(p):
    try:return loadmat(p,squeeze_me=True,struct_as_record=False,variable_names=['obj'])['obj']
    except NotImplementedError:return mat73.loadmat(p)['obj']

m=re.search(r'data_structure_(ANM\d+)_(\d{8})\.mat',KEY or PATH.name)
animal=m.group(1) if m else 'unknown'; session=m.group(2) if m else PATH.stem
try:
    obj=load(PATH); start=vec(fld(obj,'trialStartTimes')); n=len(start); mm=masks(obj,n); hp=fld(obj,'trialPropertiesHash')
    good=hfirst(hp,['GoodTrials'],n)>0; cue=hfirst(hp,['CueTime'],n)
    pstim=hfirst(hp,['PhotostimulationType','PhotostimType','Photostimulation_type'],n)
    perf=(mm['HitR']|mm['HitL']|mm['ErrR']|mm['ErrL'])&(~mm['LickEarly'])&good
    R=mm['HitR']|mm['ErrR']; L=mm['HitL']|mm['ErrL']; ctrl=perf&(pstim==0); ipsi=perf&(pstim==4); contra=perf&(pstim==5)
    U=units(obj)
    allidx=np.where(ctrl|ipsi|contra)[0]; X=matrix(U,allidx,start,cue); pos={int(t):i for i,t in enumerate(allidx)}
    def rows(mask):return np.array([pos[int(t)] for t in np.where(mask)[0] if int(t) in pos],int)
    trmask=ctrl&(mm['HitR']|mm['HitL']); tridx=np.where(trmask)[0]; y=(mm['HitR'][tridx]).astype(int); trrows=np.array([pos[int(t)] for t in tridx],int)
    pred=cvauc(X[trrows],y) if len(U)>5 and (y==0).sum()>=5 and (y==1).sum()>=5 else np.nan
    cd=fitcd(X[trrows],y) if len(tridx) and len(U)>5 else None
    score=X@cd if cd is not None else np.full(len(X),np.nan)
    def local(mask):
        ids=np.where(mask)[0]; rr=np.array([pos[int(t)] for t in ids if int(t) in pos],int)
        rlab=np.array([R[t] for t in ids if int(t) in pos],bool); llab=np.array([L[t] for t in ids if int(t) in pos],bool)
        return sep(score[rr],rlab,llab)
    s0=local(ctrl); si=local(ipsi); sc=local(contra)
    ba0=bacc(ctrl,mm['HitR'],mm['HitL'],R,L); bai=bacc(ipsi,mm['HitR'],mm['HitL'],R,L); bac=bacc(contra,mm['HitR'],mm['HitL'],R,L)
    result={'animal':animal,'session':session,'s3_key':KEY,'n_units':len(U),'predictor_cv_auc':pred,'control_separation':s0,
      'ipsi_recovery':float(si/s0) if np.isfinite(si) and np.isfinite(s0) and abs(s0)>1e-12 else None,
      'contra_recovery':float(sc/s0) if np.isfinite(sc) and np.isfinite(s0) and abs(s0)>1e-12 else None,
      'control_balanced_accuracy':ba0,
      'ipsi_behavior_recovery':float(bai/ba0) if np.isfinite(bai) and np.isfinite(ba0) and ba0>0 else None,
      'contra_behavior_recovery':float(bac/ba0) if np.isfinite(bac) and np.isfinite(ba0) and ba0>0 else None,
      'n_control':int(ctrl.sum()),'n_ipsi4':int(ipsi.sum()),'n_contra5':int(contra.sum()),
      'note':'Frozen raw spike/trial alm-3 analysis: putative pyramidal units; final 400 ms before cue; predictor from correct control trials only.'}
except Exception as e:
    result={'animal':animal,'session':session,'s3_key':KEY,'error':repr(e)}
fn=OUT/f'{animal}_{session}.json'; fn.write_text(json.dumps(result,indent=2,default=str)); print(json.dumps(result,indent=2,default=str))
