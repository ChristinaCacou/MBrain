import os,re,json,shutil,urllib.request
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.stats import rankdata,spearmanr
from remotezip import RemoteZip
import mat73

SEED=20261012
OUT=Path('analysis/chen_robust_sessions'); OUT.mkdir(parents=True,exist_ok=True)
TMP=Path('raw_replication/tmp_chen_robust.mat')
ARCH=os.environ.get('ARCHIVE_INDEX','direct'); MEMBER=os.environ['MEMBER']

def fld(x,n,d=None):return x.get(n,d) if isinstance(x,dict) else getattr(x,n,d)
def s1(x):
    if x is None:return ''
    if isinstance(x,bytes):return x.decode(errors='replace')
    if isinstance(x,str):return x
    a=np.atleast_1d(x).ravel();return '' if len(a)==0 else s1(a[0]) if a[0] is not x else str(x)
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
        b=np.full(n,np.nan if dtype is float else 0,dtype=dtype);b[:len(a)]=a;a=b
    return a
def hval(h,key):
    if isinstance(h,list):
        for z in h:
            if s1(fld(z,'keyNames','')).lower()==key.lower():return fld(z,'value')
        raise KeyError(key)
    ks=[s1(z) for z in alist(fld(h,'keyNames',[]))];vs=alist(fld(h,'value',[]))
    for i,k in enumerate(ks):
        if k.lower()==key.lower():return vs[i]
    raise KeyError(key)
def hfirst(h,names,n):
    for k in names:
        try:return vec(hval(h,k),n,float)
        except Exception:pass
    raise KeyError(names)
def masks(obj,n):
    names=[s1(z) for z in alist(fld(obj,'trialTypeStr',[]))];M=np.asarray(fld(obj,'trialTypeMat'))
    if M.ndim==1:M=M.reshape(1,-1)
    if M.shape[0]!=len(names) and M.ndim==2 and M.shape[1]==len(names):M=M.T
    def r(k):
        for i,nm in enumerate(names):
            if nm.lower()==k.lower():
                x=np.asarray(M[i]).squeeze().astype(float);return np.pad(x,(0,max(0,n-len(x))))[:n]>0
        return np.zeros(n,bool)
    return {k:r(k) for k in ['HitR','HitL','ErrR','ErrL','LickEarly']}
def extract():
    if ARCH=='direct':urllib.request.urlretrieve(f'https://zenodo.org/records/6713616/files/{MEMBER.split("/")[-1]}?download=1',TMP)
    else:
        u=f'https://zenodo.org/records/6713616/files/datafiles_{ARCH}_of_8.zip?download=1'
        with RemoteZip(u) as rz,rz.open(MEMBER) as src,open(TMP,'wb') as dst:shutil.copyfileobj(src,dst,length=8*1024*1024)
def loadobj():
    try:return loadmat(TMP,squeeze_me=True,struct_as_record=False,variable_names=['obj'])['obj']
    except NotImplementedError:return mat73.loadmat(TMP)['obj']
def getunits(obj,hemi):
    H=fld(obj,'eventSeriesHash');U=[fld(z,'value') for z in H if fld(z,'value') is not None] if isinstance(H,list) else alist(fld(H,'value',[]));out=[]
    for u in U:
        hs=s1(fld(u,'hemisphere','')).lower();c=s1(fld(u,'cellType','')).lower()
        try:q=float(np.asarray(fld(u,'manual_quality_score',np.nan)).squeeze())
        except:q=np.nan
        if hemi not in hs or 'pyramidal' not in c or 'fs' in c or q!=1:continue
        if len(np.atleast_1d(fld(u,'eventTrials',[]))):out.append(u)
    return out
def matrix(U,idx,start,cue):
    idx=np.asarray(idx,int);X=np.zeros((len(idx),len(U)),float)
    for j,u in enumerate(U):
        ev=np.asarray(fld(u,'eventTimes',[]),float).ravel();et=np.asarray(fld(u,'eventTrials',[]),int).ravel()
        for i,t0 in enumerate(idx):
            sp=ev[et==(t0+1)];lo=start[t0]+cue[t0]-0.4;hi=start[t0]+cue[t0];X[i,j]=np.sum((sp>=lo)&(sp<hi))
    return X
def fitcd(X,y):
    if (y==1).sum()<1 or (y==0).sum()<1:return None
    d=X[y==1].mean(0)-X[y==0].mean(0);z=np.linalg.norm(d);return None if not np.isfinite(z) or z==0 else d/z
def auc(sc,y):
    a=sc[y==1];b=sc[y==0]
    if len(a)<1 or len(b)<1:return np.nan
    r=rankdata(np.r_[a,b]);n1=len(a);n0=len(b);q=(r[:n1].sum()-n1*(n1+1)/2)/(n1*n0);return float(max(q,1-q))
def cvauc(X,y):
    k=int(min(5,(y==1).sum(),(y==0).sum()))
    if k<2:return np.nan
    rng=np.random.default_rng(SEED);spl=[]
    for lab in [0,1]:q=np.where(y==lab)[0].copy();rng.shuffle(q);spl.append(np.array_split(q,k))
    sc=np.full(len(y),np.nan)
    for f in range(k):
        te=np.r_[spl[0][f],spl[1][f]];tr=np.setdiff1d(np.arange(len(y)),te);cd=fitcd(X[tr],y[tr])
        if cd is not None:sc[te]=X[te]@cd
    ok=np.isfinite(sc);return auc(sc[ok],y[ok])
def bacc(mask,mm,R,L):
    rr=mask&R;ll=mask&L
    return np.nan if rr.sum()<5 or ll.sum()<5 else float((mm['HitR'][rr].mean()+mm['HitL'][ll].mean())/2)
def sep(sc,mask,R,L):
    rr=mask&R;ll=mask&L
    return np.nan if rr.sum()<5 or ll.sum()<5 else float(sc[rr].mean()-sc[ll].mean())

def analyze(obj,animal,session):
    start=vec(fld(obj,'trialStartTimes'));n=len(start);mm=masks(obj,n);hp=fld(obj,'trialPropertiesHash');good=hfirst(hp,['GoodTrials'],n)>0;cue=hfirst(hp,['CueTime','Cue_Time'],n);pstim=hfirst(hp,['PhotostimType','Photostim_type','PhotostimulationType'],n)
    perf=(mm['HitR']|mm['HitL']|mm['ErrR']|mm['ErrL'])&(~mm['LickEarly'])&good;R=mm['HitR']|mm['ErrR'];L=mm['HitL']|mm['ErrL'];ctrl=perf&(pstim==0);leftstim=perf&(pstim==1);rightstim=perf&(pstim==2)
    UL=getunits(obj,'left');UR=getunits(obj,'right');allidx=np.where(perf)[0];XL=matrix(UL,allidx,start,cue);XR=matrix(UR,allidx,start,cue);pos={int(t):i for i,t in enumerate(allidx)}
    ctrl_correct=ctrl&(mm['HitR']|mm['HitL']);cc=np.where(ctrl_correct)[0];rows=np.array([pos[int(t)] for t in cc],int);y=mm['HitR'][cc].astype(int)
    cdL=fitcd(XL[rows],y) if len(UL)>5 else None;cdR=fitcd(XR[rows],y) if len(UR)>5 else None
    local_auc_L=cvauc(XL[rows],y) if cdL is not None else np.nan;local_auc_R=cvauc(XR[rows],y) if cdR is not None else np.nan
    rhos=[];within=[];fixed={10:[],20:[]};rng=np.random.default_rng(SEED)
    rtr=np.where(ctrl&mm['HitR'])[0];ltr=np.where(ctrl&mm['HitL'])[0]
    if len(UL)>5 and len(UR)>5 and len(rtr)>=10 and len(ltr)>=10:
        for rep in range(30):
            rr=rtr.copy();ll=ltr.copy();rng.shuffle(rr);rng.shuffle(ll);nr=max(1,len(rr)//2);nl=max(1,len(ll)//2);train=np.r_[rr[:nr],ll[:nl]];test=np.r_[rr[nr:],ll[nl:]]
            trrows=np.array([pos[int(t)] for t in train],int);te=np.array([pos[int(t)] for t in test],int);ty=np.array([R[t] for t in train],int)
            a=fitcd(XL[trrows],ty);b=fitcd(XR[trrows],ty)
            if a is not None and b is not None and len(te)>=6:
                pl=XL[te]@a;pr=XR[te]@b;rho=spearmanr(pl,pr).statistic
                if np.isfinite(rho):rhos.append(float(rho))
                wr=[];labs=np.array([R[t] for t in test],int)
                for lab in [0,1]:
                    z=labs==lab
                    if z.sum()>=3:
                        q=spearmanr(pl[z],pr[z]).statistic
                        if np.isfinite(q):wr.append(float(q))
                if wr:within.append(float(np.mean(wr)))
            for nf in [10,20]:
                if XL.shape[1] < nf or XR.shape[1] < nf or len(te)<6:continue
                ur=np.random.default_rng(SEED + 10000*nf + rep)
                il=ur.choice(XL.shape[1],nf,replace=False);ir=ur.choice(XR.shape[1],nf,replace=False)
                af=fitcd(XL[trrows][:,il],ty);bf=fitcd(XR[trrows][:,ir],ty)
                if af is None or bf is None:continue
                q=spearmanr(XL[te][:,il]@af,XR[te][:,ir]@bf).statistic
                if np.isfinite(q):fixed[nf].append(float(q))
    coupling=float(np.mean(rhos)) if rhos else np.nan;coupling_within=float(np.mean(within)) if within else np.nan
    result={'animal':animal,'session':session,'member':MEMBER,'n_units_left':len(UL),'n_units_right':len(UR),'n_control':int(ctrl.sum()),'n_leftstim':int(leftstim.sum()),'n_rightstim':int(rightstim.sum()),
            'control_bilateral_coupling':coupling,'control_bilateral_within_condition_coupling':coupling_within,
            'control_bilateral_coupling_fixed10':float(np.mean(fixed[10])) if fixed[10] else None,
            'control_bilateral_coupling_fixed20':float(np.mean(fixed[20])) if fixed[20] else None,
            'local_cv_auc_left':local_auc_L,'local_cv_auc_right':local_auc_R,
            'transferable_local_cv_auc':float(np.nanmean([local_auc_L,local_auc_R])) if np.isfinite(local_auc_L) or np.isfinite(local_auc_R) else None}
    ba0=bacc(ctrl,mm,R,L);result['control_behavior_accuracy']=float(ba0) if np.isfinite(ba0) else None
    for hemi,cd,X,stim in [('left',cdL,XL,leftstim),('right',cdR,XR,rightstim)]:
        if cd is None:
            result[f'{hemi}_neuronal_recovery']=None;result[f'{hemi}_control_separation']=None;result[f'{hemi}_stim_separation']=None;result[f'{hemi}_neuronal_delta']=None
        else:
            sc=np.full(n,np.nan);sc[allidx]=X@cd;s0=sep(sc,ctrl,R,L);ss=sep(sc,stim,R,L)
            result[f'{hemi}_control_separation']=float(s0) if np.isfinite(s0) else None;result[f'{hemi}_stim_separation']=float(ss) if np.isfinite(ss) else None
            result[f'{hemi}_neuronal_recovery']=float(ss/s0) if np.isfinite(ss) and np.isfinite(s0) and abs(s0)>1e-12 else None;result[f'{hemi}_neuronal_delta']=float(ss-s0) if np.isfinite(ss) and np.isfinite(s0) else None
        bas=bacc(stim,mm,R,L);result[f'{hemi}_stim_behavior_accuracy']=float(bas) if np.isfinite(bas) else None
        result[f'{hemi}_behavior_recovery']=float(bas/ba0) if np.isfinite(bas) and np.isfinite(ba0) and ba0>0 else None;result[f'{hemi}_behavior_delta']=float(bas-ba0) if np.isfinite(bas) and np.isfinite(ba0) else None
    result['note']='Raw Chen spike/trial analysis. Fixed-10 and fixed-20 coupling and additive behavioral/neuronal endpoints are post hoc robustness analyses; no published figure coordinates used.'
    return result

m=re.search(r'data_structure_(BAYLORGC\d+)_(\d{4}_\d{2}_\d{2})\.mat',MEMBER);animal,session=m.group(1),m.group(2)
try:extract();result=analyze(loadobj(),animal,session)
except Exception as e:result={'animal':animal,'session':session,'member':MEMBER,'error':repr(e)}
finally:
    if TMP.exists():TMP.unlink()
(OUT/f'{animal}_{session}.json').write_text(json.dumps(result,indent=2,default=str));print(json.dumps(result,indent=2,default=str),flush=True)
