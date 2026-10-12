import os,json,shutil,re
from pathlib import Path
import numpy as np
import h5py
from remotezip import RemoteZip
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.metrics import accuracy_score

URL='https://zenodo.org/records/17973874/files/2026_Ghaderi_NCOMMS_Zenodo_data_code.zip?download=1'
MEMBER=os.environ['MEMBER'];SEED=20261012;NREP=100;NNEUR=50
REGIONS=['A1','wS1','wS2','wM2','ALM']
EPOCHS={'audio':(-0.05,0.30),'delay':(0.40,0.85),'whisker':(0.95,1.20)}
OUT=Path('analysis/ghaderi_raw_ephys');OUT.mkdir(parents=True,exist_ok=True)
TMP=Path('/tmp/ghaderi_raw_session.nwb')

def read_str(ds):
    try:
        a=ds.asstr()[()]
        return np.asarray(a,dtype=str).ravel()
    except Exception:
        a=ds[()];out=[]
        for z in np.asarray(a).ravel():
            if isinstance(z,bytes):out.append(z.decode(errors='replace'))
            elif isinstance(z,str):out.append(z)
            elif isinstance(z,h5py.Reference):
                try:
                    v=ds.file[z][()]
                    if isinstance(v,bytes):out.append(v.decode(errors='replace'))
                    elif isinstance(v,np.ndarray) and v.dtype.kind in 'SU':out.append(''.join(x.decode() if isinstance(x,bytes) else str(x) for x in v.ravel()))
                    else:out.append(str(v))
                except:out.append('')
            else:out.append(str(z))
        return np.asarray(out,dtype=str)
def mean_abs_dev(x):
    x=np.asarray(x,float);m=np.nanmean(x);return float(np.nanmean(np.abs(x-m)))
def window_means(t,x,starts,lo,hi):
    t=np.asarray(t,float);x=np.asarray(x,float).reshape(-1);starts=np.asarray(starts,float);out=np.full(len(starts),np.nan)
    ok=np.isfinite(x);x0=np.where(ok,x,0.0);n0=np.cumsum(ok.astype(int));s0=np.cumsum(x0)
    for i,s in enumerate(starts):
        a=np.searchsorted(t,s+lo,'left');b=np.searchsorted(t,s+hi,'right')
        if b<=a:continue
        sm=s0[b-1]-(s0[a-1] if a else 0);nn=n0[b-1]-(n0[a-1] if a else 0)
        if nn:out[i]=sm/nn
    return out
def quiet_trials(f,starts):
    wt=f['/processing/behavior/BehavioralTimeSeries/C2Whisker_Angle/timestamps'][()]
    wa=f['/processing/behavior/BehavioralTimeSeries/C2Whisker_Angle/data'][()].astype(float).reshape(-1)
    ws=np.r_[0.0,np.abs(np.diff(wa))]
    jt=f['/processing/behavior/BehavioralTimeSeries/Jaw_Coordinate/timestamps'][()]
    jaw=f['/processing/behavior/BehavioralTimeSeries/Jaw_Coordinate/data'][()].astype(float)
    jx=jaw[:,0];jy=jaw[:,1];jm=np.sqrt((jy-np.nanmedian(jy))**2+(jx-np.nanmedian(jx))**2)
    # The source MATLAB uses mode(position); median is a robust reconstruction sensitivity.
    out=np.ones(len(starts),bool)
    for t,x in [(wt,ws),(jt,jm)]:
        base=window_means(t,x,starts,-1.0,0.0);pre=window_means(t,x,starts,0.8,1.0)
        th=np.nanmedian(base)+mean_abs_dev(base)
        out &= np.isfinite(pre)&(pre<=th)
    return out
def unit_spikes(f):
    flat=f['/units/spike_times'][()].astype(float).ravel();idx=f['/units/spike_times_index'][()].astype(int).ravel();out=[];a=0
    for b in idx:out.append(flat[a:b]);a=b
    return out
def counts(spikes,starts,lo,hi):
    X=np.zeros((len(starts),len(spikes)),float)
    for j,sp in enumerate(spikes):
        a=np.searchsorted(sp,starts+lo,'left');b=np.searchsorted(sp,starts+hi,'left');X[:,j]=b-a
    return X
def cv50(X,y,session_seed):
    y=np.asarray(y,int);n0=(y==0).sum();n1=(y==1).sum()
    if X.shape[1]<NNEUR or min(n0,n1)<5:return None
    vals=[]
    for rep in range(NREP):
        rng=np.random.default_rng(session_seed+rep)
        feat=rng.choice(X.shape[1],NNEUR,replace=False)
        a=np.where(y==0)[0];b=np.where(y==1)[0];n=min(len(a),len(b));aa=rng.choice(a,n,replace=False);bb=rng.choice(b,n,replace=False);ix=np.r_[aa,bb];rng.shuffle(ix)
        xx=X[ix][:,feat];yy=y[ix];k=min(5,int((yy==0).sum()),int((yy==1).sum()))
        if k<2:continue
        cv=StratifiedKFold(n_splits=k,shuffle=True,random_state=session_seed+1000+rep);p=np.full(len(yy),-99)
        for tr,te in cv.split(xx,yy):
            model=make_pipeline(StandardScaler(),LinearSVC(C=1.0,dual='auto',max_iter=10000,random_state=session_seed+rep));model.fit(xx[tr],yy[tr]);p[te]=model.predict(xx[te])
        vals.append(accuracy_score(yy,p)*100.0)
    return float(np.mean(vals)) if vals else None

def main():
    with RemoteZip(URL) as rz,rz.open(MEMBER) as src,open(TMP,'wb') as dst:shutil.copyfileobj(src,dst,8<<20)
    rec={'member':MEMBER,'session':Path(MEMBER).stem,'rows':[],'method_note':'Post hoc raw-NWB reconstruction sensitivity. Source trial classes and 50-neuron population size are retained; fold-wise standardization and epoch-integrated spike counts are used instead of author-generated 10-ms processed decoder arrays.'}
    with h5py.File(TMP,'r') as f:
        starts=f['/intervals/trials/start_time'][()].astype(float).ravel();trial=f['/intervals/trials/trial_type'][()].astype(float).ravel();lick=f['/intervals/trials/lick_flag'][()].astype(bool).ravel();early=f['/intervals/trials/early_lick'][()].astype(bool).ravel();quiet=quiet_trials(f,starts)
        loc=read_str(f['/units/location']);ccf=read_str(f['/units/ccf_location']);sp=unit_spikes(f);n=min(len(sp),len(loc));sp=sp[:n];loc=loc[:n];ccf=ccf[:n]
        elig=(~early)&quiet;class1=elig&(trial==1)&lick;class0=elig&(trial==3)&(~lick);use=class1|class0;y=class1[use].astype(int);st=starts[use]
        rec['n_trials_class1']=int(class1.sum());rec['n_trials_class0']=int(class0.sum());rec['unit_location_unique']=sorted(set(loc.tolist()))
        for reg in REGIONS:
            ui=np.where(np.char.strip(loc.astype(str))==reg)[0]
            if len(ui)<NNEUR:continue
            rs=[sp[i] for i in ui];seed=SEED+sum(ord(c) for c in rec['session']+reg)
            for ep,(lo,hi) in EPOCHS.items():
                X=counts(rs,st,lo,hi);acc=cv50(X,y,seed+int((lo+2)*1000))
                if acc is not None:rec['rows'].append({'region':reg,'epoch':ep,'decoder_accuracy':acc,'n_units':int(len(ui)),'n_class1':int((y==1).sum()),'n_class0':int((y==0).sum())})
    try:TMP.unlink()
    except:pass
    fn=re.sub(r'[^A-Za-z0-9_.-]','_',rec['session'])+'.json';(OUT/fn).write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2),flush=True)
if __name__=='__main__':main()
