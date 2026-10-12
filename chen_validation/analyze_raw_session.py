import os, re, json, shutil, urllib.request
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.stats import rankdata
from remotezip import RemoteZip
import mat73

OUT=Path('analysis/chen_raw_sessions'); OUT.mkdir(parents=True,exist_ok=True)
TMP=Path('chen_validation/tmp_raw.mat')
ARCH=os.environ.get('ARCHIVE_INDEX','direct')
MEMBER=os.environ['MEMBER']

STANDARD={'BAYLORGC4','BAYLORGC12','BAYLORGC13','BAYLORGC15','BAYLORGC17','BAYLORGC18','BAYLORGC19','BAYLORGC21','BAYLORGC22','BAYLORGC25','BAYLORGC26','BAYLORGC50','BAYLORGC72','BAYLORGC75','BAYLORGC86','BAYLORGC87','BAYLORGC93','BAYLORGC95'}

def fld(x,name,default=None):
    if isinstance(x,dict): return x.get(name,default)
    return getattr(x,name,default)

def str1(x):
    if x is None:return ''
    if isinstance(x,bytes): return x.decode(errors='replace')
    if isinstance(x,str): return x
    if isinstance(x,(list,tuple)) and len(x): return str1(x[0])
    a=np.atleast_1d(x).ravel()
    if not len(a): return ''
    y=a[0]
    if y is x:return str(x)
    return str1(y)

def vec(x,n=None,dtype=float):
    if x is None:return np.array([],dtype=dtype)
    a=np.asarray(x).squeeze()
    try:a=a.astype(dtype)
    except:
        vals=[]
        for z in np.atleast_1d(a).ravel():
            try:vals.append(dtype(np.asarray(z).squeeze()))
            except:vals.append(np.nan if dtype is float else 0)
        a=np.array(vals,dtype=dtype)
    a=np.atleast_1d(a).ravel()
    if n is not None and len(a)<n:
        b=np.full(n,np.nan if dtype is float else 0,dtype=dtype); b[:len(a)]=a; a=b
    return a

def as_list(x):
    if x is None:return []
    if isinstance(x,list):return x
    if isinstance(x,tuple):return list(x)
    if isinstance(x,np.ndarray):return list(x.ravel())
    return [x]

def get_hash_value(h,key=None,index=None):
    vals=as_list(fld(h,'value',[]))
    if index is not None: return vals[index]
    keys=[str1(z) for z in as_list(fld(h,'keyNames',[]))]
    for i,k in enumerate(keys):
        if k.lower()==key.lower(): return vals[i]
    raise KeyError(key)

def auc(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)<2 or len(b)<2:return np.nan
    z=np.r_[a,b]; r=rankdata(z,method='average'); n1=len(a); n2=len(b)
    u=np.sum(r[:n1])-n1*(n1+1)/2
    return float(u/(n1*n2))

def extract_file():
    if ARCH=='direct':
        key=MEMBER.split('/')[-1]
        url=f'https://zenodo.org/records/6713616/files/{key}?download=1'
        urllib.request.urlretrieve(url,TMP)
    else:
        url=f'https://zenodo.org/records/6713616/files/datafiles_{ARCH}_of_8.zip?download=1'
        with RemoteZip(url) as rz, rz.open(MEMBER) as src, open(TMP,'wb') as dst:
            shutil.copyfileobj(src,dst,length=8*1024*1024)

def trial_side_masks(obj,ntr):
    names=[str1(z) for z in as_list(fld(obj,'trialTypeStr',[]))]
    mat=np.asarray(fld(obj,'trialTypeMat'))
    if mat.ndim==1: mat=mat.reshape(1,-1)
    # mat73 can preserve MATLAB transposition; orient rows to trial types
    if mat.shape[0]!=len(names) and mat.shape[1]==len(names): mat=mat.T
    def row(name):
        if name not in names:return np.zeros(ntr,bool)
        r=np.asarray(mat[names.index(name),:]).squeeze().astype(float)
        if len(r)<ntr:r=np.pad(r,(0,ntr-len(r)))
        return r[:ntr]>0
    return {k:row(k) for k in ['HitR','HitL','ErrR','ErrL','LickEarly']}

def spike_counts(unit,trial_ids,trial_start,cue,w0=-1.7,w1=-0.9):
    ev=np.asarray(fld(unit,'eventTimes',[]),float).ravel(); et=np.asarray(fld(unit,'eventTrials',[]),int).ravel()
    out=[]
    for tid in trial_ids:
        i=tid-1
        if i<0 or i>=len(trial_start) or i>=len(cue): continue
        lo=trial_start[i]+cue[i]+w0; hi=trial_start[i]+cue[i]+w1
        s=ev[et==tid]
        out.append(np.sum((s>=lo)&(s<hi)))
    return np.asarray(out,float)

def analyze(obj,animal,session):
    trial_ids=np.asarray(fld(obj,'trialIds'),int).ravel(); ntr=int(np.nanmax(trial_ids))
    masks=trial_side_masks(obj,ntr)
    tph=fld(obj,'trialPropertiesHash')
    good=vec(get_hash_value(tph,key='GoodTrials'),ntr,float)[:ntr]>0
    try: cue=vec(get_hash_value(tph,key='CueTime'),ntr,float)[:ntr]
    except:
        try: cue=vec(get_hash_value(tph,key='Cue_Time'),ntr,float)[:ntr]
        except: cue=vec(get_hash_value(tph,index=2),ntr,float)[:ntr]
    try: pstim=vec(get_hash_value(tph,key='PhotostimType'),ntr,float)[:ntr]
    except:
        try: pstim=vec(get_hash_value(tph,key='Photostim_type'),ntr,float)[:ntr]
        except: pstim=vec(get_hash_value(tph,index=4),ntr,float)[:ntr]
    trial_start=vec(fld(obj,'trialStartTimes'),ntr,float)[:ntr]
    performing=(masks['HitR']|masks['ErrR']|masks['HitL']|masks['ErrL']) & (~masks['LickEarly']) & good
    right=(masks['HitR']|masks['ErrR']); left=(masks['HitL']|masks['ErrL']); hit=(masks['HitR']|masks['HitL'])
    control=performing & (pstim==0)
    outcome={}
    for hemi,code in [('left',1),('right',2)]:
        stim=performing & (pstim==code)
        diffs=[]; detail={}
        for lab,sm in [('R',right),('L',left)]:
            ci=np.where(control&sm)[0]; si=np.where(stim&sm)[0]
            if len(ci)>=10 and len(si)>=5:
                ca=float(np.mean(hit[ci])); sa=float(np.mean(hit[si])); diffs.append(ca-sa)
                detail[lab]={'n_control':int(len(ci)),'n_stim':int(len(si)),'control_accuracy':ca,'stim_accuracy':sa,'impairment':ca-sa}
        outcome[hemi]={'impairment':float(np.mean(diffs)) if len(diffs)==2 else None,'detail':detail}
    esh=fld(obj,'eventSeriesHash')
    units=as_list(fld(esh,'value',[]))
    predictors={'left':[],'right':[]}; unit_counts={'left':0,'right':0}
    trial_R_ids=(np.where(control&right)[0]+1).astype(int); trial_L_ids=(np.where(control&left)[0]+1).astype(int)
    for unit in units:
        hemi=str1(fld(unit,'hemisphere','')).lower(); h='left' if 'left' in hemi else ('right' if 'right' in hemi else None)
        if h is None: continue
        cell=str1(fld(unit,'cellType','')).lower()
        try: qual=float(np.asarray(fld(unit,'manual_quality_score',np.nan)).squeeze())
        except: qual=np.nan
        if 'pyramidal' not in cell or 'fs' in cell or not np.isfinite(qual) or qual!=1: continue
        evtr=np.asarray(fld(unit,'eventTrials',[]),int).ravel()
        if not len(evtr):continue
        mn=max(1,int(np.min(evtr))); mx=min(ntr,int(np.max(evtr)))
        rr=trial_R_ids[(trial_R_ids>=mn)&(trial_R_ids<=mx)]; ll=trial_L_ids[(trial_L_ids>=mn)&(trial_L_ids<=mx)]
        if len(rr)<10 or len(ll)<10:continue
        cr=spike_counts(unit,rr,trial_start,cue); cl=spike_counts(unit,ll,trial_start,cue)
        a=auc(cr,cl)
        if np.isfinite(a): predictors[h].append(2*abs(a-0.5)); unit_counts[h]+=1
    pred={h:(float(np.mean(v)) if len(v)>=5 else None) for h,v in predictors.items()}
    return {'animal':animal,'session':session,'member':MEMBER,'predictor':pred,'n_eligible_units':unit_counts,'outcome':outcome,
            'n_control_trials':int(control.sum()),'n_trials':ntr,
            'analysis_note':'Raw spike/trial reanalysis. Predictor uses manually curated putative pyramidal units and control trials only in the source-native early-delay interval; outcome uses same-session unilateral photoinhibition behavioral impairment.'}

def load_obj(path):
    try:
        d=loadmat(path,squeeze_me=True,struct_as_record=False,variable_names=['obj']); return d['obj'],'scipy'
    except NotImplementedError:
        d=mat73.loadmat(path); return d['obj'],'mat73'

m=re.search(r'data_structure_(BAYLORGC\d+)_(\d{4}_\d{2}_\d{2})\.mat',MEMBER)
if not m: raise RuntimeError('cannot parse session from '+MEMBER)
animal,session=m.group(1),m.group(2)
if animal not in STANDARD: raise RuntimeError('nonstandard animal')
try:
    extract_file(); obj,loader=load_obj(TMP); result=analyze(obj,animal,session); result['loader']=loader
except Exception as e:
    result={'animal':animal,'session':session,'member':MEMBER,'error':repr(e)}
finally:
    if TMP.exists(): TMP.unlink()
fn=OUT/f'{animal}_{session}.json'; fn.write_text(json.dumps(result,indent=2,default=str)); print(json.dumps(result,indent=2,default=str),flush=True)
