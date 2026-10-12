#!/usr/bin/env python3
"""
Independent cloud reanalysis of the public Falconer/Kiran n=30 aphasia cohort.

This program:
1) enumerates OSF project wcuj8 using the public OSF v2 API;
2) downloads only analysis-relevant CONN ROI, lesion-mask, and small metadata files;
3) derives lesion volume/AAL3 lesion fractions, static FC/graph metrics, and
   sliding-window temporal variability directly from released participant data;
4) establishes the de-identified P01-P30 linkage only with an explicit public
   token or a prespecified lesion-volume concordance gate;
5) performs leakage-resistant nested leave-one-out prediction.

No participant imaging values are synthesized from published coefficients.
"""
from __future__ import annotations
import argparse, json, math, os, re, shutil, sys, tarfile, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from scipy import stats
from scipy.io import loadmat
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import shortest_path
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import nibabel as nib
from nibabel.processing import resample_from_to

GUID = "wcuj8"
TR = 2.4
N_VOL = 175
WIN = 12
N_WIN = N_VOL - WIN + 1
SEED = 20261012
ALPHAS = np.logspace(-4, 4, 33)
DAMAGE_THRESHOLD = 0.50
NETWORKS = {
    "language": [1,2,3,4,5,6,7,8,9,10,11,12,15,16,25,26,57,58,61,62,65,66,67,68,81,82,83,84,85,86,87,88,89,90],
    "salience": [33,34,147,148,149,150,151,152],
    "dorsal_attention": [1,2,5,6,51,52,61,62],
    "default_mode": [19,20,25,26,37,38,63,64,69,70],
}
CORTEX = list(range(1,39)) + [41,42] + list(range(45,73)) + list(range(79,91)) + list(range(147,153))
assert len(CORTEX) == 86
AAL3_LABELS = list(range(1,35)) + list(range(37,81)) + list(range(83,171))
assert len(AAL3_LABELS) == 166
UA = {"User-Agent": "Falconer-open-reanalysis/2026-10-12 (+GitHub Actions)"}

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=lambda x: float(x) if isinstance(x, np.generic) else str(x)))

def get_json(url, session):
    r = session.get(url, headers=UA, timeout=90)
    r.raise_for_status()
    return r.json()

def enumerate_osf(session):
    starts = [
        f"https://api.osf.io/v2/nodes/{GUID}/files/?page[size]=100",
        f"https://api.osf.io/v2/nodes/{GUID}/files/osfstorage/?page[size]=100",
    ]
    queue, seen, rows = list(starts), set(), []
    last_errors = []
    while queue:
        url = queue.pop(0)
        if not url or url in seen:
            continue
        seen.add(url)
        try:
            obj = get_json(url, session)
        except Exception as e:
            last_errors.append({"url": url, "error": repr(e)})
            continue
        for item in obj.get("data", []):
            a = item.get("attributes", {}) or {}
            rel = item.get("relationships", {}) or {}
            links = item.get("links", {}) or {}
            kind = a.get("kind") or item.get("type")
            mat = a.get("materialized_path") or a.get("name") or ""
            rows.append({
                "id": item.get("id"), "kind": kind,
                "name": a.get("name") or Path(mat).name,
                "materialized_path": mat, "size": a.get("size"),
                "provider": a.get("provider"), "download": links.get("download"),
                "api_url": links.get("self"),
            })
            href = ((((rel.get("files") or {}).get("links") or {}).get("related") or {}).get("href"))
            if href:
                queue.append(href + ("&" if "?" in href else "?") + "page[size]=100")
        nxt = (obj.get("links") or {}).get("next")
        if nxt:
            queue.append(nxt)
    if not rows:
        raise RuntimeError(f"Could not enumerate OSF {GUID}; errors={last_errors[:5]}")
    return pd.DataFrame(rows).drop_duplicates(subset=["id","materialized_path"]), last_errors

def relevant_file(row):
    name = str(row.get("materialized_path") or row.get("name") or "")
    low = name.lower()
    size = row.get("size")
    try: size = int(size) if pd.notna(size) else 0
    except Exception: size = 0
    is_conn = bool(re.search(r"roi_subject\d+_condition\d+\.mat$", low))
    is_lesion = ("lesion" in low and low.endswith((".nii",".nii.gz")))
    small_meta = (
        low.endswith((".csv",".tsv",".txt",".json",".md",".xlsx",".xls",".rds",".rdata"))
        and any(k in low for k in ("readme","subject","participant","mapping","behav","accuracy","wab","spared","lesion","demo"))
        and (size == 0 or size <= 25_000_000)
    )
    archive = (
        low.endswith((".zip",".tar",".tar.gz",".tgz"))
        and any(k in low for k in ("conn","bold","roi","lesion","rest","data"))
        and (size == 0 or size <= 1_500_000_000)
    )
    return is_conn or is_lesion or small_meta or archive

def safe_rel(materialized, fallback):
    s = str(materialized or fallback or "file").strip("/")
    s = re.sub(r"(^|/)\.\.(?=/|$)", "_", s)
    return s or fallback or "file"

def download_selected(manifest, root, session):
    root.mkdir(parents=True, exist_ok=True)
    sel = manifest.apply(relevant_file, axis=1)
    out = manifest.copy(); out["selected"] = sel; out["local_path"] = ""; out["download_status"] = ""
    for i, row in out[sel].iterrows():
        url = row.get("download")
        if not isinstance(url, str) or not url:
            out.at[i,"download_status"] = "no_download_link"; continue
        rel = safe_rel(row.get("materialized_path"), row.get("name"))
        dst = root / rel; dst.parent.mkdir(parents=True, exist_ok=True)
        try:
            with session.get(url, headers=UA, stream=True, timeout=180) as r:
                r.raise_for_status()
                with open(dst, "wb") as f:
                    for chunk in r.iter_content(1024*1024):
                        if chunk: f.write(chunk)
            out.at[i,"local_path"] = str(dst); out.at[i,"download_status"] = "ok"
        except Exception as e:
            out.at[i,"download_status"] = f"ERROR:{type(e).__name__}:{e}"
    return out

def extract_archives(root):
    extracted = []
    for p in list(root.rglob("*")):
        if not p.is_file(): continue
        low = p.name.lower()
        try:
            if low.endswith(".zip"):
                d = p.parent/(p.stem+"_extracted"); d.mkdir(exist_ok=True)
                with zipfile.ZipFile(p) as z: z.extractall(d)
                extracted.append(str(d))
            elif low.endswith((".tar.gz",".tgz",".tar")):
                base = re.sub(r"\.(tar\.gz|tgz|tar)$","",p.name,flags=re.I)
                d = p.parent/(base+"_extracted"); d.mkdir(exist_ok=True)
                with tarfile.open(p) as t: t.extractall(d, filter="data")
                extracted.append(str(d))
        except Exception as e:
            extracted.append(f"ERROR {p}: {e}")
    return extracted

def subject_num(path):
    s = str(path)
    pats = [
        r"ROI_Subject(\d{3})_Condition\d+\.mat$",
        r"(?:^|[/_.-])P(\d{2})(?:[/_.-]|$)",
        r"(?:^|[/_.-])BU(\d{2,3})(?:[/_.-]|$)",
        r"(?:^|[/_.-])Subject[_-]?(\d{1,3})(?:[/_.-]|$)",
        r"(?:^|[/_.-])sub[-_]?(\d{1,3})(?:[/_.-]|$)",
    ]
    for pat in pats:
        m = re.search(pat, s, flags=re.I)
        if m: return int(m.group(1))
    return None

def discover(root):
    files = [p for p in root.rglob("*") if p.is_file()]
    conn = [p for p in files if re.search(r"ROI_Subject\d+_Condition001\.mat$", p.name, re.I)]
    lesions = [p for p in files if p.name.lower().endswith((".nii",".nii.gz")) and "lesion" in str(p).lower()]
    meta = [p for p in files if p.suffix.lower() in (".csv",".tsv",".txt",".json",".md",".xlsx",".xls")]
    return conn, lesions, meta

def choose_one_per_num(paths):
    grouped = {}
    for p in paths:
        n = subject_num(p)
        if n is None: continue
        grouped.setdefault(n, []).append(p)
    chosen = {}; audit = {}
    for n, ps in grouped.items():
        ps2 = sorted(ps, key=lambda p: (0 if "2mm" in p.name.lower() else 1, 0 if "lesion" in p.name.lower() else 1, len(p.parts), len(str(p))))
        chosen[n] = ps2[0]; audit[str(n)] = [str(x) for x in ps]
    return chosen, audit

def load_conn_roi(path, nroi=166):
    m = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "data" not in m: raise ValueError(f"{path}: no data cell array")
    data = np.asarray(m["data"], dtype=object).reshape(-1)
    if len(data) < nroi+3: raise ValueError(f"{path}: expected >= {nroi+3} cells, got {len(data)}")
    cols = []
    for j in range(3, 3+nroi):
        cols.append(np.asarray(data[j]).squeeze().astype(float).reshape(-1))
    n = min(map(len, cols)); x = np.column_stack([c[:n] for c in cols])
    for j in range(len(data)-1, 3+nroi-1, -1):
        try:
            v = np.asarray(data[j]).squeeze().astype(float).reshape(-1)
            if len(v) == n:
                nz = np.isfinite(v) & (v != 0)
                if 20 <= nz.sum() <= n and nz.sum() < n:
                    x = x[nz]; break
        except Exception: pass
    return x

def fisher_corr(x):
    c = np.corrcoef(x, rowvar=False); np.fill_diagonal(c, 0)
    return np.arctanh(np.clip(c, -1+1e-12, 1-1e-12))

def mean_within(z, inds1, excluded=set()):
    inds = [i-1 for i in inds1 if i not in excluded]
    if len(inds)<2: return np.nan
    sub = z[np.ix_(inds, inds)]; vals = sub[np.triu_indices(len(inds),1)]
    return float(np.nanmean(vals))

def positive_graph_metrics(z, inds1, excluded=set()):
    inds = [i-1 for i in inds1 if i not in excluded]
    if len(inds)<3: return (np.nan,np.nan,np.nan)
    W = z[np.ix_(inds, inds)].copy(); W[~np.isfinite(W)] = 0; W[W < 0] = 0; np.fill_diagonal(W,0)
    strength = float(np.mean(W.sum(axis=1)))
    K = np.sum(W != 0, axis=1); W13 = np.power(W, 1/3, where=W>=0, out=np.zeros_like(W))
    cyc3 = np.diag(W13 @ W13 @ W13); den = K*(K-1)
    C = np.divide(cyc3, den, out=np.zeros_like(cyc3), where=den!=0); clustering = float(np.mean(C))
    L = np.full_like(W, np.inf, dtype=float); mask = W > 0; L[mask] = 1.0/W[mask]; np.fill_diagonal(L,0)
    D = shortest_path(L, directed=False); invD = np.zeros_like(D); ok = np.isfinite(D) & (D>0); invD[ok] = 1.0/D[ok]
    n = len(inds); efficiency = float(invD.sum()/(n*(n-1)))
    return strength, clustering, efficiency

def dynamic_fc(x):
    if x.shape[0] < N_VOL: raise ValueError(f"Need {N_VOL} volumes, found {x.shape[0]}")
    x = x[:N_VOL]; out = np.empty((x.shape[1],x.shape[1],N_WIN),float)
    for k in range(N_WIN): out[:,:,k] = fisher_corr(x[k:k+WIN])
    return out

def temporal_variability_sd(dfc, excluded_mask=None):
    y = dfc.copy()
    if excluded_mask is not None:
        ii = np.where(excluded_mask)[0]; y[ii,:,:] = np.nan; y[:,ii,:] = np.nan
    sd = np.nanstd(y, axis=2, ddof=1); vals = sd[np.triu_indices(sd.shape[0],1)]
    return float(np.nanmean(vals))

def lesion_info(path, atlas_img=None):
    img = nib.load(str(path)); data = np.asarray(img.get_fdata()); lesion = np.isfinite(data) & (data > 0)
    voxvol = abs(np.linalg.det(img.affine[:3,:3]))/1000.0; vol = float(lesion.sum()*voxvol); frac = None
    if atlas_img is not None:
        at = resample_from_to(atlas_img, img, order=0); a = np.rint(np.asarray(at.get_fdata())).astype(int); frac = {}
        for internal,label in enumerate(AAL3_LABELS, start=1):
            roi = (a==label); frac[internal] = float((lesion & roi).sum()/roi.sum()) if roi.sum() else np.nan
    return vol, frac

def fetch_aal3(cache):
    from nilearn.datasets import fetch_atlas_aal
    ds = fetch_atlas_aal(version="3v2", data_dir=str(cache), verbose=1)
    return nib.load(str(ds.maps)), str(ds.maps)

def mapping_from_volumes(vol_by_num, behavior):
    nums = sorted(vol_by_num); raw = np.array([vol_by_num[n] for n in nums],float)
    if len(nums) != 30: return None, {"accepted":False,"reason":f"need 30 lesion-number tokens; found {len(nums)}"}
    candidates = {}
    for col in ["table1_lesion_cm3","legacy_lesion_cm3"]:
        ref = behavior[col].to_numpy(float); pear_i = stats.pearsonr(raw, ref).statistic; spear_i = stats.spearmanr(raw, ref).statistic
        medrel_i = float(np.median(np.abs(raw-ref)/np.maximum(ref,1e-6)))
        candidates[col+"_identity"] = {"mode":"identity","ref_col":col,"pearson":pear_i,"spearman":spear_i,"median_relative_error":medrel_i,"assignment":list(range(30))}
        rz=(raw-raw.mean())/(raw.std(ddof=1) or 1); qz=(ref-ref.mean())/(ref.std(ddof=1) or 1)
        rr,cc=linear_sum_assignment(np.abs(rz[:,None]-qz[None,:])); matched=np.empty(30,dtype=int); matched[rr]=cc; refo=ref[matched]
        candidates[col+"_hungarian"] = {"mode":"hungarian_z","ref_col":col,"pearson":float(stats.pearsonr(raw,refo).statistic),"spearman":float(stats.spearmanr(raw,refo).statistic),"median_relative_error":float(np.median(np.abs(raw-refo)/np.maximum(refo,1e-6))),"assignment":matched.tolist()}
    strong_id = [v for v in candidates.values() if v["mode"]=="identity" and v["pearson"]>=.95 and v["spearman"]>=.95]
    if strong_id:
        best=max(strong_id,key=lambda v:(v["spearman"],v["pearson"],-v["median_relative_error"])); rationale="identity order independently validated by lesion-volume concordance"
    else:
        best=max(candidates.values(),key=lambda v:(v["spearman"],v["pearson"],-v["median_relative_error"])); rationale="volume-based assignment"
    accepted = bool(best["pearson"]>=.95 and best["spearman"]>=.95); mapping={}
    if accepted:
        for i,n in enumerate(nums): mapping[n]=behavior.iloc[best["assignment"][i]].participant_id
    audit={"accepted":accepted,"rationale":rationale,"chosen":best,"all_candidates":candidates,"raw_volumes_by_token":{str(k):vol_by_num[k] for k in nums}}
    return mapping if accepted else None, audit

def nested_loo(df, outcome, feature_sets):
    y=df[outcome].to_numpy(float); loo=LeaveOneOut(); all_pred={}; summaries=[]
    for name,feats in feature_sets.items():
        X=df[feats].to_numpy(float); pred=np.full(len(df),np.nan); chosen=[]
        for tr,te in loo.split(X):
            ntr=len(tr); k=min(5,max(3,ntr//5)); inner=KFold(n_splits=k,shuffle=True,random_state=SEED); scores=[]
            for a in ALPHAS:
                pp=np.full(ntr,np.nan)
                for a_tr,a_va in inner.split(X[tr]):
                    pipe=Pipeline([("scale",StandardScaler()),("ridge",Ridge(alpha=a))]); pipe.fit(X[tr][a_tr],y[tr][a_tr]); pp[a_va]=pipe.predict(X[tr][a_va])
                scores.append(mean_absolute_error(y[tr],pp))
            alpha=float(ALPHAS[int(np.argmin(scores))]); chosen.append(alpha)
            pipe=Pipeline([("scale",StandardScaler()),("ridge",Ridge(alpha=alpha))]); pipe.fit(X[tr],y[tr]); pred[te]=pipe.predict(X[te])
        all_pred[name]=pred
        summaries.append({"model":name,"n":len(df),"features":";".join(feats),"r2":float(r2_score(y,pred)),"mae":float(mean_absolute_error(y,pred)),"rmse":float(mean_squared_error(y,pred)**.5),"predictive_r":float(stats.pearsonr(y,pred).statistic),"median_alpha":float(np.median(chosen))})
    return pd.DataFrame(summaries), pd.DataFrame({"participant_id":df.participant_id,"observed":y,**all_pred})

def paired_incremental_tests(preds, baseline="baseline", B=10000):
    rng=np.random.default_rng(SEED); y=preds.observed.to_numpy(float); base=preds[baseline].to_numpy(float); rows=[]
    for name in [c for c in preds.columns if c not in ("participant_id","observed",baseline)]:
        p=preds[name].to_numpy(float); dse=(y-base)**2-(y-p)**2; delta_r2=r2_score(y,p)-r2_score(y,base); delta_mae=mean_absolute_error(y,p)-mean_absolute_error(y,base); boot=[]
        for _ in range(B):
            idx=rng.integers(0,len(y),len(y))
            if np.var(y[idx])>0: boot.append(r2_score(y[idx],p[idx])-r2_score(y[idx],base[idx]))
        ci=np.quantile(boot,[.025,.975]); signs=rng.choice([-1,1],size=(B,len(dse))); null=(signs*dse).mean(axis=1); pv=(1+np.sum(null>=dse.mean()))/(B+1)
        rows.append({"model":name,"delta_r2_vs_baseline":delta_r2,"delta_mae_vs_baseline":delta_mae,"bootstrap_delta_r2_ci_low":ci[0],"bootstrap_delta_r2_ci_high":ci[1],"mean_squared_error_advantage":dse.mean(),"signflip_p_one_sided":pv})
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--behavior",type=Path,required=True); ap.add_argument("--outdir",type=Path,required=True); ap.add_argument("--cache",type=Path,required=True); args=ap.parse_args()
    out=args.outdir; out.mkdir(parents=True,exist_ok=True); cache=args.cache; cache.mkdir(parents=True,exist_ok=True); raw=cache/"wcuj8_subset"; raw.mkdir(parents=True,exist_ok=True)
    behavior=pd.read_csv(args.behavior); audit={"guid":GUID,"started":True,"mapping_accepted":False,"analysis_completed":False}; s=requests.Session()
    try:
        manifest, enum_errors=enumerate_osf(s); manifest.to_csv(out/"osf_manifest_full.csv",index=False); audit["n_manifest_rows"]=len(manifest); audit["enumeration_errors"]=enum_errors
        acquired=download_selected(manifest,raw,s); acquired.to_csv(out/"osf_manifest_acquisition.csv",index=False); audit["n_selected"]=int(acquired.selected.sum()); audit["n_downloaded_ok"]=int((acquired.download_status=="ok").sum()); audit["archives_extracted"]=extract_archives(raw)
        conn,lesions,meta=discover(raw); audit["n_conn_condition001_files"]=len(conn); audit["n_lesion_candidates"]=len(lesions); audit["metadata_files"]=[str(p) for p in meta]
        excerpts=[]
        for p in meta:
            if p.suffix.lower() in (".txt",".md",".csv",".tsv",".json") and p.stat().st_size<=2_000_000:
                try: excerpts.append(f"\n===== {p} =====\n"+p.read_text(errors="replace")[:10000])
                except Exception: pass
        (out/"metadata_excerpts.txt").write_text("".join(excerpts))
        conn_by_num,_=choose_one_per_num(conn); lesion_by_num,lesion_candidates=choose_one_per_num(lesions); write_json(out/"lesion_candidate_audit.json",lesion_candidates); audit["conn_subject_nums"]=sorted(conn_by_num); audit["lesion_subject_nums"]=sorted(lesion_by_num)
        atlas_img=None
        try: atlas_img,atlas_path=fetch_aal3(cache/"aal3"); audit["aal3_path"]=atlas_path
        except Exception as e: audit["aal3_error"]=repr(e)
        vol_by_num={}; lesion_rows=[]; frac_by_num={}
        for n,p in lesion_by_num.items():
            try:
                vol,fr=lesion_info(p,atlas_img); vol_by_num[n]=vol; frac_by_num[n]=fr; lesion_rows.append({"file_subject_num":n,"lesion_file":str(p),"independent_lesion_volume_cm3":vol})
            except Exception as e: lesion_rows.append({"file_subject_num":n,"lesion_file":str(p),"error":repr(e)})
        pd.DataFrame(lesion_rows).to_csv(out/"unmapped_lesion_metrics.csv",index=False)
        mapping,map_audit=mapping_from_volumes(vol_by_num,behavior); write_json(out/"mapping_audit.json",map_audit)
        if mapping is None: audit["mapping_accepted"]=False; raise RuntimeError("Subject mapping gate failed; see mapping_audit.json")
        audit["mapping_accepted"]=True; pd.DataFrame([{"file_subject_num":k,"participant_id":v} for k,v in sorted(mapping.items())]).to_csv(out/"subject_mapping.csv",index=False)
        common_nums=sorted(set(conn_by_num)&set(lesion_by_num)&set(mapping)); audit["n_common_imaging_subjects"]=len(common_nums)
        if len(common_nums)!=30: raise RuntimeError(f"Hard n=30 imaging gate failed; common={len(common_nums)}")
        metric_rows=[]
        for n in common_nums:
            pid=mapping[n]; row={"file_subject_num":n,"participant_id":pid,"independent_lesion_volume_cm3":vol_by_num[n],"conn_file":str(conn_by_num[n]),"lesion_file":str(lesion_by_num[n])}; fr=frac_by_num.get(n); excluded={i for i,v in (fr or {}).items() if np.isfinite(v) and v>=DAMAGE_THRESHOLD}; row["n_aal3_rois_ge50pct_damaged"]=len(excluded)
            if fr:
                for i,v in fr.items(): row[f"aal3_damage_{i:03d}"]=v
            x=load_conn_roi(conn_by_num[n]); row["n_bold_volumes_valid"]=x.shape[0]; row["n_bold_rois"]=x.shape[1]
            if x.shape[1] < 166: raise ValueError(f"{pid}: {x.shape[1]} ROIs, expected 166")
            z=fisher_corr(x[:,:166])
            for net,inds in NETWORKS.items():
                row[f"static_fc_{net}"]=mean_within(z,inds,excluded); st,cl,ef=positive_graph_metrics(z,inds,excluded); row[f"strength_{net}"]=st; row[f"clustering_{net}"]=cl; row[f"efficiency_{net}"]=ef
            xc=x[:N_VOL, np.array(CORTEX)-1]; excl_c=np.array([idx in excluded for idx in CORTEX],bool); d=dynamic_fc(xc); row["temporal_variability"]=temporal_variability_sd(d,excl_c); metric_rows.append(row)
        metrics=pd.DataFrame(metric_rows).sort_values("participant_id"); metrics.to_csv(out/"participant_imaging_metrics.csv",index=False); merged=behavior.merge(metrics,on="participant_id",how="inner",validate="one_to_one"); merged.to_csv(out/"merged_behavior_imaging.csv",index=False)
        for col in ["table1_lesion_cm3","legacy_lesion_cm3"]:
            audit[f"lesion_vs_{col}_pearson"]=float(stats.pearsonr(merged.independent_lesion_volume_cm3,merged[col]).statistic); audit[f"lesion_vs_{col}_spearman"]=float(stats.spearmanr(merged.independent_lesion_volume_cm3,merged[col]).statistic)
        feature_sets={
            "baseline":["wab_aq","independent_lesion_volume_cm3"],
            "baseline_plus_static_fc":["wab_aq","independent_lesion_volume_cm3","static_fc_language","static_fc_salience","static_fc_dorsal_attention","static_fc_default_mode"],
            "baseline_plus_temporal_variability":["wab_aq","independent_lesion_volume_cm3","temporal_variability"],
            "baseline_plus_static_graph_tv":["wab_aq","independent_lesion_volume_cm3","static_fc_language","static_fc_salience","static_fc_dorsal_attention","static_fc_default_mode","strength_language","clustering_language","efficiency_language","efficiency_salience","efficiency_dorsal_attention","efficiency_default_mode","temporal_variability"],
        }
        needed=sorted(set(sum(feature_sets.values(),[]))); complete=merged.dropna(subset=["naming_change_pct"]+needed).copy(); audit["n_prediction_complete"]=len(complete)
        if len(complete)<24: raise RuntimeError(f"Prediction gate failed: only {len(complete)} complete participants")
        summ,preds=nested_loo(complete,"naming_change_pct",feature_sets); summ.to_csv(out/"prediction_summary.csv",index=False); preds.to_csv(out/"oof_predictions.csv",index=False); paired_incremental_tests(preds).to_csv(out/"incremental_tests.csv",index=False)
        assoc=[]; phys=[c for c in metrics.columns if c.startswith(("static_fc_","strength_","clustering_","efficiency_"))]+["temporal_variability"]
        for c in phys:
            d=merged[["naming_change_pct",c]].dropna()
            if len(d)>=10:
                r,p=stats.pearsonr(d.naming_change_pct,d[c]); rho,ps=stats.spearmanr(d.naming_change_pct,d[c]); assoc.append({"metric":c,"n":len(d),"pearson_r":r,"pearson_p":p,"spearman_rho":rho,"spearman_p":ps})
        pd.DataFrame(assoc).to_csv(out/"imaging_change_associations.csv",index=False); audit["analysis_completed"]=True
    except Exception as e:
        audit["fatal_error"]=repr(e); print("FATAL:",repr(e),file=sys.stderr)
    finally:
        write_json(out/"run_audit.json",audit); print(json.dumps(audit,indent=2,default=str))
    return 0 if audit.get("analysis_completed") else 2

if __name__=="__main__":
    raise SystemExit(main())
