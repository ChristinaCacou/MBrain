#!/usr/bin/env python3
"""Locked confirmatory analysis for the prospective Stentor perturbation study.

This script is intentionally narrow. It implements the pre-specified confirmatory
pipeline in Supplementary Protocol 1 and avoids exploratory model tuning.

Input
-----
metadata CSV with one row per recorded probe trial and columns:
cell_id, replicate_day, treatment_code, state_code, modality_code,
trial_id, trajectory_file, excluded, exclusion_reason

Each trajectory_file must be an .npz containing a 2-D array named `trajectory`
with shape (time, features). Segmentation/alignment occurs upstream according to
Protocol 1; the confirmatory script performs baseline subtraction and unit-energy
normalization before scale-insensitive metrics are calculated.

Outputs
-------
trial_metrics.csv, cell_state_metrics.csv, confirmatory_results.json,
primary_paired_differences.csv and confirmatory_summary.txt.

The primary confirmatory population is SHAM cells with both BASE and HAB
electrical-probe trajectories. Secondary analyses include cross-modal state
decoding and a frozen shared-state versus scalar-gain predictive log-score test.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.weightstats import ttost_paired

SEED = 260416
RNG = np.random.default_rng(SEED)
PRE_FRACTION = 0.20
BOOTSTRAPS = 5000
PERMUTATIONS = 2000
EQUIV_MARGIN_DZ = 0.20
PRIMARY_MIN_ABS_DZ = 0.35
ALLOWED_STATES = {"BASE", "HAB", "DISH", "REC"}
ALLOWED_MODALITIES = {"MECH", "ELEC"}
ALLOWED_TREATMENTS = {"SHAM", "RNAI"}
ALLOWED_EXCLUSIONS = {
    "", "none", "detached", "left_field", "cell_death", "segmentation_failure",
    "stimulus_delivery_failure", "predefined_qc_failure"
}
REQUIRED = [
    "cell_id", "replicate_day", "treatment_code", "state_code", "modality_code",
    "trial_id", "trajectory_file", "excluded", "exclusion_reason"
]


def _bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"1", "true", "t", "yes", "y"}: return True
    if s in {"0", "false", "f", "no", "n", ""}: return False
    raise ValueError(f"Cannot parse boolean value: {v!r}")


def validate_metadata(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    out = df.copy()
    out["excluded"] = out["excluded"].map(_bool)
    for c in ["cell_id", "replicate_day", "treatment_code", "state_code", "modality_code", "trial_id"]:
        out[c] = out[c].astype(str).str.strip()
    out["treatment_code"] = out["treatment_code"].str.upper()
    out["state_code"] = out["state_code"].str.upper()
    out["modality_code"] = out["modality_code"].str.upper()
    bad_state = sorted(set(out.state_code) - ALLOWED_STATES)
    bad_mod = sorted(set(out.modality_code) - ALLOWED_MODALITIES)
    bad_treat = sorted(set(out.treatment_code) - ALLOWED_TREATMENTS)
    if bad_state: raise ValueError(f"Unexpected state codes: {bad_state}")
    if bad_mod: raise ValueError(f"Unexpected modality codes: {bad_mod}")
    if bad_treat: raise ValueError(f"Unexpected treatment codes: {bad_treat}")
    reasons = out["exclusion_reason"].fillna("").astype(str).str.strip().str.lower()
    bad_reason = sorted(set(reasons[out.excluded]) - ALLOWED_EXCLUSIONS)
    if bad_reason:
        raise ValueError(f"Excluded rows contain non-preregistered reasons: {bad_reason}")
    if out.loc[~out.excluded, "trajectory_file"].isna().any():
        raise ValueError("Included rows may not have missing trajectory_file values")
    dup = out.duplicated(["cell_id", "trial_id"], keep=False)
    if dup.any():
        raise ValueError("Duplicate cell_id/trial_id rows detected")
    return out


def effective_dimensionality(x: np.ndarray) -> float:
    if x.ndim != 2 or x.shape[0] < 4:
        return np.nan
    cov = np.cov(x, rowvar=False)
    vals = np.linalg.eigvalsh(np.atleast_2d(cov))
    vals = vals[vals > 1e-12]
    if len(vals) == 0:
        return np.nan
    return float(vals.sum() ** 2 / np.square(vals).sum())


def lz76_complexity(bits: np.ndarray) -> float:
    """Normalized Lempel-Ziv phrase complexity for a binary vector.

    The implementation uses a deterministic incremental phrase dictionary. It is
    intentionally fixed for preregistration and is used only for relative,
    scale-insensitive comparisons across identically processed trials.
    """
    s = ''.join('1' if b else '0' for b in np.asarray(bits, dtype=bool).ravel())
    n = len(s)
    if n < 2:
        return np.nan
    phrases = set()
    i = 0
    c = 0
    while i < n:
        j = i + 1
        while j <= n and s[i:j] in phrases:
            j += 1
        phrases.add(s[i:min(j, n)])
        c += 1
        i = j if j <= n else n
    # Standard phrase-count normalization scale n/log2(n).
    norm = n / max(math.log2(n), 1.0)
    return float(c / norm)

def compressed_complexity(bits: np.ndarray) -> float:
    arr = np.packbits(np.asarray(bits, dtype=np.uint8).ravel())
    if arr.size == 0:
        return np.nan
    comp = zlib.compress(arr.tobytes(), level=9)
    # Normalize by byte length; fixed protocol comparison makes header cost common.
    return float(len(comp) / max(len(arr.tobytes()), 1))


def trajectory_metrics(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or min(x.shape) < 2:
        raise ValueError(f"trajectory must be 2-D time x features, got {x.shape}")
    if not np.isfinite(x).all():
        raise ValueError("trajectory contains NaN or infinite values")
    npre = max(3, int(round(PRE_FRACTION * x.shape[0])))
    baseline = x[:npre].mean(axis=0, keepdims=True)
    y = x - baseline
    energy = float(np.square(y).sum())
    if energy <= 1e-12:
        raise ValueError("trajectory has near-zero response energy")
    yn = y / math.sqrt(energy)
    # Scale-insensitive metrics use unit-energy trajectory.
    dim = effective_dimensionality(yn)
    med = np.median(yn, axis=0, keepdims=True)
    bits = yn > med
    return {
        "response_energy": energy,
        "effective_dimensionality": dim,
        "lz76_norm": lz76_complexity(bits),
        "compressed_complexity": compressed_complexity(bits),
    }


def load_trial_metrics(meta: pd.DataFrame, trajectory_dir: Path) -> pd.DataFrame:
    rows = []
    for r in meta.itertuples(index=False):
        if r.excluded:
            continue
        p = Path(r.trajectory_file)
        if not p.is_absolute():
            p = trajectory_dir / p
        if not p.exists():
            raise FileNotFoundError(f"Missing included trajectory: {p}")
        with np.load(p) as z:
            if "trajectory" not in z:
                raise ValueError(f"{p} lacks array named 'trajectory'")
            metrics = trajectory_metrics(z["trajectory"])
        rows.append({
            "cell_id": r.cell_id, "replicate_day": r.replicate_day,
            "treatment_code": r.treatment_code, "state_code": r.state_code,
            "modality_code": r.modality_code, "trial_id": r.trial_id,
            **metrics,
        })
    if not rows:
        raise ValueError("No included trials remain after pre-specified exclusions")
    return pd.DataFrame(rows)


def aggregate_cell_state(trials: pd.DataFrame) -> pd.DataFrame:
    group = ["cell_id", "replicate_day", "treatment_code", "state_code", "modality_code"]
    metrics = ["response_energy", "effective_dimensionality", "lz76_norm", "compressed_complexity"]
    return trials.groupby(group, as_index=False)[metrics].median()


def bootstrap_mean_ci(values: np.ndarray, clusters: np.ndarray | None = None, B: int = BOOTSTRAPS):
    v = np.asarray(values, float)
    if clusters is None:
        sims = [np.mean(RNG.choice(v, size=len(v), replace=True)) for _ in range(B)]
    else:
        c = np.asarray(clusters)
        uniq = np.unique(c)
        sims = []
        for _ in range(B):
            sampled = RNG.choice(uniq, size=len(uniq), replace=True)
            vv = np.concatenate([v[c == x] for x in sampled])
            sims.append(np.mean(vv))
    return np.quantile(sims, [0.025, 0.975]).tolist()


def paired_primary(cell: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    d = cell[(cell.treatment_code == "SHAM") & (cell.modality_code == "ELEC")]
    piv = d.pivot_table(index=["cell_id", "replicate_day"], columns="state_code",
                        values="effective_dimensionality", aggfunc="median")
    piv = piv.dropna(subset=["BASE", "HAB"]).reset_index()
    if len(piv) < 8:
        raise ValueError("Primary analysis requires at least 8 paired SHAM cells with BASE and HAB ELEC data")
    diff = (piv["HAB"] - piv["BASE"]).to_numpy(float)
    sd = diff.std(ddof=1)
    dz = float(diff.mean() / sd) if sd > 0 else (np.inf if diff.mean() != 0 else 0.0)
    t = stats.ttest_rel(piv["HAB"], piv["BASE"])
    try:
        w = stats.wilcoxon(piv["HAB"], piv["BASE"], zero_method="wilcox", alternative="two-sided")
        wp = float(w.pvalue)
    except ValueError:
        wp = 1.0
    pooled_ref = float(np.std(np.concatenate([piv["BASE"], piv["HAB"]]), ddof=1))
    margin = EQUIV_MARGIN_DZ * pooled_ref
    if margin > 0:
        tost = ttost_paired(piv["HAB"].to_numpy(), piv["BASE"].to_numpy(), -margin, margin)
        tost_p = float(tost[0])
    else:
        tost_p = np.nan
    ci = bootstrap_mean_ci(diff, piv["replicate_day"].to_numpy())
    res = {
        "n_paired": int(len(piv)),
        "baseline_mean": float(piv.BASE.mean()),
        "habituated_mean": float(piv.HAB.mean()),
        "mean_paired_difference_hab_minus_base": float(diff.mean()),
        "cluster_bootstrap_95ci": ci,
        "paired_dz": dz,
        "paired_t_p": float(t.pvalue),
        "wilcoxon_p": wp,
        "equivalence_margin_dz": EQUIV_MARGIN_DZ,
        "tost_equivalence_p": tost_p,
        "min_effect_for_null_rejection_abs_dz": PRIMARY_MIN_ABS_DZ,
        "passes_primary_difference_gate": bool(t.pvalue < 0.05 and abs(dz) >= PRIMARY_MIN_ABS_DZ),
        "supports_preregistered_equivalence": bool(tost_p < 0.05) if np.isfinite(tost_p) else False,
    }
    out = piv.copy(); out["difference_hab_minus_base"] = diff
    return res, out


def feature_table(cell: pd.DataFrame) -> pd.DataFrame:
    feats = ["effective_dimensionality", "lz76_norm", "compressed_complexity"]
    x = cell.copy()
    # response energy intentionally excluded from state classifier/model score.
    return x[["cell_id", "replicate_day", "treatment_code", "state_code", "modality_code", *feats]]


def cross_modal_decode(cell: pd.DataFrame, permutations: int = PERMUTATIONS) -> dict:
    d = feature_table(cell)
    d = d[(d.treatment_code == "SHAM") & (d.state_code.isin(["BASE", "HAB"]))]
    mech = d[d.modality_code == "MECH"].copy()
    elec = d[d.modality_code == "ELEC"].copy()
    feats = ["effective_dimensionality", "lz76_norm", "compressed_complexity"]
    common = sorted(set(mech.cell_id) & set(elec.cell_id))
    mech = mech[mech.cell_id.isin(common)]
    elec = elec[elec.cell_id.isin(common)]
    mech = mech.dropna(subset=feats); elec = elec.dropna(subset=feats)
    if mech.state_code.nunique() < 2 or elec.state_code.nunique() < 2 or len(elec) < 16:
        return {"available": False, "reason": "insufficient paired multimodal cells"}
    scaler = StandardScaler().fit(mech[feats])
    clf = LogisticRegression(C=1.0, solver="liblinear", random_state=SEED, max_iter=1000)
    ytr = (mech.state_code == "HAB").astype(int).to_numpy()
    clf.fit(scaler.transform(mech[feats]), ytr)
    y = (elec.state_code == "HAB").astype(int).to_numpy()
    pred = clf.predict_proba(scaler.transform(elec[feats]))[:, 1]
    auc = float(roc_auc_score(y, pred))
    # Permute state labels at cell level to preserve within-cell dependence.
    cell_labels = elec[["cell_id", "state_code"]].drop_duplicates()
    cells = np.array(sorted(cell_labels.cell_id.unique()))
    null = []
    for _ in range(permutations):
        flip = dict(zip(cells, RNG.integers(0, 2, size=len(cells))))
        yp = np.array([yy ^ flip[c] for c, yy in zip(elec.cell_id, y)])
        if len(np.unique(yp)) < 2:
            continue
        null.append(roc_auc_score(yp, pred))
    p = (1 + np.sum(np.asarray(null) >= auc)) / (1 + len(null))
    # Cell-cluster bootstrap CI.
    aucs = []
    for _ in range(2000):
        sample = RNG.choice(cells, size=len(cells), replace=True)
        idx = np.concatenate([np.flatnonzero(elec.cell_id.to_numpy() == c) for c in sample])
        yy, pp = y[idx], pred[idx]
        if len(np.unique(yy)) == 2:
            aucs.append(roc_auc_score(yy, pp))
    ci = np.quantile(aucs, [0.025, 0.975]).tolist() if aucs else [np.nan, np.nan]
    return {
        "available": True, "n_cells": int(len(cells)), "auc": auc,
        "cluster_bootstrap_95ci": ci, "permutation_p": float(p),
        "passes_cross_modal_gate": bool(p < 0.05 and ci[0] > 0.55),
    }


def mvn_logpdf(x, mu, cov):
    x = np.asarray(x, float); mu = np.asarray(mu, float)
    cov = np.asarray(cov, float) + np.eye(len(mu)) * 1e-6
    return float(stats.multivariate_normal.logpdf(x, mean=mu, cov=cov, allow_singular=True))


def frozen_model_score(cell: pd.DataFrame, permutations: int = PERMUTATIONS) -> dict:
    """Leave-replicate-day-out score: scalar-gain N0 vs shared-state-shift N1.

    N0 predicts no state change in unit-energy electrical response geometry.
    N1 predicts that the mean BASE->HAB shift learned from normalized mechanical
    trajectories transfers to electrical trajectories. No electrical HAB data from
    a held-out day are used to construct either prediction.
    """
    d = feature_table(cell)
    d = d[(d.treatment_code == "SHAM") & d.state_code.isin(["BASE", "HAB"])]
    feats = ["effective_dimensionality", "lz76_norm", "compressed_complexity"]
    days = sorted(d.replicate_day.unique())
    records = []
    for day in days:
        train = d[d.replicate_day != day]
        test = d[d.replicate_day == day]
        mp = train.pivot_table(index="cell_id", columns=["modality_code", "state_code"], values=feats)
        # Pandas layout is metric first; extract robustly through per-cell merge instead.
        m0 = train[(train.modality_code == "MECH") & (train.state_code == "BASE")].set_index("cell_id")[feats]
        m1 = train[(train.modality_code == "MECH") & (train.state_code == "HAB")].set_index("cell_id")[feats]
        ids = m0.index.intersection(m1.index)
        if len(ids) < 6:
            continue
        delta = (m1.loc[ids] - m0.loc[ids]).mean(axis=0).to_numpy()
        resid = (m1.loc[ids] - m0.loc[ids]).to_numpy() - delta
        cov = np.cov(resid, rowvar=False) if resid.shape[0] > 3 else np.eye(len(feats)) * 0.05
        e0 = test[(test.modality_code == "ELEC") & (test.state_code == "BASE")].set_index("cell_id")[feats]
        e1 = test[(test.modality_code == "ELEC") & (test.state_code == "HAB")].set_index("cell_id")[feats]
        tids = e0.index.intersection(e1.index)
        for cid in tids:
            base = e0.loc[cid].to_numpy(float)
            obs = e1.loc[cid].to_numpy(float)
            l0 = mvn_logpdf(obs, base, cov)
            l1 = mvn_logpdf(obs, base + delta, cov)
            records.append((cid, day, l1-l0))
    if len(records) < 12:
        return {"available": False, "reason": "insufficient leave-day-out multimodal pairs"}
    rr = pd.DataFrame(records, columns=["cell_id", "replicate_day", "delta_log_score_n1_minus_n0"])
    mean_d = float(rr.delta_log_score_n1_minus_n0.mean())
    # Sign-flip at the cell level is the frozen null for no model advantage.
    vals = rr.delta_log_score_n1_minus_n0.to_numpy(float)
    null = []
    for _ in range(permutations):
        signs = RNG.choice([-1, 1], size=len(vals))
        null.append(float(np.mean(vals * signs)))
    p = (1 + np.sum(np.asarray(null) >= mean_d)) / (1 + len(null))
    ci = bootstrap_mean_ci(vals, rr.replicate_day.to_numpy())
    return {
        "available": True, "n_cells": int(len(rr)), "mean_delta_log_score": mean_d,
        "cluster_bootstrap_95ci": ci, "signflip_p": float(p),
        "passes_model_score_gate": bool(mean_d > 0 and p < 0.05 and ci[0] > 0),
    }


def treatment_secondary(cell: pd.DataFrame) -> dict:
    # Secondary: compare BASE->HAB electrical dimensionality shifts between SHAM and RNAI.
    d = cell[(cell.modality_code == "ELEC") & cell.state_code.isin(["BASE", "HAB"])]
    p = d.pivot_table(index=["cell_id", "replicate_day", "treatment_code"], columns="state_code",
                      values="effective_dimensionality").dropna(subset=["BASE", "HAB"]).reset_index()
    p["delta"] = p.HAB - p.BASE
    groups = {k:g.delta.to_numpy(float) for k,g in p.groupby("treatment_code")}
    if not {"SHAM", "RNAI"}.issubset(groups) or min(len(groups["SHAM"]), len(groups["RNAI"])) < 8:
        return {"available": False, "reason": "insufficient treatment-arm pairs"}
    u = stats.mannwhitneyu(groups["SHAM"], groups["RNAI"], alternative="two-sided")
    return {
        "available": True, "n_sham": int(len(groups["SHAM"])), "n_rnai": int(len(groups["RNAI"])),
        "mean_delta_sham": float(groups["SHAM"].mean()), "mean_delta_rnai": float(groups["RNAI"].mean()),
        "mannwhitney_p": float(u.pvalue),
    }


def run_pipeline(metadata: Path, trajectory_dir: Path, outdir: Path, permutations=PERMUTATIONS):
    outdir.mkdir(parents=True, exist_ok=True)
    meta = validate_metadata(pd.read_csv(metadata))
    trials = load_trial_metrics(meta, trajectory_dir)
    cells = aggregate_cell_state(trials)
    primary, paired = paired_primary(cells)
    cross = cross_modal_decode(cells, permutations=permutations)
    score = frozen_model_score(cells, permutations=permutations)
    treatment = treatment_secondary(cells)
    overall = {
        "analysis_version": "1.0-locked",
        "seed": SEED,
        "primary": primary,
        "cross_modal": cross,
        "model_score": score,
        "treatment_secondary": treatment,
    }
    overall["decision"] = {
        "local_scalar_gain_null_rejected": bool(
            primary.get("passes_primary_difference_gate", False)
            and cross.get("passes_cross_modal_gate", False)
            and score.get("passes_model_score_gate", False)
        ),
        "interpretation_rule": (
            "The global integrative-state alternative is favored over the pre-specified scalar-gain null only "
            "if all three confirmatory gates pass. This result is not, by itself, evidence sufficient to establish consciousness."
        )
    }
    trials.to_csv(outdir / "trial_metrics.csv", index=False)
    cells.to_csv(outdir / "cell_state_metrics.csv", index=False)
    paired.to_csv(outdir / "primary_paired_differences.csv", index=False)
    with open(outdir / "confirmatory_results.json", "w") as f:
        json.dump(overall, f, indent=2)
    with open(outdir / "confirmatory_summary.txt", "w") as f:
        f.write("LOCKED PROSPECTIVE STENTOR CONFIRMATORY ANALYSIS\n")
        f.write("================================================\n")
        f.write(json.dumps(overall, indent=2))
        f.write("\n")
    return overall


def _synthetic_trajectory(rng, state, modality, integrative=True, ntime=100, nfeat=32):
    t = np.linspace(0, 1, ntime)
    # Shared response family with modality-specific temporal basis.
    basis = np.vstack([
        np.exp(-((t-0.35)/(0.09 if modality=='MECH' else 0.12))**2),
        np.sin(2*np.pi*t)*np.exp(-2*t),
        np.cos(3*np.pi*t)*np.exp(-3*t),
    ]).T
    W = rng.normal(size=(3,nfeat))
    x = basis @ W
    x += rng.normal(scale=0.08, size=x.shape)
    amp = 1.0 if state == 'BASE' else 0.55
    x *= amp
    if integrative and state == 'HAB':
        # Geometry change not removable by amplitude normalization.
        latent = np.outer(np.exp(-((t-0.55)/0.16)**2), np.linspace(-1,1,nfeat))
        x += 0.28 * latent
    return x


def self_test(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    tdir = outdir / "synthetic_trajectories"; tdir.mkdir(exist_ok=True)
    rows=[]; rng=np.random.default_rng(SEED)
    for treatment in ['SHAM','RNAI']:
        for day in range(1,5):
            for j in range(12):
                cid=f"{treatment}_D{day}_{j:02d}"
                for state in ['BASE','HAB']:
                    for mod in ['MECH','ELEC']:
                        arr=_synthetic_trajectory(rng,state,mod,integrative=True)
                        if treatment=='RNAI' and state=='HAB':
                            arr += 0.05*np.outer(np.sin(np.pi*np.linspace(0,1,arr.shape[0])), np.ones(arr.shape[1]))
                        fn=f"{cid}_{state}_{mod}.npz"
                        np.savez_compressed(tdir/fn, trajectory=arr)
                        rows.append(dict(cell_id=cid,replicate_day=f"D{day}",treatment_code=treatment,
                                         state_code=state,modality_code=mod,trial_id=f"{state}_{mod}",
                                         trajectory_file=fn,excluded=False,exclusion_reason=""))
    meta=outdir/'synthetic_metadata.csv'; pd.DataFrame(rows).to_csv(meta,index=False)
    return run_pipeline(meta,tdir,outdir/'synthetic_results',permutations=100)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--metadata', type=Path)
    ap.add_argument('--trajectory-dir', type=Path, default=Path('.'))
    ap.add_argument('--outdir', type=Path, default=Path('stentor_confirmatory_results'))
    ap.add_argument('--self-test', action='store_true')
    args=ap.parse_args()
    if args.self_test:
        result=self_test(args.outdir)
    else:
        if args.metadata is None:
            ap.error('--metadata is required unless --self-test is used')
        result=run_pipeline(args.metadata,args.trajectory_dir,args.outdir)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
