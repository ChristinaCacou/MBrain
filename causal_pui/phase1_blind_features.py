import hashlib
import json
import os
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, log_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

SEED = 20261011
N_OUTER = 4
N_INNER = 3
N_DRAWS = 3
FIXED_NEURONS = 5
TARGET_PCS = 2
OTHER_PCS = 1
WINDOW_START = 0.0
WINDOW_END = 0.2
CS = (0.03, 0.1, 0.3, 1.0, 3.0)
MIN_RESPONDED_TRIALS = 80
MIN_CLASS = 20

DATA_URLS = {
    "steinmetz_part0.npz": "https://raw.githubusercontent.com/tipsyfermion/NMA_Mapping_Brain_Networks_2020/master/code/steinmetz_part0.npz",
    "steinmetz_part1.npz": "https://raw.githubusercontent.com/tipsyfermion/NMA_Mapping_Brain_Networks_2020/master/code/steinmetz_part1.npz",
    "steinmetz_part2.npz": "https://raw.githubusercontent.com/tipsyfermion/NMA_Mapping_Brain_Networks_2020/master/code/steinmetz_part2.npz",
}

OUTDIR = Path("analysis/pui_blind")
CACHEDIR = Path(".cache/steinmetz")
OUTDIR.mkdir(parents=True, exist_ok=True)
CACHEDIR.mkdir(parents=True, exist_ok=True)

# Purely anatomical aliases. No intervention effect values are present here.
def alias_area(area):
    a = str(area).strip()
    if a.startswith("SSp"):
        return "SSp"
    if a.startswith("ACA"):
        return "ACA"
    if a.startswith("ORB"):
        return "ORB"
    return a


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_sessions():
    sessions = []
    checksums = {}
    for name, url in DATA_URLS.items():
        path = CACHEDIR / name
        if not path.exists():
            print("downloading", name, flush=True)
            urlretrieve(url, path)
        checksums[name] = sha256(path)
        with np.load(path, allow_pickle=True) as z:
            part = z["dat"]
            sessions.extend(part.tolist())
        print(name, "sessions", len(part), flush=True)
    if len(sessions) != 39:
        raise RuntimeError(f"Expected 39 sessions, got {len(sessions)}")
    return sessions, checksums


def activity_window(dat):
    spks = np.asarray(dat["spks"])
    dt = float(np.asarray(dat.get("bin_size", 0.01)).reshape(-1)[0])
    stim_onset = float(np.asarray(dat.get("stim_onset", 0.5)).reshape(-1)[0])
    i0 = int(round((stim_onset + WINDOW_START) / dt))
    i1 = int(round((stim_onset + WINDOW_END) / dt))
    i0 = max(0, min(i0, spks.shape[2] - 1))
    i1 = max(i0 + 1, min(i1, spks.shape[2]))
    return slice(i0, i1), dt


def history_matrix(dat, n_trials):
    cl = np.asarray(dat["contrast_left"], dtype=float).reshape(-1)
    cr = np.asarray(dat["contrast_right"], dtype=float).reshape(-1)
    response = np.asarray(dat["response"], dtype=float).reshape(-1)
    feedback = np.asarray(dat.get("feedback_type", np.zeros(n_trials)), dtype=float).reshape(-1)
    if not (len(cl) == len(cr) == len(response) == n_trials):
        raise ValueError("trial covariates do not align")
    signed = cr - cl
    total = cr + cl
    H = np.column_stack([cl, cr, signed, total])
    prev = np.zeros((n_trials, 4), dtype=float)
    prev[1:, 0] = response[:-1]
    prev[1:, 1] = feedback[:-1] if len(feedback) == n_trials else 0
    prev[1:, 2] = cl[:-1]
    prev[1:, 3] = cr[:-1]
    return np.column_stack([H, prev])


def eligible_trials(dat, n_trials):
    response = np.asarray(dat["response"], dtype=float).reshape(-1)
    ok = np.isfinite(response) & (response != 0)
    try:
        rt = np.asarray(dat["response_time"], dtype=float).reshape(-1)
        if len(rt) == n_trials:
            ok &= np.isfinite(rt) & (rt > WINDOW_END)
    except Exception:
        pass
    idx = np.flatnonzero(ok)
    y_raw = response[idx]
    vals = np.unique(y_raw)
    if vals.size != 2:
        return np.array([], dtype=int), np.array([], dtype=int)
    y = (y_raw == vals.max()).astype(int)
    return idx, y


def fit_transform_region(X_all, train_idx, test_idx, n_pcs):
    scaler = StandardScaler()
    Xtr = scaler.fit_transform(X_all[train_idx])
    Xte = scaler.transform(X_all[test_idx])
    ncomp = min(n_pcs, Xtr.shape[1], max(1, Xtr.shape[0] - 1))
    pca = PCA(n_components=ncomp, svd_solver="full")
    return pca.fit_transform(Xtr), pca.transform(Xte)


def choose_c_and_fit(Xtr, ytr):
    counts = np.bincount(ytr, minlength=2)
    nsplit = min(N_INNER, int(counts.min()))
    if nsplit < 2:
        return None
    cv = StratifiedKFold(n_splits=nsplit, shuffle=True, random_state=SEED + len(ytr) + Xtr.shape[1])
    best_c, best_loss = None, np.inf
    for c in CS:
        losses = []
        for a, b in cv.split(Xtr, ytr):
            ss = StandardScaler()
            A = ss.fit_transform(Xtr[a])
            B = ss.transform(Xtr[b])
            model = LogisticRegression(C=c, solver="lbfgs", max_iter=1000, random_state=SEED)
            model.fit(A, ytr[a])
            p = model.predict_proba(B)[:, 1]
            losses.append(log_loss(ytr[b], p, labels=[0, 1]))
        m = float(np.mean(losses))
        if m < best_loss:
            best_loss, best_c = m, c
    ss = StandardScaler()
    A = ss.fit_transform(Xtr)
    model = LogisticRegression(C=best_c, solver="lbfgs", max_iter=1000, random_state=SEED)
    model.fit(A, ytr)
    return ss, model, best_c


def evaluate_model(Xtr, ytr, Xte, yte):
    fit = choose_c_and_fit(Xtr, ytr)
    if fit is None:
        return np.nan, np.nan, np.nan
    ss, model, _ = fit
    B = ss.transform(Xte)
    p = model.predict_proba(B)[:, 1]
    pred = (p >= 0.5).astype(int)
    loss = log_loss(yte, p, labels=[0, 1])
    bacc = balanced_accuracy_score(yte, pred)
    auc = roc_auc_score(yte, p) if np.unique(yte).size == 2 else np.nan
    return float(loss), float(bacc), float(auc)


def condition_snr(signal, cl, cr):
    keys = list(zip(cl.tolist(), cr.tolist()))
    grand = float(np.mean(signal))
    means, within = [], []
    for key in sorted(set(keys)):
        ii = np.array([i for i, k in enumerate(keys) if k == key], dtype=int)
        if len(ii) < 2:
            continue
        means.append(float(np.mean(signal[ii])))
        within.append(float(np.var(signal[ii], ddof=1)))
    if len(means) < 2:
        return np.nan
    between = float(np.var(means, ddof=1))
    w = float(np.mean(within)) if within else np.nan
    return between / (w + 1e-12) if np.isfinite(w) else np.nan


def session_rows(session_index, dat):
    spks = np.asarray(dat["spks"], dtype=np.float32)
    n_neurons, n_trials, _ = spks.shape
    idx, y = eligible_trials(dat, n_trials)
    if len(idx) < MIN_RESPONDED_TRIALS or min(np.bincount(y, minlength=2)) < MIN_CLASS:
        return []
    H_all = history_matrix(dat, n_trials)[idx]
    win, dt = activity_window(dat)
    areas_raw = np.asarray(dat["brain_area"]).astype(str)
    aliases = np.asarray([alias_area(a) for a in areas_raw])
    invalid = {"root", "nan", "", "none"}
    region_neurons = {
        a: np.flatnonzero(aliases == a)
        for a in sorted(set(aliases))
        if a.lower() not in invalid and np.sum(aliases == a) >= FIXED_NEURONS
    }
    if len(region_neurons) < 2:
        return []
    # Sum spikes in the fixed pre-response stimulus window.
    region_activity = {
        a: spks[nids][:, idx, win].sum(axis=2).T.astype(np.float64)
        for a, nids in region_neurons.items()
    }
    cl = np.asarray(dat["contrast_left"], dtype=float).reshape(-1)[idx]
    cr = np.asarray(dat["contrast_right"], dtype=float).reshape(-1)[idx]
    mouse = str(np.asarray(dat.get("mouse_name", f"session{session_index}")).reshape(-1)[0])
    date = str(np.asarray(dat.get("date_exp", "")).reshape(-1)[0])
    outer = StratifiedKFold(n_splits=N_OUTER, shuffle=True, random_state=SEED + session_index)
    folds = list(outer.split(H_all, y))
    results = []
    for draw in range(N_DRAWS):
        rng = np.random.default_rng(SEED + session_index * 1009 + draw * 9176)
        sampled = {}
        for a, nids in region_neurons.items():
            local = rng.choice(len(nids), size=FIXED_NEURONS, replace=False)
            sampled[a] = region_activity[a][:, local]
        for target in sorted(sampled):
            fold_stats = []
            for fold, (tr, te) in enumerate(folds):
                # Build target PCs and one PC from every simultaneously observed other region.
                Ttr, Tte = fit_transform_region(sampled[target], tr, te, TARGET_PCS)
                other_tr, other_te = [], []
                for other in sorted(sampled):
                    if other == target:
                        continue
                    Otr, Ote = fit_transform_region(sampled[other], tr, te, OTHER_PCS)
                    other_tr.append(Otr)
                    other_te.append(Ote)
                Otr = np.column_stack(other_tr) if other_tr else np.empty((len(tr), 0))
                Ote = np.column_stack(other_te) if other_te else np.empty((len(te), 0))
                htr, hte = H_all[tr], H_all[te]
                base_tr, base_te = np.column_stack([htr, Otr]), np.column_stack([hte, Ote])
                full_tr, full_te = np.column_stack([htr, Otr, Ttr]), np.column_stack([hte, Ote, Tte])
                marg_tr, marg_te = np.column_stack([htr, Ttr]), np.column_stack([hte, Tte])
                target_tr, target_te = Ttr, Tte
                h_loss, _, _ = evaluate_model(htr, y[tr], hte, y[te])
                b_loss, _, _ = evaluate_model(base_tr, y[tr], base_te, y[te])
                f_loss, _, _ = evaluate_model(full_tr, y[tr], full_te, y[te])
                m_loss, _, _ = evaluate_model(marg_tr, y[tr], marg_te, y[te])
                t_loss, t_bacc, t_auc = evaluate_model(target_tr, y[tr], target_te, y[te])
                fold_stats.append((h_loss, b_loss, f_loss, m_loss, t_loss, t_bacc, t_auc))
            arr = np.asarray(fold_stats, dtype=float)
            if not np.all(np.isfinite(arr[:, :5])):
                continue
            mean = np.nanmean(arr, axis=0)
            # population quality controls are computed without behavior labels
            raw_signal = sampled[target].mean(axis=1)
            snr = condition_snr(raw_signal, cl, cr)
            results.append({
                "session": session_index,
                "mouse": mouse,
                "date": date,
                "region": target,
                "draw": draw,
                "n_trials": len(idx),
                "n_regions_conditioned": len(sampled) - 1,
                "n_neurons_recorded": len(region_neurons[target]),
                "n_neurons_matched": FIXED_NEURONS,
                "bin_size_s": dt,
                "window_start_s": WINDOW_START,
                "window_end_s": WINDOW_END,
                "pui_nats_per_trial": float(mean[1] - mean[2]),
                "marginal_gain_nats_per_trial": float(mean[0] - mean[3]),
                "target_only_loss": float(mean[4]),
                "target_only_bacc": float(mean[5]),
                "target_only_auc": float(mean[6]),
                "condition_snr": float(snr) if np.isfinite(snr) else np.nan,
                "mean_window_spikes": float(np.mean(raw_signal)),
            })
    return results


def main():
    sessions, raw_checksums = load_sessions()
    rows = []
    for i, dat in enumerate(sessions):
        out = session_rows(i, dat)
        rows.extend(out)
        print(f"session {i+1}/39: {len(out)} feature rows", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No feature rows produced")
    df = df.sort_values(["region", "session", "draw"]).reset_index(drop=True)
    session_summary = (df.groupby(["region", "session"], as_index=False)
        .agg(pui_nats_per_trial=("pui_nats_per_trial", "mean"),
             marginal_gain_nats_per_trial=("marginal_gain_nats_per_trial", "mean"),
             target_only_bacc=("target_only_bacc", "mean"),
             target_only_auc=("target_only_auc", "mean"),
             condition_snr=("condition_snr", "mean"),
             mean_window_spikes=("mean_window_spikes", "mean"),
             n_neurons_recorded=("n_neurons_recorded", "first"),
             n_trials=("n_trials", "first"),
             n_regions_conditioned=("n_regions_conditioned", "first"),
             mouse=("mouse", "first")))
    region_summary = (session_summary.groupby("region", as_index=False)
        .agg(pui_mean=("pui_nats_per_trial", "mean"),
             pui_median=("pui_nats_per_trial", "median"),
             pui_sd=("pui_nats_per_trial", "std"),
             marginal_mean=("marginal_gain_nats_per_trial", "mean"),
             local_bacc_mean=("target_only_bacc", "mean"),
             local_auc_mean=("target_only_auc", "mean"),
             snr_mean=("condition_snr", "mean"),
             mean_window_spikes=("mean_window_spikes", "mean"),
             neurons_recorded_mean=("n_neurons_recorded", "mean"),
             conditioned_regions_mean=("n_regions_conditioned", "mean"),
             n_sessions=("session", "nunique"),
             n_mice=("mouse", "nunique")))
    df.to_csv(OUTDIR / "blind_features_draws.csv", index=False)
    session_summary.to_csv(OUTDIR / "blind_features_session.csv", index=False)
    region_summary.to_csv(OUTDIR / "blind_features_region.csv", index=False)
    config = {
        "seed": SEED,
        "outcome_blind": True,
        "perturbation_labels_loaded": False,
        "primary_metric": "pui_mean",
        "target": "responded-trial choice direction",
        "window_s": [WINDOW_START, WINDOW_END],
        "history_covariates": ["contrast_left", "contrast_right", "signed_evidence", "total_contrast", "prev_response", "prev_feedback", "prev_contrast_left", "prev_contrast_right"],
        "fixed_neurons_per_region": FIXED_NEURONS,
        "target_pcs": TARGET_PCS,
        "other_region_pcs": OTHER_PCS,
        "outer_folds": N_OUTER,
        "inner_folds": N_INNER,
        "draws": N_DRAWS,
        "C_grid": list(CS),
        "raw_sha256": raw_checksums,
        "script_sha256": sha256(Path(__file__)),
    }
    lock_files = [OUTDIR / "blind_features_draws.csv", OUTDIR / "blind_features_session.csv", OUTDIR / "blind_features_region.csv"]
    config["feature_sha256"] = {p.name: sha256(p) for p in lock_files}
    with open(OUTDIR / "blind_lock.json", "w") as f:
        json.dump(config, f, indent=2, sort_keys=True)
    report = [
        "# Blinded neural feature lock",
        "",
        "Perturbation outcomes were not present in this phase or script.",
        f"Sessions loaded: {len(sessions)}",
        f"Draw rows: {len(df)}",
        f"Session-region rows: {len(session_summary)}",
        f"Regions: {len(region_summary)}",
        "",
        "Primary frozen predictor: `pui_mean` = held-out base log-loss minus held-out full log-loss, in nats/trial.",
        "Positive values mean the target region adds predictive information after conditioning on current stimulus, trial history, and other simultaneously recorded regions.",
        "",
        "## SHA-256 lock",
        "```json",
        json.dumps(config["feature_sha256"], indent=2, sort_keys=True),
        "```",
    ]
    (OUTDIR / "BLIND_LOCK.md").write_text("\n".join(report))
    print(region_summary.sort_values("pui_mean", ascending=False).head(20).to_string(index=False), flush=True)

if __name__ == "__main__":
    main()
