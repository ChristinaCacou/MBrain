import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

SEED = 20261011
N_SPLITS = 3000
ROOT = Path('analysis/pui_blind_fixedC')
INFILE = ROOT / 'blind_features_session.csv'
OUTDIR = ROOT / 'diagnostics'
OUTDIR.mkdir(parents=True, exist_ok=True)

METRICS = {
    'pui_nats_per_trial': 'conditional predictive gain (locked PUI variable)',
    'marginal_gain_nats_per_trial': 'marginal predictive gain',
    'target_only_auc': 'local decoder AUC',
    'target_only_bacc': 'local decoder balanced accuracy',
}


def fast_spearman(x, y):
    rx = rankdata(np.asarray(x, float))
    ry = rankdata(np.asarray(y, float))
    rx -= rx.mean(); ry -= ry.mean()
    den = np.sqrt((rx @ rx) * (ry @ ry))
    return float((rx @ ry) / den) if den > 0 else np.nan


def split_half(df, unit_col, min_units=4):
    rng = np.random.default_rng(SEED + (0 if unit_col == 'session' else 1000))
    base = df.groupby(['region', unit_col], as_index=False).agg(
        **{m: (m, 'mean') for m in METRICS}
    )
    groups = {}
    for region, g in base.groupby('region'):
        if g[unit_col].nunique() >= min_units:
            groups[region] = {m: g[m].to_numpy(float) for m in METRICS}
    regions = sorted(groups)
    vals = {m: np.empty(N_SPLITS, float) for m in METRICS}
    for s in range(N_SPLITS):
        left = {m: [] for m in METRICS}; right = {m: [] for m in METRICS}
        for r in regions:
            n = len(groups[r][next(iter(METRICS))])
            order = rng.permutation(n)
            h = n // 2
            a = order[:h]; b = order[h:2*h]
            for m in METRICS:
                left[m].append(np.nanmean(groups[r][m][a]))
                right[m].append(np.nanmean(groups[r][m][b]))
        for m in METRICS:
            x = np.asarray(left[m]); y = np.asarray(right[m]); ok = np.isfinite(x) & np.isfinite(y)
            vals[m][s] = fast_spearman(x[ok], y[ok]) if ok.sum() >= 4 else np.nan
    rows=[]
    for m, arr in vals.items():
        arr = arr[np.isfinite(arr)]
        rows.append({
            'split_unit': unit_col,
            'metric': m,
            'label': METRICS[m],
            'eligible_regions': len(regions),
            'n_valid_splits': int(len(arr)),
            'median_spearman_rho': float(np.median(arr)),
            'ci_low': float(np.percentile(arr, 2.5)),
            'ci_high': float(np.percentile(arr, 97.5)),
        })
    return rows


def main():
    df = pd.read_csv(INFILE)
    rows = split_half(df, 'session', 4) + split_half(df, 'mouse', 4)
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / 'split_half_reliability.csv', index=False)
    summary = {
        'seed': SEED,
        'n_splits': N_SPLITS,
        'analysis_type': 'post-reveal diagnostic using frozen outcome-blind features only',
        'outcome_values_loaded': False,
        'results': rows,
    }
    (OUTDIR / 'split_half_reliability.json').write_text(json.dumps(summary, indent=2))
    print(out.to_string(index=False))

if __name__ == '__main__':
    main()
