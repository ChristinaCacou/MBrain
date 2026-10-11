import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 20261011
N_SPLITS = 20000
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


def split_half(df, unit_col, min_units=4):
    rng = np.random.default_rng(SEED + (0 if unit_col == 'session' else 1000))
    base = df.groupby(['region', unit_col], as_index=False).agg(
        **{m: (m, 'mean') for m in METRICS}
    )
    counts = base.groupby('region')[unit_col].nunique()
    eligible = sorted(counts[counts >= min_units].index)
    base = base[base.region.isin(eligible)].copy()
    values = {m: [] for m in METRICS}
    n_used = []
    for _ in range(N_SPLITS):
        left = {}
        right = {}
        for region, g in base.groupby('region'):
            order = rng.permutation(len(g))
            h = len(g) // 2
            if h < 1:
                continue
            a = g.iloc[order[:h]]
            b = g.iloc[order[h:2*h]]
            left[region] = a
            right[region] = b
        common = sorted(set(left) & set(right))
        n_used.append(len(common))
        for m in METRICS:
            x = np.array([left[r][m].mean() for r in common], float)
            y = np.array([right[r][m].mean() for r in common], float)
            ok = np.isfinite(x) & np.isfinite(y)
            rho = spearmanr(x[ok], y[ok]).statistic if ok.sum() >= 4 else np.nan
            values[m].append(rho)
    rows = []
    for m, arr in values.items():
        arr = np.asarray(arr, float)
        arr = arr[np.isfinite(arr)]
        rows.append({
            'split_unit': unit_col,
            'metric': m,
            'label': METRICS[m],
            'eligible_regions': len(eligible),
            'median_regions_per_split': float(np.median(n_used)),
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
