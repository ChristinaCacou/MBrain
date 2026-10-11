import ast
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LinearRegression

SEED = 20261011
N_PERM = 200000
N_BOOT = 50000
ROOT = Path('analysis/pui_blind_fixedC')
FEATURE = ROOT / 'blind_features_region.csv'
LOCK = ROOT / 'blind_lock.json'
OUT = ROOT / 'reveal'
OUT.mkdir(parents=True, exist_ok=True)
OUTCOME_COMMIT = 'a17677ef2ab38c198761ad66f8a8548e4b857901'
OUTCOME_URL = f'https://raw.githubusercontent.com/elliottower/neural-geometry-reliability/{OUTCOME_COMMIT}/experiments/shared_bundle.py'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load_outcomes():
    text = urlopen(OUTCOME_URL).read().decode('utf-8')
    tree = ast.parse(text)
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == 'SILENCING_EFFECTS':
                    obj = ast.literal_eval(node.value)
                    return {str(k): float(v) for k, v in obj.items()}, hashlib.sha256(text.encode()).hexdigest()
    raise RuntimeError('SILENCING_EFFECTS not found')


def spearman(x, y):
    r, p = spearmanr(np.asarray(x, float), np.asarray(y, float))
    return float(r), float(p)


def perm_p(x, y, observed, rng, alternative='greater'):
    x = np.asarray(x, float); y = np.asarray(y, float)
    count = 0
    for _ in range(N_PERM):
        rp, _ = spearmanr(x, rng.permutation(y))
        if alternative == 'greater':
            count += (rp >= observed - 1e-15)
        else:
            count += (abs(rp) >= abs(observed) - 1e-15)
    return float((count + 1) / (N_PERM + 1))


def boot_ci(x, y, rng):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x); vals=[]
    for _ in range(N_BOOT):
        ii = rng.integers(0, n, size=n)
        if len(np.unique(x[ii])) < 2 or len(np.unique(y[ii])) < 2:
            continue
        r, _ = spearmanr(x[ii], y[ii])
        if np.isfinite(r): vals.append(r)
    if not vals: return [np.nan, np.nan]
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def partial_rank(x, y, controls):
    xr = rankdata(np.asarray(x, float)); yr = rankdata(np.asarray(y, float)); Z = np.column_stack([rankdata(np.asarray(z, float)) for z in controls])
    X = np.column_stack([np.ones(len(xr)), Z])
    bx = np.linalg.lstsq(X, xr, rcond=None)[0]; by = np.linalg.lstsq(X, yr, rcond=None)[0]
    rx = xr - X @ bx; ry = yr - X @ by
    r, p = spearmanr(rx, ry)
    return float(r), float(p)


def loo_linear_predictions(x, y):
    x=np.asarray(x,float).reshape(-1,1); y=np.asarray(y,float); pred=np.zeros(len(y))
    for i in range(len(y)):
        keep=np.arange(len(y))!=i
        m=LinearRegression().fit(x[keep],y[keep]); pred[i]=m.predict(x[i:i+1])[0]
    r,p=spearman(pred,y); mse=float(np.mean((pred-y)**2)); base=float(np.mean((y-y.mean())**2)); r2=1-mse/base if base>0 else np.nan
    return pred, {'loo_spearman_rho':r,'loo_spearman_p':p,'loo_r2':float(r2),'loo_mse':mse}


def main():
    lock=json.loads(LOCK.read_text())
    expected=lock['feature_sha256']['blind_features_region.csv']; actual=sha256(FEATURE)
    if actual != expected:
        raise RuntimeError(f'BLIND LOCK VIOLATION: expected {expected}, observed {actual}')
    if lock.get('perturbation_labels_loaded') is not False:
        raise RuntimeError('phase1 lock does not certify outcome blindness')
    features=pd.read_csv(FEATURE)
    outcomes, outcome_source_sha=load_outcomes()
    outcome_df=pd.DataFrame({'region':list(outcomes.keys()),'silencing_effect':list(outcomes.values())})
    merged=outcome_df.merge(features,on='region',how='inner').sort_values('region').reset_index(drop=True)
    if len(merged)!=len(outcomes):
        missing=sorted(set(outcomes)-set(merged.region)); raise RuntimeError(f'Missing outcome regions: {missing}')
    y=merged.silencing_effect.to_numpy(float)
    metrics={
        'pui_mean':'PUI (primary)',
        'pui_median':'PUI median (sensitivity)',
        'marginal_mean':'marginal conditional-on-history information',
        'local_bacc_mean':'local choice decoder balanced accuracy',
        'local_auc_mean':'local choice decoder AUC',
        'snr_mean':'condition SNR',
        'mean_window_spikes':'mean early-window spikes',
        'neurons_recorded_mean':'recorded-neuron yield',
        'n_sessions':'number of recording sessions'
    }
    rng=np.random.default_rng(SEED); results=[]
    for key,label in metrics.items():
        x=merged[key].to_numpy(float); ok=np.isfinite(x)&np.isfinite(y); xx=x[ok]; yy=y[ok]
        rho,p=spearman(xx,yy); one=perm_p(xx,yy,rho,rng,'greater'); two=perm_p(xx,yy,rho,rng,'two'); ci=boot_ci(xx,yy,rng); pred,loo=loo_linear_predictions(xx,yy)
        results.append({'metric':key,'label':label,'n':int(len(xx)),'spearman_rho':rho,'asymptotic_p_two_sided':p,'permutation_p_one_sided_positive':one,'permutation_p_two_sided':two,'bootstrap_ci_low':ci[0],'bootstrap_ci_high':ci[1],**loo})
    res=pd.DataFrame(results)
    # Primary partial-rank sensitivity controlling sampling and signal quality.
    prho,pp=partial_rank(merged.pui_mean,y,[merged.neurons_recorded_mean,merged.snr_mean])
    # Direct predeclared comparison: does PUI beat marginal information?
    pui_rho=float(res.loc[res.metric=='pui_mean','spearman_rho'].iloc[0]); marg_rho=float(res.loc[res.metric=='marginal_mean','spearman_rho'].iloc[0]); delta=pui_rho-marg_rho
    xp=merged.pui_mean.to_numpy(float); xm=merged.marginal_mean.to_numpy(float)
    count=0; rng_delta=np.random.default_rng(SEED+77)
    for _ in range(N_PERM):
        yp=rng_delta.permutation(y); dp=spearmanr(xp,yp).statistic-spearmanr(xm,yp).statistic
        count += (dp >= delta - 1e-15)
    delta_p=float((count+1)/(N_PERM+1))
    # Leave-one-region-out rho stability of primary association.
    loo_rhos=[]
    for i in range(len(merged)):
        keep=np.arange(len(merged))!=i; loo_rhos.append(float(spearmanr(xp[keep],y[keep]).statistic))
    primary={
        'feature_sha256_verified':actual,
        'blind_lock_script_sha256':lock['script_sha256'],
        'outcome_source_url':OUTCOME_URL,
        'outcome_source_sha256':outcome_source_sha,
        'n_regions':int(len(merged)),
        'regions':merged.region.tolist(),
        'primary_metric':'pui_mean',
        'primary_spearman_rho':pui_rho,
        'primary_permutation_p_one_sided_positive':float(res.loc[res.metric=='pui_mean','permutation_p_one_sided_positive'].iloc[0]),
        'primary_permutation_p_two_sided':float(res.loc[res.metric=='pui_mean','permutation_p_two_sided'].iloc[0]),
        'primary_bootstrap_95_ci':[float(res.loc[res.metric=='pui_mean','bootstrap_ci_low'].iloc[0]),float(res.loc[res.metric=='pui_mean','bootstrap_ci_high'].iloc[0])],
        'partial_rank_rho_controlling_neuron_yield_and_snr':prho,
        'partial_rank_asymptotic_p':pp,
        'marginal_rho':marg_rho,
        'delta_rho_pui_minus_marginal':delta,
        'delta_rho_permutation_p_one_sided':delta_p,
        'leave_one_region_out_rho_min':float(np.min(loo_rhos)),
        'leave_one_region_out_rho_median':float(np.median(loo_rhos)),
        'leave_one_region_out_rho_max':float(np.max(loo_rhos)),
        'success_criterion':'positive primary association with one-sided permutation p<0.05 AND PUI rho > marginal rho',
    }
    primary['success']=bool(primary['primary_permutation_p_one_sided_positive']<0.05 and delta>0)
    merged.to_csv(OUT/'revealed_region_table.csv',index=False); res.to_csv(OUT/'metric_comparison.csv',index=False); (OUT/'primary_result.json').write_text(json.dumps(primary,indent=2,sort_keys=True))
    lines=['# Exact blinded causal-prediction test','',f"**Result:** {'PASS' if primary['success'] else 'FAIL'}",'',f"Frozen feature SHA-256: `{actual}`",f"Outcome source: `{OUTCOME_COMMIT}` / `SILENCING_EFFECTS`",'',f"Primary PUI vs silencing: Spearman rho = {pui_rho:.4f}; one-sided permutation p = {primary['primary_permutation_p_one_sided_positive']:.6g}; two-sided permutation p = {primary['primary_permutation_p_two_sided']:.6g}.",f"Bootstrap 95% CI: [{primary['primary_bootstrap_95_ci'][0]:.4f}, {primary['primary_bootstrap_95_ci'][1]:.4f}].",f"Marginal-information baseline rho = {marg_rho:.4f}; delta rho = {delta:.4f}; delta permutation p = {delta_p:.6g}.",f"Partial-rank rho controlling neuron yield + SNR = {prho:.4f}.",'','## All baselines','',res.to_markdown(index=False)]
    (OUT/'RESULT.md').write_text('\n'.join(lines))
    print('\n'.join(lines[:12])); print(res.to_string(index=False)); print(json.dumps(primary,indent=2))
if __name__=='__main__': main()
