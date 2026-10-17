"""
Supplementary Code 3
Retrospective public-data reanalysis of Stentor coeruleus learning/forgetting RNA-seq.

Source study:
Rajan DH et al. Molecular pathways for learning in the single-cell Stentor coeruleus.
Current Biology 36, 2367-2381.e9 (2026). DOI: 10.1016/j.cub.2026.03.080

Public data:
Dryad DOI: 10.5061/dryad.z34tmpgrm
Processed analysis outputs: https://github.com/aralbright/2025_stentor_learning

This analysis is retrospective and is NOT part of the prospectively locked Stentor
perturbation experiment described in the associated manuscript.

Python 3.11+
Dependencies: numpy, pandas, scipy, scikit-learn
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest
from sklearn.decomposition import PCA

SEED = 20261017
N_BOOT = 3000
N_PERM = 1000
OUT = Path("retrospective_stentor_omics")
OUT.mkdir(exist_ok=True)

BASE = "https://raw.githubusercontent.com/aralbright/2025_stentor_learning/main/exp2/"
OBS = BASE + "exp2_rna_no0_sleuth_obsnorm.csv"
TRAIN = BASE + "exp2_rna_t19ut19_sleuth_wald_v2.csv"
FORGET = BASE + "exp2_rna_t19f90_sleuth_wald.csv"

EXPECTED = {
    "t19": ["DR4876_01", "DR4876_04", "DR4876_07"],
    "u19": ["DR4876_03", "DR4876_06", "DR4876_09"],
    "f10": ["DR4876_10", "DR4876_12", "DR4876_14"],
    "f90": ["DR4876_16", "DR4876_18", "DR4876_20"],
}

def map_sample_columns(columns):
    columns = [str(c) for c in columns]
    def find(sample):
        if sample in columns:
            return sample
        suffix = sample[-2:]
        candidates = [c for c in columns if sample in c or c.endswith("_" + suffix) or c.endswith(suffix)]
        if len(candidates) != 1:
            raise RuntimeError(f"Could not uniquely map {sample}; columns={columns}")
        return candidates[0]
    return {g: [find(s) for s in samples] for g, samples in EXPECTED.items()}

def load_matrix():
    raw = pd.read_csv(OBS, index_col=0).apply(pd.to_numeric, errors="coerce")
    raw = raw.loc[:, raw.notna().any(axis=0)]
    groups = map_sample_columns(raw.columns)
    used = sum(groups.values(), [])
    raw = raw[used].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    keep = (raw.mean(axis=1) >= 1.0) & (raw.var(axis=1) > 0)
    x = np.log1p(raw.loc[keep])
    return raw, x, groups, used

def learning_axis_analysis(x, groups, rng):
    mu_u = x[groups["u19"]].mean(axis=1).to_numpy()
    mu_t = x[groups["t19"]].mean(axis=1).to_numpy()
    delta = mu_t - mu_u
    denom = float(delta @ delta)
    def group_index(g):
        mu = x[groups[g]].mean(axis=1).to_numpy()
        return float((mu - mu_u) @ delta / denom)
    full = {g: group_index(g) for g in groups}
    boot = []
    for _ in range(N_BOOT):
        idx = {g: rng.integers(0, 3, 3) for g in groups}
        ub = x[[groups["u19"][i] for i in idx["u19"]]].mean(axis=1).to_numpy()
        tb = x[[groups["t19"][i] for i in idx["t19"]]].mean(axis=1).to_numpy()
        d = tb - ub
        q = float(d @ d)
        if q <= 1e-15:
            continue
        vals = {}
        for g in ("f10", "f90"):
            m = x[[groups[g][i] for i in idx[g]]].mean(axis=1).to_numpy()
            vals[g] = float((m - ub) @ d / q)
        boot.append((vals["f10"], vals["f90"], vals["f10"] - vals["f90"]))
    b = np.asarray(boot)
    ci = {
        "f10": np.quantile(b[:, 0], [0.025, 0.975]).tolist(),
        "f90": np.quantile(b[:, 1], [0.025, 0.975]).tolist(),
        "f10_minus_f90": np.quantile(b[:, 2], [0.025, 0.975]).tolist(),
    }
    cross = []
    for r in range(3):
        train = [j for j in range(3) if j != r]
        ub = x[[groups["u19"][j] for j in train]].mean(axis=1).to_numpy()
        tb = x[[groups["t19"][j] for j in train]].mean(axis=1).to_numpy()
        d = tb - ub
        q = float(d @ d)
        for g in groups:
            score = float((x[groups[g][r]].to_numpy() - ub) @ d / q)
            cross.append((r + 1, g, score))
    cross_df = pd.DataFrame(cross, columns=["replicate", "condition", "memory_index"])
    cross_df.to_csv(OUT / "crossfitted_memory_axis.csv", index=False)
    return full, ci, cross_df

def pca_for_visualization(x, groups, used):
    nvar = min(1000, x.shape[0])
    variable = x.var(axis=1).nlargest(nvar).index
    a = x.loc[variable, used].T
    z = (a - a.mean(axis=0)) / a.std(axis=0, ddof=1).replace(0, 1)
    pca = PCA(n_components=2, random_state=SEED)
    pcs = pca.fit_transform(z)
    rev = {sample: group for group, samples in groups.items() for sample in samples}
    df = pd.DataFrame({"sample": used, "condition": [rev[s] for s in used], "PC1": pcs[:, 0], "PC2": pcs[:, 1]})
    df.to_csv(OUT / "pca_coordinates.csv", index=False)
    return pca.explained_variance_ratio_.tolist()

def tidy_wald(url):
    d = pd.read_csv(url)
    return d[["target_id", "b", "qval"]].rename(columns={"target_id": "gene"}).dropna()

def differential_alignment(rng):
    train = tidy_wald(TRAIN)
    forget = tidy_wald(FORGET)
    j = train.merge(forget, on="gene", suffixes=("_train", "_forget"))
    x = j["b_train"].to_numpy(float)
    y = j["b_forget"].to_numpy(float)
    cosine = float(x @ y / (np.linalg.norm(x) * np.linalg.norm(y)))
    rho, rho_p = spearmanr(x, y)
    sig = j[j["qval_train"] < 0.05].copy()
    n = len(sig)
    same = int((np.sign(sig["b_train"]) == np.sign(sig["b_forget"])).sum())
    sign_p = float(binomtest(same, n, 0.5, alternative="greater").pvalue)
    ratio = float(np.median(np.abs(sig["b_forget"].to_numpy(float)) / (np.abs(sig["b_train"].to_numpy(float)) + 1e-12)))
    w2 = sig["b_train"].to_numpy(float) ** 2
    n_eff = float(w2.sum() ** 2 / (w2 @ w2))
    ordered = np.sort(w2)[::-1]
    cumulative = np.cumsum(ordered) / ordered.sum()
    n50 = int(np.searchsorted(cumulative, 0.50) + 1)
    n80 = int(np.searchsorted(cumulative, 0.80) + 1)
    nx, ny = np.linalg.norm(x), np.linalg.norm(y)
    exceed = 0
    for _ in range(N_PERM):
        yp = rng.permutation(y)
        exceed += (x @ yp / (nx * ny) >= cosine)
    cosine_perm_p = float((exceed + 1) / (N_PERM + 1))
    j.to_csv(OUT / "gene_level_training_forgetting_alignment.csv", index=False)
    sig.sort_values("qval_train").to_csv(OUT / "significant_training_genes.csv", index=False)
    effect = sig.copy()
    effect["effect_energy"] = effect["b_train"].astype(float) ** 2
    effect = effect.sort_values("effect_energy", ascending=False)
    effect["effect_energy_fraction"] = effect["effect_energy"] / effect["effect_energy"].sum()
    effect["cumulative_effect_energy"] = effect["effect_energy_fraction"].cumsum()
    effect[["gene", "effect_energy_fraction", "cumulative_effect_energy"]].to_csv(OUT / "training_effect_energy_curve.csv", index=False)
    return {
        "common_genes": int(len(j)), "training_q_lt_0.05": int(n), "cosine": cosine,
        "cosine_gene_identity_permutation_p": cosine_perm_p, "spearman_rho": float(rho),
        "spearman_p": float(rho_p), "same_direction_n": same,
        "same_direction_fraction": float(same / n), "same_direction_binomial_p": sign_p,
        "median_abs_t19_minus_f90_over_t19_minus_untrained": ratio,
        "effective_number_of_training_effect_genes": n_eff,
        "genes_for_50pct_training_effect_energy": n50, "genes_for_80pct_training_effect_energy": n80,
    }

def main():
    rng = np.random.default_rng(SEED)
    raw, x, groups, used = load_matrix()
    memory, bootstrap_ci, cross = learning_axis_analysis(x, groups, rng)
    pca_var = pca_for_visualization(x, groups, used)
    alignment = differential_alignment(rng)
    results = {
        "source": {"study": "Rajan et al., Current Biology (2026)", "dryad_doi": "10.5061/dryad.z34tmpgrm", "processed_data_repository": "aralbright/2025_stentor_learning"},
        "status": "retrospective public-data reanalysis; not part of prospective preregistration",
        "original_matrix_shape": [int(raw.shape[0]), int(raw.shape[1])],
        "genes_after_expression_filter": int(x.shape[0]), "mapped_columns": groups,
        "full_learning_axis_memory_index": memory, "bootstrap_95pct_CI": bootstrap_ci,
        "crossfitted_condition_means": cross.groupby("condition")["memory_index"].mean().to_dict(),
        "pca_explained_variance_ratio": pca_var, "differential_alignment": alignment,
    }
    (OUT / "results.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
