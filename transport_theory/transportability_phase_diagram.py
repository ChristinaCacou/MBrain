"""Deterministic simulation accompanying the sufficient-conditions theorem.

The analytic result is derived in the manuscript. This script illustrates how rank
transport degrades when three invariances are violated: feature/circuit alignment,
intervention-strength invariance, and observational measurement precision.

No empirical manuscript dataset is used here.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 20261012
N_TARGETS = 80
N_REP = 1000
SIGMA_EPS2 = 1.0
MEASUREMENT_NOISE = 0.10
ALIGNMENTS = np.array([0.0, 0.25, 0.50, 0.75, 1.0])
LOG_ALPHA_SD = np.array([0.0, 0.20, 0.40, 0.60, 0.80])
OUT = Path("analysis/transport_theory")
OUT.mkdir(parents=True, exist_ok=True)


def one_rep(rng, alignment, log_alpha_sd, measurement_noise=MEASUREMENT_NOISE):
    z_source = rng.normal(size=N_TARGETS)
    if alignment >= 1.0:
        z_target = z_source.copy()
    elif alignment <= 0.0:
        z_target = rng.normal(size=N_TARGETS)
    else:
        z_target = alignment * z_source + np.sqrt(1.0 - alignment**2) * rng.normal(size=N_TARGETS)

    # q = beta^2 Var(X_i | Z): the unique signal term in the theorem.
    q_source = np.exp(z_source)
    q_target = np.exp(z_target)

    # Gaussian conditional predictive gain: 1/2 log(1 + q/sigma_eps^2).
    cpg = 0.5 * np.log1p(q_source / SIGMA_EPS2)
    if measurement_noise > 0:
        cpg = cpg + rng.normal(scale=measurement_noise * np.std(cpg), size=N_TARGETS)

    # Multiplicative intervention strength; log-mean chosen so E(alpha)=1.
    s = float(log_alpha_sd)
    if s == 0:
        alpha = np.ones(N_TARGETS)
    else:
        alpha = np.exp(rng.normal(loc=-0.5 * s * s, scale=s, size=N_TARGETS))

    causal_impact = alpha**2 * q_target
    return float(spearmanr(cpg, causal_impact).statistic)


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for a in ALIGNMENTS:
        for s in LOG_ALPHA_SD:
            vals = np.array([one_rep(rng, a, s) for _ in range(N_REP)], float)
            rows.append({
                "alignment_r": float(a),
                "log_intervention_strength_sd": float(s),
                "measurement_noise_relative_sd": MEASUREMENT_NOISE,
                "mean_spearman": float(np.mean(vals)),
                "median_spearman": float(np.median(vals)),
                "p025": float(np.quantile(vals, 0.025)),
                "p975": float(np.quantile(vals, 0.975)),
            })
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "phase_diagram.csv", index=False)

    # Exact theorem sanity check: perfect alignment, invariant intervention,
    # and noiseless measurement must preserve the target ranking.
    sanity_rng = np.random.default_rng(SEED + 99)
    sanity = [one_rep(sanity_rng, 1.0, 0.0, measurement_noise=0.0) for _ in range(200)]
    result = {
        "seed": SEED,
        "n_targets_per_rep": N_TARGETS,
        "n_repetitions_per_cell": N_REP,
        "sigma_eps2": SIGMA_EPS2,
        "measurement_noise_relative_sd": MEASUREMENT_NOISE,
        "perfect_invariance_sanity_min_spearman": float(np.min(sanity)),
        "perfect_invariance_sanity_max_spearman": float(np.max(sanity)),
        "interpretation": "Illustrative, not inferential: rank transport is exact under the theorem's invariances and degrades as feature/circuit alignment falls or intervention strength varies across targets.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
