import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

SEED=20261012
OUT=Path('analysis/transportability_theory'); OUT.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(SEED)
sigma_grid=np.linspace(0,1.5,21)
align_grid=np.linspace(0,1.5,21)
R=np.zeros((len(sigma_grid),len(align_grid)))
N_TARGETS=400; REPS=30

for i,sig_sd in enumerate(sigma_grid):
    for j,align_sd in enumerate(align_grid):
        vals=[]
        for _ in range(REPS):
            # D = beta^2 Var(X|Z), the idealized squared conditional-mean
            # displacement induced by clamping the residual feature component.
            D=np.exp(rng.normal(0,1.0,N_TARGETS))
            sigma2=np.exp(rng.normal(0,sig_sd,N_TARGETS))
            cpg=0.5*np.log1p(D/sigma2)
            # Target-specific effective perturbation/alignment scale.
            a=np.exp(rng.normal(0,align_sd,N_TARGETS))
            observed_perturbation_effect=D*a
            vals.append(float(spearmanr(cpg,observed_perturbation_effect).statistic))
        R[i,j]=float(np.mean(vals))
np.savetxt(OUT/'rank_agreement.csv',R,delimiter=',')
result={
 'status':'ILLUSTRATIVE_SUFFICIENT_CONDITIONS_SIMULATION',
 'seed':SEED,'n_targets_per_cell':N_TARGETS,'replicates_per_cell':REPS,
 'residual_noise_log_sd_grid':sigma_grid.tolist(),
 'intervention_alignment_log_sd_grid':align_grid.tolist(),
 'theorem':'For Y=mu(Z)+beta X+epsilon with epsilon~N(0,sigma^2) and Var(X|Z)=v, conditional predictive gain is 0.5*log(1+beta^2*v/sigma^2). Under an ideal intervention that clamps the residual component of X to E[X|Z], the expected squared displacement of the conditional mean is D=beta^2*v. If sigma^2 is invariant across targets, CPG is strictly increasing in D and therefore preserves causal-impact rank. Heterogeneity in residual noise or effective perturbation/alignment scale can break this ordering.',
 'interpretation':'Illustrative theoretical benchmark, not a fitted biological model and not a causal-identification result.'
}
(OUT/'result.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
