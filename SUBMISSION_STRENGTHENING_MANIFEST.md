# Submission strengthening manifest

Frozen 2026-10-12 after completion of the analyses listed below. This file distinguishes prospectively locked numerical tests from post hoc diagnostics.

## Sufficient-conditions theorem and phase diagram
- Script: `theory/transportability_sufficient_conditions.py`
- Result specification: `analysis/transportability_theory/result.json`
- Seed: 20261012
- Status: illustrative theoretical benchmark, not fitted to biological data and not a causal-identification result.
- Core result: in the linear-Gaussian benchmark, CPG = 0.5 log(1 + beta^2 v / sigma^2). Under an ideal residual-feature clamp, squared conditional-mean displacement is D = beta^2 v. With invariant sigma^2, CPG ranks D exactly; heterogeneity in residual-noise scale or intervention/alignment scale can break the ordering.

## Chen 2021 to Chen 2024 zero-refit audit
- Training rows: `analysis_inputs/chen_mouse_rows.csv` (18 Chen 2021 mice).
- Test rows: `analysis_inputs/chen2024_right_behavior_mouse_rows.csv` (8 Chen 2024 right-behavior secondary mice).
- Script: `external_replication/chen_zero_refit_prediction.py`
- Archived result: `analysis/chen_zero_refit_prediction/result.json`
- 2021 model: recovery = 0.83604771 + 0.20404770 * coupling.
- Applied to Chen 2024 without refitting or recalibration: external R2 = 0.47929563; MAE = 0.02235645; Spearman rho = 0.66666667; 200,000-permutation P for external R2 = 0.01878991; bootstrap 95% R2 = [-2.50511943, 0.74636826].
- Status: post hoc diagnostic. The Chen 2024 right-side outcome was already known when this audit was conceived, so it does not convert the prespecified secondary result into confirmatory replication.

## Frozen Guo/Inagaki matched-system test
- Dataset: DANDI 000009 version 0.220126.1903.
- Lock: `external_replication/GUO_INAGAKI_2017_LOCK.md`, committed before numerical perturbation-outcome extraction.
- Analysis: `external_replication/guo_inagaki_analysis.py`.
- Workflow run: 38327399003.
- Artifact: `guo-inagaki-locked-results`, SHA-256 `f770f3bcabb1fc7af6a2c19749797ddf57ce888485b731ba14fc463607b478d5`.
- 116 ecephys+ogen NWB assets inspected; 18 sessions from 5 mice met all frozen eligibility/matched-site criteria.
- Session-level primary: rho = 0.46130031; one-sided permutation P = 0.02800794; bootstrap 95% [-0.00104058, 0.78458857]; leave-one-session-out R2 = -0.06788105.
- Mouse-mean sensitivity: n = 5, rho = -0.10; P = 0.60857478; leave-one-mouse-out R2 = -2.4206.
- Frozen success criterion required rho>0, P<0.05, positive held-out R2, >=10 eligible sessions, and >=6 mice. Criterion was not met because only five mice were eligible and held-out R2 was negative.
- ALM and THAL strata are sensitivity analyses only and are not promoted as anatomical specificity.

## Interpretation hierarchy
1. The theorem supplies sufficient conditions under which observational predictive gain and idealized perturbation impact are rank-equivalent; it does not identify causality from observation alone.
2. The zero-refit Chen audit gives encouraging cross-study predictive evidence but is explicitly post hoc.
3. The Guo/Inagaki analysis is a genuinely frozen numerical external test and is reported as a failed/ambiguous test under its prespecified success criterion despite a nominally positive session-level association.

No endpoint, eligibility rule, predictor, or success threshold was changed after the Guo/Inagaki numerical outcome extraction to improve its result.
