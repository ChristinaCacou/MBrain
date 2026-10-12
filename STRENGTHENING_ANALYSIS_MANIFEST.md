# Strengthening analysis manifest

This manifest records the additional analyses used in the strengthened manuscript. It distinguishes formal results, outcome-frozen numerical tests, post hoc audits, and illustrative simulations.

## 1. Sufficient-conditions theorem and transportability simulation

**Status:** mathematical result plus illustrative simulation; not empirical inference.

The manuscript derives a linear-Gaussian sufficient condition under which conditional predictive gain and an innovation-attenuation perturbation have identical target rankings. `transport_theory/transportability_phase_diagram.py` illustrates violations of feature/circuit alignment and intervention-strength invariance.

- Workflow run: `38322346154`
- Artifact: `transport-theory-phase-diagram`
- Artifact SHA-256: `520f58933c2e13b89b4338849a1bef5a31f66adf4035589996a3e966a7ddba9e`
- Archived result: `analysis/transport_theory/result.json`
- Archived phase table: `analysis/transport_theory/phase_diagram.csv`
- Perfect-invariance/noiseless sanity check: Spearman rho = 1.0 in every check replicate.

## 2. Chen 2021 -> Chen 2024 zero-recalibration transport audit

**Status:** explicitly post hoc. The audit was designed after the Chen 2024 outcomes were known and does not change the inferential status of the frozen 2024 secondary endpoint.

The raw-scale linear relation between bilateral ALM coupling and right-ALM behavioral recovery was fit once in the 18 Chen 2021 mice and applied without slope or intercept refitting to the eight Chen 2024 mice with finite right-side outcomes.

- Workflow run: `38322938472`
- Artifact: `chen-no-refit-transport`
- Artifact SHA-256: `b1b4bfe18043587eae8abc262f6f65d4326468cf1b128a0b627f44d3d7eff84d`
- Source fit: slope 0.2040477017; intercept 0.8360477118
- Target rank association: rho = 0.6666667; exact one-sided 8! permutation P = 0.0415427
- External R2 relative to target-mean variance: 0.4792956
- MAE: 0.0223564; RMSE: 0.0261275
- MAE improvement versus the source-mean baseline: 51.15%
- Calibration slope: 1.23093; calibration intercept: -0.21878
- Archived result: `analysis/chen_no_refit_transport/result.json`
- Archived mouse predictions: `analysis/chen_no_refit_transport/predictions.csv`

## 3. Guo/Inagaki 2017 reciprocal frontal-thalamic loop test

**Status:** frozen numerical source-data test. The source publication and its qualitative conclusions were public, so this is not described as a blinded preregistration. Numerical outcome extraction occurred only after the analysis lock.

### Outcome-blind structural inventories

Before the numerical lock, inventory code read asset paths, subject IDs, recording-region metadata, optogenetic target labels, trial/stimulation timing metadata, and trial/unit counts, but not the trial `response` column or spike times.

- Schema inventory run `38322382356`; artifact SHA-256 `c088508d93cb54445c605ebdf7c9aea40a489b3aa35a5e346c66e40b8b28b1e6`
- Design probe run `38322479892`; artifact SHA-256 `3e3f55b6a329a914ec5f6be11b3473cda92d028c2da78cb2af738fb43ee3b817`
- Full design inventory run `38322585289`; artifact SHA-256 `aea94aa82e3ddd05710c1ab321f3bdd05fc84cb510b5d18a5fb8f7460a0ec34a`
- Optogenetic-target inventory run `38323157040`; artifact SHA-256 `1cf3a4f0391b40c0cf2975ab11913e6b9de10411506b04ec6d225b089e17e1fe`
- Structural design: 39 ALM-recording/left-thalamus-stimulation sessions from 7 mice and 40 thalamus-recording/left-ALM-stimulation sessions from a disjoint 7 mice.

### Numerical lock

`external_replication/GUO2017_CROSSLOOP_LOCK.md` was committed before behavioral responses or spike times were read.

- Lock commit: `561215ed975cff693af2e93795d81dbbeae7aafc`
- Frozen predictor: control-trial late-delay population ROC AUC in the final 0.5 s before the cue.
- Frozen outcome: contralateral behavioral impairment during reciprocal-node delay photoinhibition.
- Frozen primary: pooled within-stratum rank association across the two reciprocal directions, with mouse as the inferential unit.
- Success required >=12 mice total, >=5 per direction, positive association, one-sided P<0.05, and positive leave-one-mouse-out R2.

### Numerical result

The successful numerical workflow used the locked definitions. A preceding workflow attempt failed on an import-only coding error before numerical analysis; the correction changed `roc_auc_score` to its proper scikit-learn import and did not alter any frozen scientific definition.

- Successful workflow run: `38323571805`
- Code commit: `c6da3d6b73603522a4193cf98ec4638667d8e3ec`
- Artifact: `guo2017-crossloop-result`
- Artifact SHA-256: `30866193943a30ccc2ec609becd5ed583f734ad1f3a15257bb0e32d3d92dab07`
- Eligible source data: 78 sessions, 14 mice (7 ALM-recording and 7 thalamus-recording mice)
- Primary pooled within-stratum rank association: -0.1785714
- One-sided within-stratum permutation P: 0.7339365
- Stratified bootstrap 95% interval: -0.7333 to 0.5429
- Leave-one-mouse-out prediction: rho = -0.05934; R2 = -0.04812
- Direction-specific secondary: ALM recording / thalamus perturbation rho = -0.500, exact P = 0.8821; thalamus recording / ALM perturbation rho = 0.1429, exact P = 0.3913
- Prespecified early-window temporal control pooled association: 0.3928571
- Primary success criterion: **not met**
- Archived result: `analysis/guo2017_crossloop/result.json`
- Archived mouse table: `analysis/guo2017_crossloop/mouse_results.csv`

The Guo result is retained as a negative transport stress test. No alternative time window, direction, predictor, outcome, or pooling rule is substituted after reveal.

## Interpretation boundaries

1. The sufficient-conditions theorem is a sufficient, not necessary, result and is specific to the stated linear-Gaussian innovation-attenuation model.
2. The phase diagram is an illustration of theorem violations, not an empirical confidence analysis.
3. The Chen zero-recalibration audit is post hoc and cannot upgrade the frozen 2024 right-side secondary to a confirmatory replication.
4. The Guo test is numerically frozen before outcome extraction but is not a blinded preregistration because the source publication was already public.
5. Negative frozen results are retained without endpoint or metric substitution.
