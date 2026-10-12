# Submission analysis manifest — 2026-10-12

This manifest freezes the computational record used by the revised manuscript. It distinguishes locked/prespecified tests from post hoc robustness analyses and records the workflow artifacts needed to audit the reported numbers.

## Outcome-blind Steinmetz -> Zatka-Haas CPG test
- Branch: `causal-pui-blinded-fixedC-20261011`
- Phase-1 feature SHA-256: `93813e3e0232d0dd7aac93d7322fe0aef94f7ce0e33ca32ff87f188020c6d236`
- Phase-1 script SHA-256: `def269ec1743df0f5cbc78f1f6b9847ac6f94c6dd157a132dad7354f5358816f`
- Pinned secondary outcome source commit: `a17677ef2ab38c198761ad66f8a8548e4b857901`
- Primary locked result: n=12 regions; Spearman rho=0.384615; one-sided permutation P=0.109; LORO R2=-0.285.
- Inferential status: outcome-blind computational lock, not formal preregistration.

## Post-reveal Zatka-Haas primary-data provenance audit
- Branch: `zatka-primary-spatial-chen-contrast-20261012`
- Workflow run: `38313250735`
- Artifact: `zatka-primary-reconstruction`
- Artifact SHA-256: `551ce9ba0aa36c370e4f42acc8fe26ac95de8bbb1407df24002a43eb7c59bad2`
- Primary trial table: 47,002 trials, 91 sessions, 52 laser coordinates; 6,045 equal nonzero bilateral-contrast trials used to reproduce the source Fig. 3 coordinate-level choice-change calculation.
- Anatomy-only Beryl mapping recovered 7/12 frozen regional targets from 48/52 coordinates. On those seven, CPG versus reconstructed raw effect: rho=-0.643, one-sided P=0.945. The independently reconstructed regional values did not reproduce the pinned secondary 12-region table (Spearman=0.000; MAE=0.1860).
- Inferential status: post-reveal provenance sensitivity. It does not replace the locked secondary outcome table and is not an independent confirmatory test.

## Post-reveal spatial sensitivity of the 12-region locked test
- Workflow run: `38312782135`
- Artifact: `zatka-spatial-sensitivity-final`
- Artifact SHA-256: `47113d22d72ab6ee6b732e6b4353429e1d6518cbd990d68ccd2ef2ffe7000710`
- Coordinate-adjusted partial-rank rho=0.358; one-sided permutation P=0.171.
- Two-Moran-mode-adjusted rho=0.399; one-sided permutation P=0.126.
- Exponential spatial-GLS rank slopes remained nonsignificant across length scales (one-sided P=0.267, 0.396, 0.409).
- Inferential status: post-reveal descriptive robustness.

## Ghaderi matched-system analysis and source-data reconstruction
- Branch: `ghaderi-matched-validation-20261011`
- Prespecified processed-data result: pooled within-epoch rank association=0.733; exact one-sided P=1/120=0.0083.
- Raw/source-data matrix workflow: `38282285751`.
- Forced aggregate workflow: `38285040160`.
- Artifact: `ghaderi-raw-final-forced`; SHA-256 `dac8d3939c6784da45e5162a6c33ea8581e169faf22126bd304a3d90c4ab0ce1`.
- Post hoc source-data reconstruction: 35 ephys NWB files, 267 session-region-epoch decoder rows, 226 optogenetic rows, all 15 region-epoch cells; association=0.667, exact one-sided P=0.0167.
- Inferential status: processed-data numerical specification prespecified before automated extraction but not outcome-blind to the publication; raw reconstruction is post hoc technical robustness.

## Raw Chen bilateral-ALM analysis
- Branch: `raw-chen-alm3-replication-20261012`
- Session matrix workflow: `38270386016`.
- Forced aggregate workflow: `38278814249`.
- Artifact: `chen-direct-final-forced`; SHA-256 `8e7a67bba3c5766ee911cccd620f56e56a88c770a37d101a01203b2eefa8af4a`.
- Right-behavior endpoint: n=18 mice; rho=0.552; one-sided P=0.00977; bootstrap 95% CI 0.105–0.819; LOMO R2=0.046; Bonferroni/BH adjusted P/q=0.0391.
- Other direct endpoints did not meet the positive threshold.
- Inferential status: de novo raw-session analysis but exploratory because source publication and an earlier processed bridge were known.

## Chen paired specificity audit
- Workflow run: `38312829941`.
- Artifact: `chen-specificity-final`; SHA-256 `ca956da215c3e81876d7cb6899bd242ec273016ed5f982256cd6f3f21a114a44`.
- Right behavior minus left behavior: delta rho=0.576; exact one-sided P=0.0751; bootstrap 95% CI -0.027 to 1.163.
- Right behavior minus right neuronal: delta rho=0.223; P=0.245; CI -0.379 to 0.836.
- Right neuronal minus left neuronal: delta rho=0.250; P=0.257; CI -0.401 to 0.908.
- Inferential status: post hoc. These tests do not establish statistically significant side or endpoint specificity.

## Frozen CRCNS alm-3 transport test
- Branch: `raw-chen-alm3-replication-20261012`
- Common lock: `raw_replication/analysis_lock.json`, frozen before numerical alm-3 outcome extraction.
- Definitive corrected aggregate workflow: `38273421403`.
- Artifact: `alm3-definitive-corrected-final`.
- Primary: n=10 mice; local control-only AUC versus ipsilateral neuronal recovery rho=-0.455; one-sided P=0.911; bootstrap 95% CI -0.862 to 0.220; LOMO R2=-0.636.
- Inferential status: frozen numerical external transport test, not an outcome-blind preregistration and not a test of Chen bilateral coupling.

## Submission interpretation lock
The manuscript may claim that observational predictors can acquire intervention-predictive value under specific feature/circuit/intervention alignments and that such predictiveness does not automatically transport. It must not claim that the Chen side or endpoint specificity is statistically established, that the Zatka-Haas primary regional mapping was independently reproduced, or that alm-3 replicated the Chen bilateral-coupling feature.

The parent commit immediately before this manifest was `405b3cfc371d0673026c6e44107287a2350a3124`. The commit containing this manifest should be treated as the submission analysis snapshot for this branch.
