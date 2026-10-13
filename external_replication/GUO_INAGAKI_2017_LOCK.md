# Frozen numerical specification: Guo/Inagaki frontal-thalamic recording-perturbation test

Frozen 2026-10-12 before numerical perturbation-outcome extraction from DANDI 000009.

## Scientific question
Does control-trial delay-period population choice coding predict behavioral impairment produced by photoinhibition of the same recorded anatomical structure in an independent frontal-thalamic circuit?

## Dataset
DANDI 000009, published version 0.220126.1903 (Guo et al., Nature 2017). Public qualitative conclusions are known from the paper; this lock concerns numerical extraction and the predictor-outcome association.

## Inferential unit
Recording session. Mouse-clustered sensitivity will be reported when multiple eligible sessions occur in the same subject.

## Eligibility
Use extracellular NWB sessions only. A session is eligible when: (1) recorded units are in ALM or thalamus; (2) the session contains good, correct, non-photostimulation control trials for both instructed left and right responses; (3) the session contains delay-period photostimulation trials with valid behavioral responses; (4) the photostimulation anatomical structure matches the recorded structure after the prespecified aliases ALM->ALM and VM/VAL/thalamus->THAL; (5) at least 5 quality-passed recorded units and at least 8 control trials per choice are available; and (6) at least 8 matched-site photostimulation trials and 8 comparable no-stimulation trials are available for the behavioral outcome.

## Predictor
For each eligible session, using control trials only, estimate delay-period population choice coding. The delay window is the final 600 ms before the response/go event if explicit task-event timestamps permit this; if the NWB exposes an author-defined delay epoch instead, use that delay epoch without tuning. For every quality-passed unit calculate trial spike count in the frozen window. Standardize unit counts using control trials. Use deterministic 5-fold stratified cross-validation with logistic regression (L2 penalty, C=1, random seed 20261012) to predict left versus right instructed response. The session predictor is cross-validated balanced accuracy minus 0.5. No perturbation trials enter feature construction or model fitting.

A prespecified diagnostic predictor is the mean absolute single-unit left/right standardized count difference in the same window. It cannot replace the primary predictor.

## Outcome
The behavioral perturbation outcome is matched-site delay photoinhibition impairment: absolute reduction in correct-response probability under matched-site delay photostimulation relative to comparable no-stimulation trials, calculated within session. If trial structure distinguishes instructed side, compute side-specific impairment first and average across instructed sides with finite estimates. No sign flipping or endpoint selection after extraction.

## Primary test
Across eligible sessions, Spearman correlation between frozen control-trial population-code strength and matched-site behavioral impairment. Directional alternative rho > 0. Exact/permutation P value from 500,000 permutations of outcome labels across sessions, random seed 20261012. Bootstrap 95% interval from 50,000 session resamples. Also report leave-one-session-out linear prediction R2 against a training-mean baseline.

## Cluster sensitivity
If any mouse contributes more than one eligible session, repeat the association on mouse-mean predictor/outcome values and perform mouse-label permutation. This is the preferred robustness analysis but does not replace the frozen session-level primary.

## Stratified sensitivity
Report ALM and THAL associations separately when each stratum has at least 5 eligible sessions. Do not claim anatomical specificity from one significant and one nonsignificant stratum without a formal interaction/permutation contrast.

## Success criterion
The external matched-system test is counted as successful only if the primary association has rho > 0, one-sided permutation P < 0.05, leave-one-session-out R2 > 0, and at least 10 eligible sessions from at least 6 mice. If eligibility falls below these minima, label the test structurally underpowered/ineligible rather than successful.

## Multiplicity and stopping
Only the population decoder association is confirmatory for this frozen numerical test. The diagnostic single-unit predictor, anatomical strata, alternate reasonable trial windows required by schema constraints, and cluster analyses are sensitivity analyses. No predictor, endpoint, anatomical matching rule, window, or success criterion may be changed after numerical outcome extraction to improve the result.
