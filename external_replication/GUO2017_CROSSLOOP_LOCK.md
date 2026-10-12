# Frozen numerical analysis: Guo/Inagaki 2017 cross-loop transport test

**Lock date:** 2026-10-12  
**Dataset:** DANDI 000009, version `0.220126.1903`  
**Qualitative source result is public:** this is a frozen *numerical* source-data test, not a blinded preregistration.

## Outcome-separation statement

Before this lock was committed, structural inventory code read asset paths, subject IDs, recording-region metadata, optogenetic-site labels, trial-type labels, stimulation timing/coordinates, and trial/unit counts. It did **not** read the trial `response` column or unit spike times. The structural inventory established two cross-loop strata with disjoint mouse sets: ALM recording during left-thalamus delay photoinhibition (39 sessions, 7 mice) and thalamus recording during left-ALM delay photoinhibition (40 sessions, 7 mice). No numerical behavioral impairment or neural coding statistic has been extracted before this lock.

## Scientific question

Does the strength of a control-trial late-delay directional population code in one node of the frontal-thalamic loop predict the behavioral impact of delay-period photoinhibition of the reciprocal node?

## Frozen strata

Only extracellular NWB sessions satisfying one of these two exact source-exported metadata patterns are eligible:

1. electrode `location` contains `brain_region: ALM` and the NWB optogenetic site is exactly `left-Thalamus`; or
2. electrode `location` contains `brain_region: Thalamus` and the NWB optogenetic site is exactly `left-ALM`.

The session must contain at least one trial with `stim_present=True` and `photo_stim_period == "delay"`. No other recording region, stimulation target, or whole-cell session may be substituted after outcome extraction.

## Frozen observational predictor

The predictor is computed from **non-photostimulation trials only** and never uses the trial response.

- Eligible trial instruction labels: `lick left` and `lick right`.
- Neural window: the final 0.500 s before `cue_start_time`.
- Units: all source-exported extracellular units in the session; no post hoc cell-type filtering.
- Feature: spike count of each unit in the frozen window.
- Class label: instructed trial type (`lick left` vs `lick right`).
- Decoder: scikit-learn pipeline `StandardScaler` + L2 `LogisticRegression(C=1, solver="liblinear", max_iter=2000, random_state=20261012)`.
- Cross-validation: stratified shuffled folds with seed 20261012; 5 folds when both classes contain at least five trials, otherwise the number of folds equals the smaller class count; sessions with fewer than three trials in either class are ineligible.
- Session score: concatenated held-out ROC AUC across all eligible control trials.
- Mouse predictor: unweighted mean session AUC within its frozen cross-loop stratum.

A prespecified temporal-control predictor repeats the identical pipeline in the first 0.500 s after `pole_in_time`. It is descriptive and does not govern primary success.

## Frozen causal outcome

Both frozen perturbation sites are in the **left** hemisphere, so the prespecified contralateral instructed response is `lick right`.

For each eligible session:

- control trials: `type == "lick right"` and `stim_present == False`;
- perturbation trials: `type == "lick right"`, `stim_present == True`, and `photo_stim_period == "delay"`;
- a trial is scored correct iff the source-exported `response` string, after lowercase/strip, equals `correct`;
- session behavioral impairment = control correct fraction minus delay-photostimulation correct fraction;
- sessions require at least 10 control and 10 perturbation `lick right` trials with a nonmissing response.

Mouse outcome is the unweighted mean eligible-session impairment within its frozen cross-loop stratum. No normalization by control accuracy will replace this endpoint after reveal.

## Primary inferential statistic

Mouse is the independent unit. The two strata contain disjoint mice.

Within each stratum separately, mouse predictor values and mouse impairment values are converted to average ranks and centered by the stratum mean rank. The primary statistic is the Pearson correlation of the concatenated centered predictor ranks and centered outcome ranks across the two strata (a pooled within-stratum rank association).

Primary P value: 500,000 one-sided permutations, seed 20261012, permuting outcome ranks **within each stratum independently** while keeping predictor ranks fixed. P = (number of permuted statistics >= observed + 1)/(500000 + 1).

Uncertainty: 50,000 stratified mouse bootstrap resamples, resampling mice with replacement within each stratum and recomputing the pooled within-stratum rank statistic.

## Frozen held-out prediction test

Leave one mouse out. On each training set, fit ordinary least squares to raw mouse impairment using late-delay AUC plus a binary stratum indicator. Predict the held-out mouse without refitting on it. Report Spearman correlation of predictions with observed outcomes and R2 relative to a fold-specific baseline that predicts the training mean outcome of the held-out mouse's stratum.

## Prespecified secondary results

- exact one-sided Spearman permutation P values within each 7-mouse cross-loop stratum (all 7! permutations when all seven remain eligible);
- the temporal-control early-window predictor using the same pooled within-stratum statistic;
- difference between late-delay and early-window pooled rank associations, reported descriptively.

No secondary result can rescue a failed primary criterion.

## Success criterion

The primary cross-loop test is considered passed only if all conditions hold:

1. at least 12 total eligible mice and at least 5 eligible mice in each frozen direction;
2. pooled within-stratum rank association > 0;
3. one-sided within-stratum permutation P < 0.05; and
4. leave-one-mouse-out R2 > 0.

If structural eligibility falls below the minimum after response/spike extraction, the primary is reported as structurally ineligible rather than redefined.
