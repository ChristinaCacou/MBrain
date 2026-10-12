# Frozen external test of bilateral ALM coupling in Chen et al. 2024 / DANDI 000363

Frozen: 2026-10-12, before numerical perturbation-outcome extraction from DANDI 000363.

## Scientific question
Does the control-trial bilateral ALM coupling feature that was positive in the 2021 Chen raw reanalysis predict behavioral recovery during unilateral ALM photoinhibition in an independent experiment?

## Dataset
DANDI 000363, version 0.230822.0128 (Chen et al., Cell 2024). Only source NWB sessions are eligible. No result from the perturbation outcome may be used to alter eligibility, feature construction, time window, endpoint, direction, or statistical test.

## Eligibility
A session must contain:
1. simultaneous left-ALM and right-ALM spiking data;
2. source-defined quality-passed units in both ALM hemispheres;
3. correct non-photostimulation trials for both behavioral choices;
4. unilateral ALM photoinhibition trials with enough trials to estimate choice-balanced behavioral performance.

An animal contributes only if at least one session meets all criteria. Eligible sessions are averaged equally within animal. If fewer than 8 animals are eligible, the result is reported as underpowered/descriptive and cannot be labeled a successful external replication even if nominal P<0.05.

## Observational predictor
The feature is fixed to the successful 2021 raw-analysis construction as closely as the independent NWB schema permits:
- final 400 ms before the go cue / movement instruction;
- source quality-passed ALM units only;
- separate left- and right-ALM population choice directions trained on correct control trials;
- 30 deterministic split-halves, seed 20261012;
- within each split, train each hemisphere's coding direction on half of the correct control trials and project held-out correct control trials;
- compute Spearman correlation between held-out left- and right-ALM projections;
- average the 30 correlations to obtain session-level bilateral coupling.

If the independent NWB files do not contain pyramidal/FS labels, no post hoc waveform threshold will be invented: all source quality-passed ALM units will be used, and that deviation from the 2021 source feature will be stated explicitly.

## Primary perturbation outcome
Because the current manuscript does not establish hemisphere specificity in the 2021 Chen result, the primary endpoint is side-averaged unilateral ALM behavioral recovery rather than choosing a side after inspection.

For each unilateral stimulation side, behavioral performance is balanced accuracy across the two instructed choices. Recovery = photoinhibition balanced accuracy / matched control balanced accuracy. When both unilateral sides are estimable in a session, the primary session endpoint is their arithmetic mean. If only one unilateral side is present by experimental design, that session is ineligible for the side-averaged primary analysis but may contribute to the prespecified side-specific secondary analysis.

## Secondary outcomes
1. right-ALM behavioral recovery;
2. left-ALM behavioral recovery;
3. neuronal recovery along the corresponding control-trained ALM choice direction, if the NWB schema permits the same definition without adding outcome-informed preprocessing;
4. within-choice bilateral coupling as a control predictor.

Secondary results cannot substitute for the primary endpoint.

## Inferential unit and test
Animal is the inferential unit. The primary test is Spearman correlation between animal-mean bilateral coupling and animal-mean side-averaged behavioral recovery.

Inference:
- 500,000 one-sided animal-label permutations, seed 20261012;
- 50,000 animal bootstrap resamples for the 95% interval;
- leave-one-animal-out linear prediction with predictive Spearman rho, MAE, and R2 relative to a training-mean baseline.

Primary success requires all of:
1. rho > 0;
2. one-sided permutation P < 0.05;
3. held-out R2 > 0;
4. at least 8 eligible animals.

Failure of any component is failure of this external transport/replication test. No feature, time-window, side, endpoint, exclusion rule, or success criterion will be changed after numerical outcome extraction.

## Allowed schema adaptation before outcome extraction
Purely structural inspection may identify the NWB field names for: trial start/stop, go cue, instructed choice, correctness, photostimulation side, unit quality, electrode location and hemisphere. Such schema translation is allowed only before numerical perturbation outcomes are calculated and must be documented in code. It may not change the definitions above.
