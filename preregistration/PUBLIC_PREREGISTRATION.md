# Public preregistration: prospective Stentor perturbation study

**Title:** Energy-matched multimodal perturbation of *Stentor coeruleus* to discriminate local adaptation from a global integrative-state model  
**Author:** Christina Cacoulidis, University of Missouri-Kansas City  
**Protocol version:** 1.0  
**Lock date:** 16 October 2026  
**Public timestamp:** the GitHub commit containing this file is the public time stamp for this preregistration snapshot.

## Scope and claim boundary

This is a prospective model-discrimination experiment, not a direct test of phenomenal consciousness. The experiment asks whether mechanical habituation in *Stentor coeruleus* is adequately explained by a local scalar change in response gain/threshold, or whether it leaves a modality-general signature in the geometry of an orthogonal evoked response after response magnitude is removed.

A positive result may reject one explicit non-conscious null model. It must not be reported as demonstrating that *Stentor* is conscious.

If confirmatory observations were collected before this public commit outside this workflow, this deposit must be treated as retrospective rather than pre-data registration.

## Frozen competing models

### N0: local adaptive-gain/threshold model

For cell i, state s and modality m,

`Y_i,s,m(t) = a_i,s,m F_i,m(t) + epsilon`.

Habituation may change response probability, scalar gain, threshold and latency in the mechanically trained pathway, but after baseline subtraction, rigid alignment and unit-energy normalization the modality-specific response geometry should remain stable. N0 therefore predicts no reproducible BASE-to-HAB change in normalized electrical-probe geometry and no above-chance cross-modal state decoding from mechanical to electrical responses.

### N1: global integrative-state alternative

N1 adds a shared latent state `z_i,s`:

`Y_i,s,m(t) = a_i,s,m F_i,m(t) + L_m z_i,s(t) + epsilon`.

N1 predicts a BASE-to-HAB change in unit-energy electrical response geometry, above-chance transfer of a mechanically trained state decoder to electrical trials, and a directional return toward baseline after dishabituation/recovery.

## Primary endpoint

The primary endpoint is the within-cell BASE-to-HAB change in effective dimensionality of the amplitude-normalized electrical-probe response trajectory in SHAM-treated cells. Effective dimensionality is the covariance participation ratio:

`(sum lambda)^2 / sum(lambda^2)`.

Trajectories are baseline-subtracted and divided by the square root of total response energy before geometry metrics are computed, so a pure scalar-amplitude effect cannot by itself generate the confirmatory result.

## Three conjunctive confirmatory gates

The local scalar-gain null is rejected only if **all three** gates pass.

1. **Primary geometry gate.** Two-sided paired P < 0.05 and absolute paired standardized effect |d_z| >= 0.35 for SHAM HAB minus BASE electrical effective dimensionality. If the 90% confidence interval lies wholly within +/-0.20 SD and paired TOST P < 0.05, the result is instead interpreted as equivalence at the registered margin.

2. **Cross-modal decoder gate.** A fixed L2 logistic regression (C = 1) is trained only on unit-energy MECH BASE versus MECH HAB features and then applied without refitting to ELEC BASE versus ELEC HAB. The gate requires a cell-cluster bootstrap lower 95% confidence bound for AUC > 0.55 and a cell-level permutation P < 0.05. Response energy is excluded from the decoder.

3. **Frozen predictive-model gate.** In leave-recording-day-out analysis, N0 predicts zero normalized electrical BASE-to-HAB feature shift. N1 predicts transfer of the mechanical BASE-to-HAB shift learned on training days to held-out electrical responses. The gate requires positive mean held-out log-score difference for N1 versus N0, a day-cluster bootstrap 95% confidence interval excluding zero, and sign-flip P < 0.05.

Failure of any gate leaves N0 viable for the registered protocol.

## Sample size and replication

Target enrollment is 96 cells per treatment arm, 192 cells total, blocked as 24 cells per arm on each of four independent culture/recording days. The SHAM arm is primary. The pre-specified causal secondary arm is SteCoe_6763 RNAi, provided the exact validated reagent/procedure is available before acquisition; otherwise that arm must be omitted rather than replaced post hoc.

For a paired standardized effect d_z = 0.40, approximately 68 analyzable pairs provide 90% power at two-sided alpha = 0.05. Enrollment of 96 per arm allows approximately 25% pre-specified attrition. There is no efficacy stopping rule and no interim significance testing.

Culture/recording day is the biological replicate block. Cells are the primary paired unit, with uncertainty additionally estimated by resampling recording days.

## Experimental sequence

A separate non-confirmatory calibration cohort of up to 20 cells may be used to fix hardware settings before confirmatory acquisition. Calibration cells never enter confirmatory analyses.

For each confirmatory cell:

- baseline unstimulated recording followed by randomized-order low-force mechanical and orthogonal electrical probes separated by at least 3 min;
- electrical probe fixed at 0.5 s and 6.5 V unless a pre-data safety/viability amendment is publicly timestamped;
- high-speed bright-field video from 5 s before to 10 s after each probe, recommended >=100 frames/s;
- low-force mechanical training at 60-s intervals;
- habituation criterion: three consecutive failures to contract after at least three initially responsive training trials, with a maximum of 30 training stimuli;
- immediately after criterion, one electrical HAB probe;
- after 3 min, a fixed high-force mechanical dishabituating pulse followed 30 s later by an electrical DISH probe;
- after 30 min without stimulation, an electrical REC probe followed by one low-force mechanical recovery test.

A cell that fails to reach habituation is coded HAB_NOT_REACHED and is not excluded from the intention-to-measure dataset.

## Locked trajectory representation

Each frame is segmented to recover the body contour, resampled to 32 arc-length landmarks, centroid-translated and rigidly aligned to a common oral-aboral axis. The 64-dimensional x/y landmark vector is interpolated to 100 equally spaced post-stimulus time bins. Segmentation must not use condition labels.

Locked confirmatory metrics are:

- effective dimensionality (primary);
- Lempel-Ziv phrase complexity;
- zlib compressed complexity;
- response energy, reported descriptively but excluded from cross-modal decoding.

No post hoc feature engineering may replace or redefine the primary endpoint.

## Randomization, blinding and exclusions

Within each recording day, treatment allocation is balanced and randomized before imaging. Baseline MECH/ELEC order is independently randomized. Cell identifiers contain no condition information.

Where feasible, acquisition is treatment-blinded. Segmentation and confirmatory analysis use coded treatment/state labels. The treatment key is revealed only after the frozen analysis script has completed integrity checks and generated the confirmatory outputs.

A trial may be excluded only for: detachment, leaving the field of view, cell death, segmentation failure, stimulus-delivery failure, or pre-defined image-quality/QC failure. Failure to habituate is an outcome, not an exclusion. Statistical outliers are not excluded for weakening the predicted result.

## Statistical analysis

The report must include paired means, paired differences, paired d_z, day-cluster bootstrap 95% confidence intervals, paired t-test P, Wilcoxon signed-rank P, and the registered equivalence test at +/-0.20 SD.

Cross-modal decoding uses exactly three unit-energy features: effective dimensionality, Lempel-Ziv complexity and compressed complexity. Standardization and classifier coefficients are fit only to MECH BASE/HAB data before application to ELEC BASE/HAB.

The frozen predictive score is computed leave-recording-day-out with no hyperparameter search. Secondary endpoint P values are controlled with Benjamini-Hochberg FDR at q = 0.05. Exploratory analyses cannot rescue a failed confirmatory gate.

## Frozen computational implementation

The confirmatory implementation is `Prospective_Stentor_Validation_Analysis.py` with fixed seed 260416, fixed metrics, fixed classifier and fixed decision thresholds. It includes a synthetic software self-test only; synthetic data are not biological evidence.

Each biological probe trial is represented by one metadata row and one `.npz` file containing a two-dimensional array named `trajectory`. The script rejects unknown codes, duplicate identifiers, non-preregistered exclusion reasons, missing included trajectories and non-finite values.

## Cryptographic lock manifest

The following SHA-256 digests identify the exact files frozen before this public registration:

```text
05ed3f2bb4284ef2c63227362528bc4a16487644482418d7f42fd5de72c34c05  Aneural_Adaptive_Dynamics_Consciousness_Article.docx
f2c7c0b7c329a45d532a6c96bcae6417a5a4f9cc4772f28928b1ab832bbc308c  Stentor_Prospective_Perturbation_Preregistration.docx
7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151  Prospective_Stentor_Validation_Analysis.py
a6b08df58cc94973659099396cf48db7c01c6f16adb86eabfb2ad4b1a7b57ab2  Stentor_Validation_Data_Template.csv
```

Synthetic software dry-run output frozen before registration:

```text
08c8ca533d59656511113480c4fee39d8e567a0f4b787d2f578100efe29bca66  prospective_selftest/synthetic_results/confirmatory_results.json
```

Any post-registration change to the protocol, thresholds, analysis code, sample plan or exclusions must be accompanied by a separately timestamped amendment that states whether confirmatory outcomes had been viewed.

## Interpretation matrix

- **All three gates pass:** reject the registered local scalar-gain null for this protocol; favor the global integrative-state alternative for measured response geometry. Do not infer consciousness.
- **Primary gate fails and equivalence passes:** normalized electrical geometry is approximately invariant at the registered margin; N0 remains viable.
- **Primary differs but cross-modal gate fails:** a state-dependent electrical change exists, but the registered evidence for modality-general access is absent; N0 is not rejected on the full criterion.
- **Cross-modal decoding passes but predictive-model gate fails:** classifier performance alone is insufficient; N0 remains viable.
- **RNAi modifies behavior but not cross-modal geometry:** molecular learning machinery affects habituation without the registered global-state signature.
- **RNAi modifies both behavior and cross-modal geometry:** supports a causal link between known learning machinery and the registered distributed state signature, but still does not establish consciousness.

## Registration requirements before biological acquisition

Before the first confirmatory cell, record the exact low-force and high-force actuator settings from the non-confirmatory calibration cohort, confirm the exact SteCoe_6763 RNAi and SHAM procedures or omit that arm, record microscope/camera model and frame rate, generate the treatment randomization, preserve the blinded code key separately, and archive the locked analysis self-test output.

The frozen full protocol is identified above by SHA-256. The purpose of this public commit is to establish an externally visible timestamp for the registered design and its cryptographic file identities before confirmatory analysis is represented as prospective.