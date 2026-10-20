# Aneural causal-specificity study — reproducibility snapshot v18

**Status:** computational and secondary-data research draft; NOT peer reviewed; no new *Stentor* wet-lab data were collected.

This snapshot consolidates the executable code underlying *Shared physiological drive can mimic cross-modal integration* (v18). It documents an adversarial demonstration that the **particular frozen three-gate Stentor-oriented analysis** may favor an integration-like model when two response modules are influenced by a shared exogenous physiological variable, with **no directed path between modules**. This is **not** a claim about the validated human TMS–EEG Perturbational Complexity Index, nor proof that any biological organism lacks consciousness.

## Core reproduction

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python verify_release.py --root . --out ./audit_run
```

The reproducibility runner:
1. verifies SHA-256 of `code/Prospective_Stentor_Validation_Analysis.py` is `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151` (original prospectively frozen code);
2. runs the original synthetic common-driver and no-change negative-control benchmark (`code/Supplementary_Code_6.py`);
3. runs the observational equivalence / intervention analysis over 12,000 synthetic campaigns (`code/Supplementary_Code_7.py`);
4. checks numerical consistency with the reported synthetic metrics, and independently checks the covariance identity of the two structural models.

For reanalysis of independent public living-cell observations, see `code/Supplementary_Code_16_source_reproducer.py`. The 2026 Doan et al. source is `https://github.com/nhithuydoan/associative_learning_Stentor` (`data/dataset_followup_2.csv`; Git blob `04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc`). Refer to authors' repository license for reuse; original raw biological data are **not** included here.

Other historical illustrative simulations are included as Supplementary Codes 1, 2, 3 and 5. The exploratory analyses are not the preregistered biological study. Any future biological analysis must retain the original preregistration and never treat synthetic data as real observations.

## Verified results (deterministic seeds)

* Synthetic common driver: geometry gate P = 2.289466e−17; mechanical-to-electrical state AUC = 0.8421223958; leave-day-out model log-score advantage = 0.6795835992; all three software gates pass.
* Matched synthetic no-change control: combined decision rule does not pass.
* Independent repeated seeds 2026101710 and 2026101720: common-driver passes and matched no-change does not pass the combined decision rule in either run. This is a limited random-seed sensitivity check, **not** independently sampled biological validation.
* Intervention simulation: pure common cause 5.075%; directed transfer 81.375%; off-target common cause 69.225% of 4,000 synthetic campaigns per condition rejected at nominal 5%.
* External *Stentor* behavior: observed first-five contraction rate 741/2590 (28.61%), last-five 331/2598 (12.74%), across eight dates. This demonstrates repeated-stimulus response attenuation, not consciousness or directed cross-modal transfer.

## Caution and novelty

The literature already notes that observational complexity and apparent integration can be confounded by correlated input, common drivers and noise. This study's narrower contribution is an **executable counterexample to a specific multistage cross-modal scoring pipeline** and a worked exclusion-restriction failure. The mathematical principle that observationally equivalent causal models may differ under intervention is not new (Pearl, 2009; Hernán & Robins, 2020). Clinical human PCI uses direct cortical perturbation and should not be collapsed into our bespoke synthetic gate battery. See Sarasso et al., 2021, *Neuroscience of Consciousness*, doi:10.1093/nc/niab023; Hansen, 2024, *Biology & Philosophy*, doi:10.1007/s10539-024-09953-1.

## Known limitations

* Four simulated recording days for the three-gate result; the cell-level P value does not constitute a robust four-cluster inferential guarantee.
* The adversarial generator was chosen to demonstrate possibility, not sampled to estimate a model-class false-positive rate.
* Reanalysis of previously acquired living-cell responses does not directly validate the synthetic common-drive architecture.
* The published processed RNA-seq matrix omits untreated recovery-time control pools; no raw RNA-seq recovery quantification was completed.
* No human-independent audit of source code or claims is implied by a rerun using the same scripts.

## Author and disclosure status

Before submission the author must verify contribution roles, funding, competing financial/nonfinancial interests, full references, corresponding email, permissions and disclosure of generative-AI assistance. No such factual declaration is inferred by this archive.

## Version and integrity

This branch is a versioned public Git snapshot. Use a **full commit permalink**, not only the mutable branch URL, in the final submitted manuscript. SHA-256 manifest and numerical verification receipt accompany the downloadable package. The original prospective study registration is distinct: `07949d983493186d63d3232a546247476adae4bf`.
