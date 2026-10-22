# Aneural causal specificity audit — version 18 (21 October 2026)

This is the code and numerical-fixture release supporting **“Shared physiological drive can mimic cross-modal integration: a causal specificity audit.”** The accompanying Word manuscript, technical supplement, larger historical code archive and declaration checklist are distributed separately in the ChatGPT Library. **This is a public source snapshot, not a peer-reviewed study, independent biological validation, or a DOI-minted release.**

## Core files
- `Supplementary_Code_6.py`: generates two virtual response pathways driven by a shared exogenous gain (no mechanical-to-electrical transfer) and reruns the **unchanged, pinned** three-gate analysis.
- `Prospective_Stentor_Validation_Analysis.py`: frozen original analysis; SHA256 `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`. The published Git blob matches the original byte sequence.
- `Supplementary_Code_7.py`: observationally equivalent Gaussian models and 12,000 synthetic intervention campaigns.
- `Supplementary_Code_16_source_reproducer.py`: independent secondary count analysis from Doan et al. public **video-scored** `Stentor` dataset, with source Git blob validation.
- `results/ExtendedData_Table8_common_driver_gates.csv`: original locally reproduced comparison.
- `results/ExtendedData_Table9b_intervention_assurance.csv`: original locally reproduced intervention summary.

## Reproduce
Use Python 3.10+ with `numpy pandas scipy scikit-learn statsmodels matplotlib` in a clean environment. Download these source files in one directory, retaining `results/` for the fixtures.

```bash
python Supplementary_Code_6.py --locked Prospective_Stentor_Validation_Analysis.py --outdir work/common_driver
python Supplementary_Code_7.py --outdir work/causal_interventions
```

Expected: the common-driver case passes all three analysis gates, the no-change case fails all three. The observational intervention simulations yield rejection fractions 0.05075 (common cause), 0.81375 (directed transfer), 0.69225 (off-target shared cause). These are **operating characteristics of artificial generative models**, not empirical false-positive rates for biological consciousness.

To reanalyse independent `Stentor` behavioral scores:

```bash
python Supplementary_Code_16_source_reproducer.py --output work/behavior
```

This retrieves and checks `https://raw.githubusercontent.com/nhithuydoan/associative_learning_Stentor/main/data/dataset_followup_2.csv`, pinned to Git blob `04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc`. Internet access is required for that source unless `--csv` is supplied.

## What this does *not* establish
The common-driver model was *selected* to challenge a rule, not calibrated to observed cell biology. The retrospective tests are not an independent preregistered prospective Stentor experiment. Six untreated RNA-seq control profiles have not been quantified and independent Miao gene-level results have not been reprocessed here. Distinct response modules sharing a physiological influence need not be truly causally independent given that shared state, and a non-phenomenal computational description cannot prove an organism lacks experience.

## Provenance
The code was reproduced locally (not externally independently reviewed) against saved results with differences below 1e-10. Original prospective analysis SHA-256 checking is enforced by Supplementary Code 6. All supplementary model codes not needed for the main computational claim and a complete zipped source/data/checksum bundle accompany the offline manuscript package.

## DOI / archiving
This GitHub directory has no DOI. For durable citation, archive a tagged release using Zenodo or a comparable repository and add the accession to the manuscript *before submission*. Do not cite this as a registered new prospective experiment.
