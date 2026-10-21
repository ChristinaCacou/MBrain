# Aneural causal specificity: consolidated computational archive (v18)

Research code and verification for the manuscript **Shared physiological drive can mimic cross-modal integration**. This is a *computational and secondary-data* project. **No new living-cell experiments, direct consciousness measurements, or prospective Stentor confirmatory biological recordings were performed.**

## Reproducibility entry points

- `code/Supplementary_Code_6.py`: generate the adversarial passive-common-driver and negative-control synthetic cohorts; run the frozen three-gate Stentor analysis.
- `code/Prospective_Stentor_Validation_Analysis.py`: byte-for-byte copy of the originally locked analysis script, SHA-256 `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`. Its prospective protocol/thresholds are **not** amended by this archival copy.
- `code/independent_audit.py`: separate numerical verification of the positive-case geometry, cross-modal AUC, held-out Gaussian log score, negative-control geometry, and day-level exact sign-flip reference values.
- `code/Supplementary_Code_7.py`: synthetic observational equivalence and intervention-selectivity experiment.
- `code/Supplementary_Code_16_source_reproducer.py`: biological habituation reanalysis from author-released `dataset_followup_2.csv` pinned Git blob `04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc`.
- `code/Supplementary_Code_1.py`, `2.py`, `3.py`, `5.py`: earlier ancillary synthetic/robustness analyses; historically exploratory.

## Commands

Create Python environment containing `numpy pandas scipy scikit-learn statsmodels matplotlib` and `python-docx` if rebuilding manuscript files.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg \
python code/Supplementary_Code_6.py \
  --locked code/Prospective_Stentor_Validation_Analysis.py --outdir replication
python code/independent_audit.py
```

The independent verifier expects `replication/common_driver/locked_outputs/cell_state_metrics.csv`, and both cases' generated synthetic trajectory folders in the same parent as the script. This repository stores the source and an audit summary, not the approximately 22 MB temporary synthetic NPZ trajectory corpus; regenerate locally.

**Audit status:** the positive synthetic case completed under the unchanged frozen analysis and reproduced all previously reported scores exactly. The independent verifier matched dimensionality means 6.6107553 and 6.0613777, cross-modal AUC 0.8421224 and held-out mean Δ log score 0.6795836. The independent negative-control trajectory geometry means 6.6107553 and 6.6123322 with paired t P=0.9650. A complete second frozen negative-control run timed out in the available session. Do **not** describe it as fully independently rerun.

**Critical inferential qualification:** the original frozen primary P=2.289e-17 and its predictive-model sign-flip P=0.00050 use **cells**, whereas synthetic data contain only **4 recording-day units**. Recalculation of two-sided *exact day-level sign-flips* gives P=0.125 for geometry and predictive log score. The frozen gates can pass without proving robust day-level significance; this is a model-design vulnerability, not a confirmed empirical false-positive rate.

## Public datasets and provenance
- Doan et al. video-scored behavioral data: https://github.com/nhithuydoan/associative_learning_Stentor
- Rajan et al. original analyses: https://github.com/aralbright/2025_stentor_learning
- Rajan raw archives: https://doi.org/10.5061/dryad.z34tmpgrm
- Original Stentor preregistration: https://github.com/ChristinaCacou/MBrain/commit/07949d983493186d63d3232a546247476adae4bf

The original RNA-seq recovery effect lacks quantification of six archived untreated time-matched controls. This archive does not claim to have reprocessed them or the independent Miao 21-library cohort.

## Scientific status and limitations
This is a constructive counterexample for a specified statistical validation design, not an inference about the presence or absence of phenomenal experience. Shared physiological drive and directed information transfer are not mutually exclusive in all conceivable systems, and neither mechanism is directly identified in Stentor by the available datasets. The source code needs external investigator review, environment pinning and archive DOI issuance for full journal-level computational reproducibility. The scripts were historically generated/refined with AI assistance; human authors must review and certify all analyses.

Version v18 is a post-registration archival and interpretation update. It does not modify the original preregistered Stentor intervention.
