# Aneural causal specificity: consolidated reproducibility snapshot (v18)

Public sources and summary outputs for the article *Shared physiological drive can mimic cross-modal integration*.

## What the archive does and does not establish

- `Supplementary_Code_6.py` synthesizes 96 virtual SHAM-labelled cells over four recording days, with local mechanical and electrical response modules receiving a common exogenous gain state but no modeled M-to-E or E-to-M pathway. It calls the **byte-identical frozen analysis** and regenerates the three gate values and matched no-change control.
- `frozen_Prospective_Stentor_Validation_Analysis.py` is the previously registered pipeline and is NOT modified for this archive. SHA-256: `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`.
- `independent_audit_v18.py` independently analyzes saved synthetic features, estimating independent-day summaries and mechanical-train/electrical-test leave-one-day-out diagnostics. Its day-level exact sign-flip tests are exploratory.
- `Supplementary_Code_7.py` tests observationally matched Gaussian common-cause and directed-transfer models in 12,000 synthetic campaigns, including an off-target-intervention failure case.
- `Supplementary_Code_16_source_reproducer.py` recalculates independent public *Stentor* contraction data and verifies the original Git object SHA: `04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc`.
- `stentor_public_reanalysis.py` reuses authors' *processed* RNA-seq counts and reported differential-expression tables, with Git blob SHA checks. It does NOT quantify original FASTQs.
- `Supplementary_Code_10.py` performs sensitivity analysis for missing time-matched untreated recovery pools using the extracted `previous_14genes.tsv`. It does NOT infer unmeasured control profiles.
- Codes 1, 2, 3 and 5 preserve historical supporting synthetic simulations. There is no representation here of new wet-lab measurements or completed prospective Stentor validation.

## Run (offline synthetic pipeline)

```bash
python -m pip install -r requirements.txt
python reproduce_all.py --out results
```

This regenerates the common-driver and no-change cases, validates the frozen SHA-256, performs an independent four-day diagnostic, and runs the selective-intervention benchmark. Approximately 0.4 GB+ disk and several minutes may be useful; results are deterministic on the pinned environment but exact floating point text can vary slightly across platforms.

If the original `dataset_followup_2.csv` is downloaded from the Doan et al. repository, run:

```bash
python reproduce_all.py --out results --behavior-csv path/to/dataset_followup_2.csv
```

For Rajan et al. processed expression matrices, see `stentor_public_reanalysis.py --help` and use files named exactly as specified there. Quantification of six omitted recovery control libraries and independent Miao cohort reads is **not** completed.

## Interpretation constraints

1. The original *cell-level* gate P values are not the same as strong inference over four recording days. With four day summaries the minimal one-sided exact sign-flip P value is 1/16 = 0.0625.
2. A classifier discriminating baseline from posttraining responses cannot uniquely identify directed information transfer when both pathways share a driver.
3. The intervention simulation is a conceptual benchmark only. Intervention specificity and exclusion restrictions must be validated in any biological application.
4. Public behavior is single-cell scored data, but it is from an independent study and is **not** the same experiment as the Rajan molecular recovery analysis.
5. RNA-seq analyses use previously published **pooled** sample profiles, not new recordings, and cannot establish molecular forgetting without untreated time-matched controls.
6. The code describes modeled causal structures; the absence of a consciousness variable in a model is not a proof that any corresponding physical system lacks experience.

## Upstream sources and author rights

- Doan et al. source: https://github.com/nhithuydoan/associative_learning_Stentor (raw published `dataset_followup_2.csv` is **not** mirrored here).
- Rajan/Albright processed source: https://github.com/aralbright/2025_stentor_learning.
- Raw Rajan Dryad archive: https://doi.org/10.5061/dryad.z34tmpgrm.
- Related Miao dataset: NCBI BioProject PRJNA1489669 (not independently reprocessed here).
- Preregistered analysis lock: https://github.com/ChristinaCacou/MBrain/commit/07949d983493186d63d3232a546247476adae4bf.

Private treatment assignment keys, personal correspondence and unpublished wet-lab data are not included.

## Provenance and reporting

`independent_audit_v18.json` is output of exploratory checks. The package is not a claim of peer-reviewed independent replication. Some methods and manuscript text were assisted by OpenAI ChatGPT, with author verification required. Funding, conflicts and final attribution need corresponding-author sign-off.
