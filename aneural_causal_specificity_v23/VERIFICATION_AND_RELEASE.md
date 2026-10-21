# Aneural causal-specificity: Version 23 computational verification

This is a **post-registration, synthetic-data** verification note for the manuscript *Shared physiological drive can mimic cross-modal integration*. It is not evidence that any organism is conscious or that a living Stentor experiment has been carried out.

## Public code archive

The complete consolidated research source is preserved in this repository at [aneural_causal_specificity_v18/code](../aneural_causal_specificity_v18/code). In particular:
- `Prospective_Stentor_Validation_Analysis.py` is the original frozen three-gate analysis, expected SHA-256 `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`.
- `Supplementary_Code_6.py` generates the passive common-driver and no-change synthetic cases and runs the **unchanged** analysis against both.
- `Supplementary_Code_7.py` reproduces the Gaussian common-cause/directed-transfer/off-target intervention scenarios.
- `independent_audit.py` is the earlier independent reconstruction; the refreshed two-case independent checks are included in the local versioned submission package.

## Fresh rerun (21 October 2026)

A clean local execution from the versioned source included with the manuscript reproduced:

| Condition | BASE effective dimension | HAB effective dimension | Cross-modal AUC | Held-out mean log-score advantage | Frozen decision |
|---|---:|---:|---:|---:|---|
| Common-driver | 6.610755314 | 6.0613777 | 0.842122 | +0.679584 | 3/3 gates pass |
| No change | 6.610755314 | 6.612332191 | 0.493164 | -0.001868 | 0/3 gates pass |

An independent script recomputed the participation ratios from the generated `.npz` trajectories, cross-modal AUC from the feature matrices, and held-out predictive log scores from the cell-level output table. Numerical values matched the frozen summaries for both conditions. Four independent artificial day summaries yielded exact two-sided sign-flip P=0.125 for the common-driver geometry and predictive advantage; this is not a conventionally significant cluster-level biological result.

Run locally:

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pandas==2.2.3 scikit-learn==1.8.0 statsmodels==0.14.6 matplotlib==3.10.8
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg python aneural_causal_specificity_v18/code/Supplementary_Code_6.py --locked aneural_causal_specificity_v18/code/Prospective_Stentor_Validation_Analysis.py --outdir fresh_rerun
```

Regenerate rather than relabel the trajectories as organism observations. The repository includes **synthetic generator code**, not actual Stentor electrical traces.

## Novelty limits

Causal non-identifiability of common inputs, intervention exclusion restrictions, and theory-specific testing are established concepts. The added contribution is a *design-specific executable counterexample*: an explicit passive-driver generator satisfies a previously frozen three-gate statistical procedure. It does not show that all markers of consciousness are invalid or that Stentor is/is not conscious.

## Disclosure and reproducibility status

This new file is an archival verification note, **not** a new preregistration or experimental amendment. Prospective biological experiments and raw time-matched RNA-seq quantification remain unperformed. Research code and manuscript were developed with AI assistance and require final human author audit and approval. These repository files are public and addressable by immutable Git commit after this addition. A Zenodo DOI was not obtained. The prior v18 README records an older, incomplete negative-control verification; the fresh run documented above supersedes **that specific status claim**, not the original preregistered protocol.

Public biological sources belong to their original authors. Do not infer that their raw data are included in this archive.
