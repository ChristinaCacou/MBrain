# Shared physiological drive can mimic cross-modal integration — reproducibility snapshot v18

**Status: computational verification complete, biological validation incomplete.**

Source provenance: frozen prospective algorithm SHA-256 `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`. The synthetic stress test is a **post-registration audit**, not a biological preregistration result. The author did not run new *Stentor* cultures in this package.

## Run the numerical analyses

Python >=3.10, `numpy`, `pandas`, `scipy`, `scikit-learn`, `statsmodels`, `matplotlib`.

```bash
python run_reproduction.py
```

This generates `central_reproduction/`, `causal_reproduction/`, and `AUDIT_NUMERICAL_RESULTS.json`. Execution reruns synthetic examples, including the unchanged frozen statistical pipeline. The biological observations are **not regenerated**, because the 13,546-row source dataset is maintained by the original authors. The source reproducer retrieves it and checks its recorded Git blob hash; run it separately with access to the public GitHub source.

Note: historical scripts are provided solely for traceability; they are not required to reproduce the focused paper's primary results. The RNA table is a secondary extract; incomplete time-matched untreated profiles mean no time-corrected forgetting claim.

## Key status distinctions

1. Numeric output re-executed from unchanged frozen code and independently rechecked against a separate set of equations and rank/eigenvalue computations. This is **internal independent-implementation verification**, not independent human peer review.
2. 96 synthetic cells on four synthetic days; common-drive test passed all three software gates while no-change test failed all three. These do not establish biological false-positive rates.
3. Monte Carlo intervention model: 4,000 campaigns in each of 3 synthetic mechanisms, not prospective Stentor intervention data.
4. The observed Stentor early/late contraction result comes from real previously recorded cells; the associated paper by Doan et al. remains the primary experimental contribution.
5. Publication venue, corresponding author declarations, funding, ethics, and AI editorial assistance require author approval.

## Integrity

All files have a SHA-256 manifest. Source snapshots are research artifacts; cite original biological data producers appropriately. Do not alter the publicly frozen prospective Stentor protocol based on this post-registration analysis. No DOI or external independent peer review is claimed.