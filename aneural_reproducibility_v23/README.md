# Aneural causal specificity: reproducibility snapshot v23

**Scope:** methodologic/theoretical manuscript, synthetic Gaussian models, and secondary analysis of external published biological recordings. No prospectively recorded Stentor observations or independent RNA-seq control quantification.

## One-command minimal verification

```bash
python -m pip install -r requirements.txt
bash verify_core.sh
```

The verification: (a) checks the *unmodified, preregistered* pipeline SHA-256 at runtime, (b) generates synthetic shared-driver and matched no-change cohorts, (c) runs all three frozen gates, (d) independently calculates participation ratios using singular values and exact four-day sign flips, and (e) independently verifies hypothetical negative-control proxy interval bounds. **It does not imply external peer verification or biological validity.**

## Results and interpretation

- Shared-driver synthetic cohort: 96 virtual cells over 4 modeled days. Frozen gates all pass; AUC approximately 0.8421; held-out log score difference approximately +0.6796. Matched no-change synthetic condition fails all three gates.
- Day-level exact two-sided sign flip for the four synthetic days: P = 0.125; cell-level gate P values are not independent-day replication evidence.
- SVD participation-ratio check versus frozen covariance eigenspectrum gives maximum absolute difference below 2e-14 for all 768 synthetic waveforms.
- Proxy bounds use hypothetical, uncalibrated parameters (rho=0.60, observed intervention contrasts e=0.25, m=1.0); no empirical inference about living Stentor is warranted.
- See `results/` for saved summaries. Additional ancillary scripts preserve v19-v22 stress tests.

## Source traceability

`code/Prospective_Stentor_Validation_Analysis.py` must hash to `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`, consistent with the existing prospective preregistration. It is not changed.

The living-cell observational source is Doan et al., `dataset_followup_2.csv` (original Git blob `04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc`). Third-party biological observations are NOT republished here; obtain at https://github.com/nhithuydoan/associative_learning_Stentor and comply with source terms.

RNA author releases: https://github.com/aralbright/2025_stentor_learning and https://doi.org/10.5061/dryad.z34tmpgrm . Six untreated recovery-time profiles were not quantified here. The independent Miao RNA-seq count matrix was not analyzed here.

## Limitations and release status

The source-code archive documents exactly one reproducible computational argument within a restricted model family. This is not proof of subjective experience, no subjective experience, integrated information theory, or general invalidity of neural PCI. The public GitHub commit is a versioned snapshot, not a DOI-minted immutable deposit. No license for original scripts is asserted pending explicit author review.
