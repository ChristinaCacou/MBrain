# Shared physiological drive can mimic cross-modal integration — reproducibility release v18

**Content types:** synthetic data, secondary biological summaries derived from cited public material, scripts, audit reports. No new biological measurements or unpublished human data.

## 1. Environment
Python 3.10+ with `numpy`, `scipy`, `pandas`, `matplotlib`, `scikit-learn`, `statsmodels`, `python-docx` (only to regenerate documents). Simulation scripts are standalone and use fixed RNG seeds.

## 2. Original three-gate simulation
`python src/Supplementary_Code_6.py --locked src/Prospective_Stentor_Validation_Analysis.py --outdir results/fresh_common_driver`

This checks SHA256 of the frozen original analytic program `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`, generates 96 virtual cells across four simulated days, and runs the original analysis on both a shared physiological driver and a no-change control. It is post-registration exploratory simulation, **not** observed Stentor data or a completed prospectively registered experiment.

## 3. Independent numerical recomputation
`python src/independent_verification.py`

The independent implementation recomputes dimension from stored trajectories in `data/common_driver` and `data/no_change` (for its original relative-path expectation, copy/symlink `data/` to `frozen_audit/`; see `verify_release.py`). The published audit finds a common-driver reduction from 6.610755 to 6.061378 and an exact two-sided four-day sign-flip P=0.125; the cell-level P value is far smaller due to nested pseudo-independent observations. An in-workflow second implementation is not independent human peer review.

## 4. Analytic intervention benchmark
`python src/Supplementary_Code_7.py --outdir results/fresh_causal_audit`

This generates synthetic observationally equivalent structures and intervention outcomes. Expected observed rejection rates under its fixed-seed design: common cause 0.05075; directed transfer 0.81375; off-target common cause 0.69225. These are design-specific Monte Carlo frequencies, not empirical false positives.

## 5. Independent biological dataset
`src/Supplementary_Code_16_source_reproducer.py` checks the original public Doan et al. dataset and regenerates source-derived run/date analyses. Network access is needed to retrieve `https://github.com/nhithuydoan/associative_learning_Stentor/blob/main/data/dataset_followup_2.csv` (Git blob 04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc). The reduced `data/ExtendedData_Table13*` summaries are reanalyses of source data, **not new wet-lab observations**.

## 6. Provenance & limitations
The original prospective Stentor protocol (public commit `07949d983493186d63d3232a546247476adae4bf`) remains unchanged and **unexecuted**. The larger 21-sample Rajan raw RNA-seq validation is **not** included: six time-matched untreated controls were not quantified. The separate independent Miao cohort was not quantified gene by gene. Figures in documents were generated from numerical inputs, not AI text-to-image generation.

## 7. Claims and release policy
The manuscript's novelty is limited to an executable adversarial failure of a **specified restricted comparator workflow**. Observational equivalence, pseudoreplication, and broad skepticism about plant-consciousness indicators are established prior art. All scripts were created/edited with assistance from OpenAI ChatGPT and require human scientific verification before submission. No external expert review is claimed. See `docs/Author_Submission_Checklist_v18.docx` for mandatory author signoff.

## Archive integrity
Every file listed in `MANIFEST_SHA256.txt` is hash-pinned. Do not include the separate private blinded-treatment key in a public source release.
## Public GitHub branch execution

Clone branch `aneural-causal-specificity-v18-reproducibility`, then run:

```bash
cd research/aneural-v18
python -m pip install numpy scipy pandas matplotlib scikit-learn statsmodels
python run_public_replication_v18.py
```

This command generates both synthetic cohorts into `src/frozen_audit/`, checks the frozen source hash, reruns the same three-gate analysis, recomputes effective dimensionality without importing that analysis, and runs the causal-intervention benchmark. The saved summaries under `results/` are source-audit receipts. Raw generated trajectories are reproduced by code rather than checked into the Git branch. The complete offline ZIP, including saved trajectories and immutable manifest, is provided as a separate manuscript reproducibility package. The Git branch is public but is **not itself peer review**.
