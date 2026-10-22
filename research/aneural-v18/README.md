# Shared physiological drive can mimic cross-modal integration
## Version 18 — public computational reproducibility snapshot

This folder is part of a PUBLIC source branch. All new simulation datasets are SYNTHETIC. There are no new Stentor biological measurements, and the registered prospective Stentor experiment remains unperformed.

### Reproduce the exact executable counterexample
Clone the `aneural-causal-specificity-v18-reproducibility` branch of this repository and run:

```bash
cd research/aneural-v18
python -m pip install numpy scipy pandas matplotlib scikit-learn statsmodels
python run_public_replication_v18.py
```

The runner generates 96 virtual cells across four simulated days for both common-driver and no-change conditions, verifies the unchanged analysis script against SHA-256 `7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151`, runs the original three-gate evaluation, recalculates dimensionality from saved synthetic trajectories with an independent numerical implementation, applies a four-day exact sign-flip sensitivity test, and repeats a separate 12,000-campaign causal-intervention simulation. Generated trajectories are *not* tracked in this Git branch, since they are deterministically rebuilt.

### Selected expected numerical outputs (synthetic only)

| Metric | Common-driver result |
|---|---:|
| Electrical participation ratio, BASE | 6.610755 |
| Electrical participation ratio, HAB | 6.061378 |
| Cross-modal AUC | 0.842122 |
| Leave-day-out model advantage | 0.679584 |
| Frozen cell-level paired P | 2.289×10⁻¹⁷ |
| Four-day exact two-sided sign-flip P | 0.125 |
| No-change cell-level P | 0.965 |

The original three-gate procedure declares its *restricted scalar-gain null* rejected for the common-driver model and not rejected for the no-change control. Because there are only four independent simulated recording days, the exact sign-flip test cannot provide a two-sided P below 0.125. This is not biological validation and does not establish conscious or unconscious experience.

The selective-intervention benchmark gives rejection rates 0.05075 under a valid common-cause null, 0.81375 under directed transfer, and 0.69225 under common cause with an off-target perturbation. These are MODEL-SPECIFIC Monte Carlo operating characteristics, not diagnostic rates in Stentor.

### Contents

- `src/Prospective_Stentor_Validation_Analysis.py` — unchanged hash-locked analysis from prior public preregistration
- `src/Supplementary_Code_6.py` — common-driver and no-change generation plus frozen gate evaluation
- `src/Supplementary_Code_7.py` — observational equivalence and intervention simulation
- `src/independent_verification.py` — separate numerical participation-ratio implementation and exact four-day sign-flip check
- `src/Supplementary_Code_16_source_reproducer.py` — published Stentor contraction data reanalysis, requires internet for public source CSV
- `results/*.csv` — cached synthetic numerical results
- `run_public_replication_v18.py` — reproducible one-command driver

The full **offline author package** separately contains source files, generated trajectories, trial tables, figures, manuscripts, SHA-256 manifest and author-signoff documents. That ZIP is not hosted in this Git branch.

### Data provenance and limits

Original Stentor data are from Doan and colleagues' public repository [`dataset_followup_2.csv`](https://github.com/nhithuydoan/associative_learning_Stentor), Git blob `04e271fab29e9f826f7d3a9dca9c4e9040c1dbfc`. The original prospective study has prior public registration commit `07949d983493186d63d3232a546247476adae4bf`. The separately published Rajan RNA-seq tables do NOT contain the six time-matched untreated 10- and 90-minute controls used for molecular forgetting identification; those controls were not quantified in this project. Independent gene-level Miao replication is also unperformed.

### Scientific novelty and attribution

Observational equivalence, the need for selective interventions, pseudoreplication, and broad criticism of behavior-only consciousness indicators are established prior art. This project claims, narrowly, to construct an executable counterexample to mechanistic interpretation of a **particular predeclared three-gate workflow** with an insufficient comparator class. It does NOT prove that all markers of consciousness fail, or that Stentor is phenomenally unconscious.

OpenAI ChatGPT assisted with source code, mathematical checking, writing and data processing. A separate computational reimplementation is not outside expert peer review. Human authors remain responsible for independently verifying the scientific content and making any journal submission. Do not redistribute a private blinded-treatment key.
