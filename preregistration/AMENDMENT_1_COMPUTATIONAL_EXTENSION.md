# Amendment 1: post-registration computational extension

**Date:** 16 October 2026  
**Original public preregistration:** commit `07949d983493186d63d3232a546247476adae4bf`  
**Prospective protocol version:** 1.0, unchanged

## Purpose

This amendment records a post-registration strengthening of the computational/theoretical manuscript. It does **not** change the prospective *Stentor coeruleus* confirmatory protocol, its sample plan, acquisition sequence, exclusions, endpoints, thresholds, three conjunctive decision gates, randomization/blinding plan, or the locked `Prospective_Stentor_Validation_Analysis.py` identified in the original registration.

No confirmatory biological outcomes were added or represented as collected in this amendment.

## Added computational analyses

1. **Correlation-adjusted marker multiplicity.** For `k` standardized Gaussian markers with equal marginal mean shift and exchangeable pairwise correlation `rho`, the analytic effective independent-marker count is `k_eff = k/[1+(k-1)rho]`. This formalizes why nominally multiple correlated consciousness-linked markers need not provide independent evidence.

2. **Cross-topology regime classification.** A fixed L2 logistic classifier was trained only on four energy-excluded perturbational features from 30 paired Erdos-Renyi networks and applied without refitting to 30 paired Watts-Strogatz and 30 paired Barabasi-Albert networks. The analysis tests whether the marker bundle transfers as a classifier of a generic dynamical gain regime rather than as a consciousness-specific signature.

3. **Stable linear recurrent null.** One hundred stable sparse linear recurrent systems, with no nonlinear activation and no slow adaptation variable, were tested at recurrent gains 1.0 and 0.45. This analysis asks whether the apparent perturbational-complexity collapse requires nonlinear adaptive dynamics.

4. **Figure-quality correction.** Figure 4 was re-laid out to restore clipped x-axis labels; no underlying numerical result or prospective biological endpoint was changed.

## Frozen hashes for the strengthened computational snapshot

```text
7e9d031fb74cfbd7f25294e424a8662aaee1726f339749803160cf94b4c97fe5  Aneural_Adaptive_Dynamics_Consciousness_Article.docx
b42c76d0390addce54e9861332daa552235b23df097ca0e2a0be40703ccc328b  Supplementary_Code_2.py
af4497d05ed82eea30750b995193e088e91b93bda2a55d65b925454972f5dd22  ExtendedData_Figure3_regime_classifier.png
719eaa8bde2c7a8dbe71fea7b09ee59ffd097bf359e591eb43fa799c0e6d3a00  ExtendedData_Table3_cross_topology_classifier.csv
47fec2ad375fb3e60d5e6cbd8fd4826ef30a534736ff6260fa1a9c63131b8885  ExtendedData_Table3b_classifier_predictions.csv
a6e8bae262db9dcce7f4818af850b061d728714d8e6b2826ba9bb5402ba397bc  ExtendedData_Table4_linear_nulls.csv
5a06f532cf41a8233bbfa052e187e4ce073274bef17ce5f919c6255b98665a5b  ExtendedData_Table4b_linear_null_raw.csv
```

`Supplementary_Code_2.py` passed Python syntax compilation for this snapshot. The numerical outputs listed above are computational/synthetic outputs and are not biological validation.

## Prospective lock remains intact

The original *Stentor* confirmatory model comparison remains the registered biological test. In particular, rejection of the local adaptive-gain null still requires all three preregistered gates to pass, and even such a rejection is not to be interpreted as proof of phenomenal consciousness.
