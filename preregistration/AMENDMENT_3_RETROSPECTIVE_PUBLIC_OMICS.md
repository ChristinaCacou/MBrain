# Amendment 3: retrospective public Stentor omics empirical extension

**Date:** 17 October 2026  
**Original prospective preregistration:** commit `07949d983493186d63d3232a546247476adae4bf`  
**Prospective protocol version:** 1.0, unchanged

## Purpose and registration boundary

This amendment records a retrospective analysis of an independently generated, publicly available *Stentor coeruleus* learning/forgetting transcriptomic dataset that was conceived and executed after the prospective perturbation protocol had already been publicly registered.

The analysis is an empirical extension of the computational/theoretical manuscript only. It does **not** modify the prospective sample plan, treatment arms, acquisition sequence, exclusions, endpoints, thresholds, randomization/blinding plan, analysis code, or three conjunctive confirmatory gates. It does not count toward any prospective gate and cannot be used to rescue a failed future prospective experiment.

No new biological observations were collected for this retrospective extension.

## Public biological source

- Rajan DH et al., *Molecular pathways for learning in the single-cell Stentor coeruleus*, Current Biology (2026).
- DOI: `10.1016/j.cub.2026.03.080`
- Public data DOI: `10.5061/dryad.z34tmpgrm`
- Public processed-analysis repository: `aralbright/2025_stentor_learning`

The reanalysis used three biological replicates each for untrained 19 h, trained 19 h, trained + 10 min recovery, and trained + 90 min recovery. The processed matrix contained 34,494 transcripts; 28,296 passed a conservative mean-expression/nonzero-variance filter.

## Retrospective analyses and results

### 1. Transcriptomic learning axis

The mean untrained centroid was defined as 0 and the mean trained centroid as 1. Recovery samples were projected onto that fixed axis.

- 10-min recovery index: **0.4191**, biological-replicate bootstrap 95% CI **0.2051–0.7360**.
- 90-min recovery index: **0.4503**, 95% CI **0.2140–0.7791**.
- 10-min minus 90-min recovery: bootstrap 95% CI **−0.3437 to 0.2386**; the available data therefore do not resolve monotonic decay between those sampled recovery times.

Leave-one-replicate-pair-out mean scores were 0.3535 for untrained, 0.6840 for trained, 0.4697 for 10-min recovery and 0.4533 for 90-min recovery.

### 2. Independent gene-level alignment test

Across 27,229 genes shared by the trained-versus-untrained and trained-versus-90-min-recovery contrasts:

- cosine alignment = **0.5948**;
- gene-identity permutation **P = 0.000999** (1,000 permutations);
- Spearman **ρ = 0.4503**;
- 83 genes met training-contrast q < 0.05;
- **78/83 (93.98%)** retained the same effect direction in the trained-versus-90-min-recovery contrast, one-sided binomial **P = 3.20 × 10⁻18**;
- median `|trained − 90-min recovery| / |trained − untrained|` = **0.5787**.

The training-associated signal was distributed rather than dominated by one transcript: the squared-effect participation ratio corresponded to an effective **35.44 genes**; 14 genes accounted for 50% and 32 for 80% of total squared training-effect energy among the 83 FDR-significant transcripts.

## Interpretation and claim boundary

These real biological data show that *Stentor* learning is accompanied by a distributed, history-dependent and partially reversible molecular state. That strengthens the manuscript's central specificity argument: distributed or globally coordinated state changes can be genuine properties of aneural cellular learning while remaining mechanistically compatible with ordinary molecular and cellular processes. Therefore such signatures are not, by themselves, specific evidence for phenomenal consciousness.

This retrospective analysis does **not** establish consciousness, validate the recurrent-network simulations as a literal biological mechanism, or substitute for the prospectively registered perturbation experiment.

## Frozen strengthened-manuscript snapshot

```text
e7e3505a9a3bcda177945db288cba97c1185703930fd1e3a9fed0bda5d2224bb  Aneural_Adaptive_Dynamics_Consciousness_Article.docx
6d1a26ea58126cab982aba6768d68873ad944bd84952ae5404110328a999b1ab  Supplementary_Code_3_Public_Stentor_Omics.py
3514c5fc49c0816b4b7f45ba7d24cd245a1ec54531c3d084e8d6f3138d281fff  Figure5_public_Stentor_transcriptomics.png
fe35295d2c2e3db2ff6962f6bcec84abcfef2868c02aa7f7d25c5e2828232729  ExtendedData_Table5_public_Stentor_omics.csv
68bc0d16768622a6dc99f073fde95c17600532865b8acadab8232836a24a9de9  Retrospective_Stentor_Omics_Results.json
49872d61a94a6fce82a864a9c2e2afeca6699f317baf167daa335be5c1097932  Figure5a_memory_axis.csv
e1119666396ec444bfa3747bf1456c48b3cb24a9029c19ff442c0acfebf92698  Figure5b_crossfitted_memory_axis.csv
f601d8b4c54e2568f0d639814d06086260b50427244fade6c4009e64c1af4a1b  Figure5c_alignment_quantile_bins.csv
83757d4d3c7063948fe4880dd77ce0c16f07e18c5ed8a90dc766e1c8a118caa5  Figure5d_effect_energy_curve.csv
```

The prospective *Stentor* confirmatory protocol remains unchanged and prospectively locked.