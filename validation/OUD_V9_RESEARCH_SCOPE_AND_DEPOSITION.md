# OUD–EGR1 research release candidate: evidence boundaries (2026-10-25)

**Status: public reproducibility documentation, not a DOI, institutional approval, or peer-reviewed validation.** This folder holds analysis scripts for public RNA-seq datasets. The executable model runs and downloadable result archives were reviewed separately.

## What is independently supported
- GSE174409: 40 paired donors in DLPFC and NAc; real robust edgeR QL full-cohort DLPFC EGR1 negative log2FC (−0.864), region-wide BH q=0.0269; across 40 separately refitted edgeR leave-one-donor-out analyses 40/40 negative direction, and 37/40 with q<0.05. These are *in-cohort* sensitivity checks, not independent replications.
- GSE182321: published independent BA9 brain series; independent public-count-based edgeR run returned a negative EGR1 effect but not genome-wide 5% FDR. A second software version changed q, so do not pool the q-values as though they were independent datasets.
- Other *independent donor cohorts*, Corradin and Sosnowski, have published EGR1 signed effects in Carter et al. (2026) Supplementary Table 6. This provides **published cross-cohort directional context**, NOT a new independent rerun by these authors or novel discovery. Carter already meta-analyzed GSE174409 and GSE182321.
- GSE198121: seven living male OUD cases and seven controls in blood. Published VST matrix, exploratory Welch p=0.229 and BH q=0.866 for EGR1; negative direction but NOT replication of a brain-tissue effect.
- HLA-DRB1/HLA-DRB5 and MTRNR2L12 exploratory sex×OUD results are model-, donor- or locus-sensitive. Candidate HLA inference must remain exploratory; no personalized HLA alignment/read-level validation was run.

## Independent brain cohort newly identified (NOT analyzed here)
Green et al. (bioRxiv preprint 2026-03-03), human prefrontal cortex BA10 single-nucleus multiomic data, 44 donors across control, OUD, HIV and OUD+HIV cohorts. NeMO Collection nemo:col-c9h2ph2 has restricted access to some data, and raw protected data use dbGaP phs004026.v1.p1. The preprint does not provide a separately audited EGR1 OUD/control human BA10 effect from our analysis. Do not call it completed replication.
Source https://pmc.ncbi.nlm.nih.gov/articles/PMC13001470/

## Phenotype identifiability
If all source OUD cases died from overdose while controls did not, terminal-overdose and OUD case status are perfectly coincident. Their separate causal/associational coefficients are not identifiable from this sample without an independently varying exposure. This cannot be fixed statistically by merely adding regression covariates. Therefore write **postmortem OUD/fatal-overdose case definition** rather than causally isolated chronic OUD.

## Manuscript decisions for revision
- Make DLPFC EGR1 inferential robustness the primary finding.
- Put HLA and sex-interaction candidates into exploratory supplement, with no validated mechanism claim.
- Keep published cross-cohort signed evidence in context, without pretending it was newly re-executed.
- Treat living-blood data as a negative, underpowered, different-tissue follow-up.
- Separately confirm all author contributions, COI, funding, and institutional research determination. AI drafting/code assistance was used and must be disclosed honestly.

## Archiving
The Git commit SHA is a public, version-specific code reference but is **not an archival DOI**. Zenodo archival DOI requires a linked account and successful deposit/release by the repository owner. Do not claim that an archival DOI, Zenodo release or journal submission exists before confirmation.

## Public data and papers
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174409
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182321
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198121
- https://doi.org/10.1038/s41398-026-03977-9
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13001470/

**This document is a research disclosure and reproducibility guide, not institutional sign-off.**
