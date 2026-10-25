# Opioid-use-disorder EGR1 robustness: public analysis code and provenance
Status: author-review research tools, 25 October 2026. Not peer-reviewed and not an archival DOI.

## Reproducible computations
- Original GSE174409 paired DLPFC/NAc source: NCBI GEO GSE174409; analysis scripts: validation/run_discovery_edgeR.R, validation/donor22_edgeR_QL.R, validation/edgeR_EGR1_all_donors.R. The CI workflows in .github/workflows/oud-edger-modern-20261025.yml and oud-egr1-full-LOO.yml document versions and inputs. Reported DLPFC EGR1 log2FC negative 40/40, BH q<0.05 in 37/40 donor-deletion fits. This is sensitivity on one dataset, not 40 replications.
- Independent postmortem BA9 GSE182321: original GEO RNA-seq; different R/edgeR versions yield distinguishable p/q values. Does not meet genome-wide replication q<0.05.
- Independently recruited living male clinical blood cohort GSE198121: validation/GSE198121_independent_blood_EGR1.py with .github/workflows/oud-living-blood-gse198121.yml. The source is the GEO-published variance-stabilized normalized 7 OUD + 7 control matrix (no raw count DE). Exploratory EGR1 difference -0.940288751 VST units; Welch p=0.228976257; BH q=0.866165399 across 55133 finite rows. Cross-tissue and non-significant; cannot replicate cortical effects.
- Reference HLA locus diagnostic: .github/workflows/oud-hla-reference-mappability.yml. Ensembl canonical reference cDNAs HLA-DRB1 ENST00000360004 and HLA-DRB5 ENST00000374975. Exact identical unique kmer counts at 35,50,75,100,125,150 bases: 198,102,28,0,0,0. This is NOT patient-allele-aware FASTQ mapping and does not adjudicate assignment quality in original reads.

## Structural nonidentifiability of original OUD/fatal-overdose group
The matched 40-donor GEO metadata used by the reanalysis includes diagnosis, sex, age, PMI and RIN, but does NOT include an independent actual per-donor terminal-overdose variable. Under the published/discussed scenario in which terminal overdose T equals case diagnosis D, any regression including both terms has perfect collinearity: beta_D and beta_T cannot be estimated separately; only their sum is identified. Do not use D as though it were two independently measured exposures. A living peripheral-blood analysis does not identify chronic-brain-specific effects.

## Integrity and journal limitations
- The full results archive and original public source data are made available to the investigator as a separate downloadable research ZIP. These tools are PUBLIC code, not a journal-accepted deposition or immutable Zenodo DOI.
- CI artifacts expire after short retention; do not cite temporary artifact links as permanent data archives.
- HLA genotypes, original read-level alignments, independent prospective brain phenotypes, causality and terminal-overdose separation have not been established.
- Final author approval, CRediT contributions, funding, competing-interest declarations, permissions, institution-specific ethics research determination and any journal submission must be confirmed by human authors and competent authorities.
- AI assistance with drafting and code should be accurately disclosed; do not claim no AI was used.

## Sources and journal
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174409 ;
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182321 ;
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198121 .
Dai et al. (2022) doi:10.3389/fpsyt.2022.893303.
Translational Psychiatry requirements: https://www.nature.com/tp/authors-and-referees/preparation-of-articles .
