# Independent 2023 human ventral-midbrain EGR1 audit — 28 October 2026

**Status:** completed *source-table extraction* from a previously published independent human brain cohort; **not** a new donor-level differential-expression fit or a validated DLPFC replication.

## Source and exact rows
Wei J, Lambert TY, Valada A, et al. *Nature Communications* 2023;14:5610. doi:10.1038/s41467-023-41455-8; GEO GSE240457. This cohort studied 45 opioid-overdose decedents and 50 drug-free controls in ventral midbrain.

Original journal Supplementary Data 3: `41467_2023_41455_MOESM6_ESM.xlsx`. Downloaded and validated as a real workbook by the executed workflow at `.github/workflows/oud-wei-egr1-independent-brain.yml`. SHA-256: `3f29793b1338bdf97927470e3be5143dcf3b1a933032adea463a7c7bc48c6053`.

The original row-2 published labels are: `gene, baseMean, estimate, stderror, statistic, p.value, p.adjusted, MCls, contrast, celltype`.

| Original Excel row | Contrast | Cell type | EGR1 estimate | Std error | P value | FDR |
|---|---|---|---:|---:|---:|---:|
| 463 | Opioid | ODC (oligodendrocytes) | +1.03015565 | 0.369698 | 0.00532841 | 0.05505898 |
| 5551 | Opioid | OPC (oligodendrocyte precursors) | +1.12736382 | 0.459083 | 0.01406154 | 0.08912135 |

Both were positive in ventral-midbrain glial cell types, **opposite direction to the negative DLPFC bulk-tissue EGR1** from GSE174409. Neither satisfies FDR<0.05, but both satisfy the published study's exploratory FDR<0.10 standard. This does not directly contradict the DLPFC result in the same tissue; it shows that the EGR1 direction should not be generalized to every brain region or cell type.

The original authors' published inferential estimates are copied **verbatim to numerical precision**, not re-estimated with raw donor-level count data. Therefore treat them as independent donor context, not independently executed differential-expression validation by the current authors.

## BA10 source access
Green et al. 2026 bioRxiv (doi:10.64898/2026.03.01.708876) reports 7 controls vs 4 OUD-only donors in its BA10 prefrontal OUD comparison. NeMO collection nemo:col-c9h2ph2 currently labels access restricted without a public downloadable donor-level count matrix in the collection metadata; raw data access must follow dbGaP phs004026.v1.p1. We did not perform BA10 pseudobulk analysis.

## Complete public HLA read inputs (not remapped)
GSE174409 SRA project SRP319708: 80 paired-end sample libraries from 40 paired DLPFC/NAc donors; 160 compressed FASTQ mate files total **359,472,774,795 bytes** (359.473 GB decimal) before BAM intermediates or reference. A complete cohort-wide patient-personalized DRB1/DRB5 alignment / re-count / edgeR re-analysis was **not** executed. The 2-library exact 31-mer prefixes from earlier v10 are a preliminary ambiguity screen, not genotyping or full remapping. Validated donor HLA genotypes or explicitly qualified RNA-derived typing, computing resources, and complete alignment QC are still needed.

Full publication package v11 includes a source workbook, 80-row accession/byte manifest, machine-readable audit, and analysis status documentation; no protected human genotype or restricted NeMO files are included.
