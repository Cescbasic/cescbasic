# Mtb alveolar macrophage transcriptomic reanalysis

This folder contains a consolidated, reproducible implementation of the final analysis workflow used for the manuscript:

**Temporal inflammatory-to-interferon reprogramming of human alveolar macrophages after Mycobacterium tuberculosis infection: cross-platform validation and splice-junction analysis**

## Public datasets
- GSE189996: Ion Torrent Proton AmpliSeq, 124 deposited libraries, 22 donors. Three libraries were excluded before primary modelling after convergent processed-expression and raw-read QC:
  - SRR17096572 (D12, Control 24 h)
  - SRR17096600 (D17, Mtb 72 h)
  - SRR17096618 (D22, Control 2 h)
- GSE223863: Illumina HiSeq 3000 single-end RNA-seq, 30 libraries, 8 donors.

## Final locked GSE189996 model
- raw reconstruction: STAR 2.7.11b, GRCh37.p13, GENCODE v19
- formal universe: 10,879 genes after edgeR::filterByExpr
- TMM normalisation
- design: ~ donor + time * condition
- voomWithQualityWeights
- limma + robust empirical Bayes
- six prespecified contrasts:
  1. Mtb_vs_Control_2h
  2. Mtb_vs_Control_24h
  3. Mtb_vs_Control_72h
  4. Response_24h_vs_2h
  5. Response_72h_vs_2h
  6. Response_72h_vs_24h
- BH FDR < 0.05

## GSE223863 gene-level model
- STAR/featureCounts, reverse-strand counting
- 16,027 genes after filtering/TMM
- limma-voom with donor blocking using duplicateCorrelation
- consensus within-donor correlation ~0.339
- the same six contrasts as GSE189996

## Pathway analysis
MSigDB 2025.1.Hs Hallmark and Reactome collections.
- GSE189996: limma::camera, inter.gene.cor = 0.01, use.ranks = FALSE
- GSE223863: cameraPR on model statistics

## Sex-interaction model
The GSE189996 sex analysis used the locked 121-sample cohort. Because sex is nested within donor, sex main effects were not estimable with donor fixed effects. The full-rank design added five interaction terms to the primary donor-aware design:
- sexM_time24
- sexM_time72
- sexM_Mtb
- sexM_Mtb_time24
- sexM_Mtb_time72

Positive contrasts indicate a stronger response in males; negative contrasts indicate a stronger response in females.

## Strict splice-junction validation
GSE223863 STAR junctions were filtered to junctions with >=3 reads in >=3 samples and genes with >=2 supported junctions. Junctions were then restricted to exact matches in the GENCODE v50 STAR sjdb list. The final strict matrix contained 105,122 junctions across 9,496 genes. Donor-aware edgeR quasi-likelihood differential-splicing inference was used. Only HLA-DQB2 and RPL13A remained genome-wide significant after final GENCODE ownership and junction-level validation. Salmon-only MX2/STAT1/STAT2 transcript signals were not accepted as validated switches.

## Scripts
- `01_GSE189996_final_analysis.R`: filtering, TMM, primary DE, sensitivity, pathways, sex interaction and sex pathway sensitivity.
- `02_GSE223863_gene_level_analysis.R`: donor-blocked gene-level replication and pathway analysis.
- `03_GSE223863_strict_junction_DTU.R`: strict GENCODE-only junction differential-usage inference from the prepared STAR junction matrix.
- `04_cross_cohort_concordance.R`: gene- and Hallmark-level cross-cohort concordance and manuscript summary tables.

## Reproducibility note
The original CHPC analysis was built interactively in numbered R/PBS stages. These scripts consolidate the final locked statistical specifications, thresholds, exclusions and contrasts used in the completed reanalysis into a clean reproducible form. Raw public FASTQs are not duplicated here.

## Manuscript code-availability URL
https://github.com/Cescbasic/cescbasic/tree/main/Mtb_alveolar_macrophage_reanalysis
