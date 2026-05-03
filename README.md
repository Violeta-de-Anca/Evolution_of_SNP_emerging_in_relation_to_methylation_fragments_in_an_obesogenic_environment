## What is this repository?

This repository contains the scripts used in the analyses presented in the manuscript:

**“Germline genomic and methylomic dynamics following three generations of early-life metabolic challenges”**  
V. de Anca Prado, F. Pértille, D. Andersson, M. Mourin-Fernandez, M. Gòdia, J.C. Jiménez-Chillarón, J. Rüegg, C. Guerrero-Bosagna (corresponding author)

The project investigates how a **multigenerational early-life metabolic challenge** impacts the **sperm genome** and **methylome** across generations, and whether there is any dynamic between **genetic novelties (SNP/CNV)** and **methylomic signals**.

---

## Study design

- **Model:** outbred ICR-CD1 mice
- **Exposure:** litter size reduction during neonatal period  
  - **OG (metabolic challenged group):** 4 pups/litter  
  - **CG (control group):** 8 pups/litter
- **Duration:** repeated across **three generations (F0–F2)**
- **Lineage tracked:** paternal line (males mated to naïve external females each generation)
- **Tissue:** sperm isolated from testes

---

## Data types

This study uses **GBS-MeDIP**, which produces:
- **GBS** libraries (genotyping-by-sequencing, for reduced representation of the genome)
- **GBS-MeDIP** libraries (methylation-enriched fragments in the same reduced representation library)

With these methods we obtained:
- On average with the GBS libraries, we mapped ~ 0.3 % of the mouse genome per individual
- On average with the GBS-MeDIP libraries, we mapped ~ 0.03 % of the mouse genome per individual

The BAM files related to both GBS and GBS-MeDIP libraries are in ENA, in the project with number: PRJEB105468


---

## Main findings

- **Genome:** SNP-based clustering suggested a detectable treatment-related signal (unrelated OG families clustered together in the SNP space).
- **Genome instability:** CNV events were detected in OG, with enrichment in **LINEs** and **LTRs** in several events, consistent or expanding events across generations.
- **Methylome:** PCA of methylation windows did not show a strong global treatment separation, but downstream analyses identified:
  - pathway-level differences (including metabolic/developmental pathways),
  - regional differences involving transposable/repetitive elements,
  - reduced methylation in **LINE** elements in OG across generations in parts of the dataset.
- **Genome–epigenome dynamics:** OG showed reduced emergence of novel SNPs in offspring compared with CG in the analyzed fragments, and the methylation–novel SNP relationship differed between groups.

> Note: This repository focuses on reproducing the computational analyses; interpretation belongs to the manuscript.

---

## Repository layout

At the top level you will find:

- `GBS_violeta/`  
  Scripts for **GBS genetic analyses**, including SNP calling/filters and downstream population-genetic analyses (PCA, hierarchical clustering/IBS, SFS, etc).

- `GBS-MEDIP/`  
  Scripts for **GBS-MeDIP methylation analyses**, including window definition, counting, normalization, differential methylation, annotation, and pathway analyses.


---

## Analysis workflow

```mermaid
flowchart TD
  A[Raw FASTQ] --> B[Demultiplex and trimming with QC]
  B --> C1[GBS alignment]
  B --> C2[GBS-MeDIP alignment]

  C1 --> D1[SNP calling and filtering]
  D1 --> E1[PCA based on SNP, formation of IBS tree, SFS test]

  C1 --> F1[CNV calling]
  F1 --> G1[Repeat enrichment and variant annotation (VEP)]

  C2 --> D2[Window definition (MACS3)]
  D2 --> E2[Count matrix creation (featureCounts) and normalization]
  E2 --> F2[PCA methylome and density profiles]
  E2 --> G2[DMA with intragenerational and intergenerational comparisions]
  E2 --> H2[Annotation of genes, repeat elements enrichment and pathway enrichment]
  E2 --> I[Genome-epigenome dynamics: level of methylation in father vs novel SNP and CNV emergence in offspring]
```
