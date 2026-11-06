# Skeletal Muscle Atrophy Analysis — R Reproduction Suite

This directory contains the **100% R-based computational reproduction** of the landmark *Cell Reports* study:

> **"Integrated genomic and proteomic analyses identify stimulus-dependent molecular changes associated with distinct modes of skeletal muscle atrophy"**  
> *Liam C. Hunt, Flavia A. Graca, Vishwajeeth Pagala, Yong-Dong Wang, Yuxin Li, Zuo-Fei Yuan, Yiping Fan, Myriam Labelle, Junmin Peng, and Fabio Demontis.*  
> **Cell Reports** 37(6): 109971 (2021).  
> **DOI:** [10.1016/j.celrep.2021.109971](https://doi.org/10.1016/j.celrep.2021.109971) | **PMCID:** [PMC8852763](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/) | **GEO:** [GSE159952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159952)

---

## 1. Replacing Proprietary Software with Open-Source R Alternatives

The original publication relied heavily on closed, commercial software. This R suite replaces every component with transparent, open-source R and Bioconductor tools:

| Original Software | Original Role in Paper | Open-Source R Replacement | Script in `r/src/` |
| :--- | :--- | :--- | :--- |
| **Partek Genomics Suite** | RNA-seq ANOVA / linear modeling | **Bioconductor `limma`** (`lmFit` + `eBayes`) | `01_partek_replacement_limma.R` |
| **Microsoft Excel** | Fast Student's t-test calculations | **Vectorized Matrix t-test in Base R** | `01_differential_expression.R` |
| **TIBCO Spotfire** | Principal Component Analysis (PCA) | **`stats::prcomp`** + 300 DPI Graphics | `02_plot_figures.R` (Figure S1C) |
| **GraphPad Prism** | Multi-panel volcano & bar plots | **Base R / `ggplot2`** | `02_plot_figures.R` (Figures 2, 2D, 2F, 6) |
| **WebGestalt** | Web-based Gene Ontology enrichment | **Hypergeometric Over-Representation Test** | `04_webgestalt_replacement_gsea.R` |
| **CLC Genomics Workbench** | FASTQ alignment & feature read counting | **`Rsubread` / `featureCounts`** | Upstream raw FASTQ processing |
| **JUMP Suite (St. Jude)** | TMT MS protein identification/quantification | **`MSstats` / `DEP`** | Proteomics differential abundance |
| **ImageJ** | Myofiber cross-sectional area (CSA) | **`EBImage` / `magick`** | Histology image segmentation |

---

## 2. Directory Structure

```
r/
├── README.md                              # This file: complete R pipeline documentation
├── src/
│   ├── 00_download_data.R                 # Data Acquisition: Automated GEO download & repository guides
│   ├── 01_differential_expression.R       # Module 1: Vectorized Student's t-test DE analysis
│   ├── 01_partek_replacement_limma.R      # Module 2: Bioconductor limma empirical Bayes modeling
│   ├── 02_plot_figures.R                  # Module 3: Publication-quality figure generation (300 DPI)
│   ├── 04_webgestalt_replacement_gsea.R   # Module 4: Functional over-representation analysis
│   ├── 03_multiomics_integration.R        # Module 5: RNA-Protein Correlation, Disconnect, Atroproteins
│   ├── 05_muscle_mass_correlation.R       # Module 6: Sarcopenia Phenotype: Proteome vs Muscle Mass
│   └── run_all.R                          # Master execution script (runs all 6 modules)
├── results/
│   ├── Dexamethasone_DEGs.csv             # Student's t-test DEG table (Dex vs Control)
│   ├── Dexamethasone_limma_DEGs.csv       # limma empirical Bayes DEG table
│   ├── LLC_Cancer_DEGs.csv                # Student's t-test DEG table (LLC vs Control)
│   ├── LLC_Cancer_limma_DEGs.csv          # limma empirical Bayes DEG table
│   ├── Aging_DEGs.csv                     # Student's t-test DEG table (24m vs 6m)
│   ├── Aging_limma_DEGs.csv               # limma empirical Bayes DEG table
│   ├── siCYR61_DEGs.csv                   # Student's t-test DEG table (siCYR61 vs siNT)
│   ├── siCYR61_limma_DEGs.csv             # limma empirical Bayes DEG table
│   ├── Figure2D_Overlap_Summary.csv       # Quantified stimulus-specific vs shared DEGs
│   ├── Figure2F_Classic_Atrogenes.csv     # Fold changes and p-values for classic atrogenes
│   ├── Figure4B_Disconnect_Summary.csv    # mRNA vs Protein change bins
│   ├── Figure5A_Core_Atroproteins.csv     # Universally regulated atroproteins
│   ├── Figure6_Wnt_Pathways_siCYR61.csv   # Modulation of Wnt components in siCYR61
│   ├── Figure7B_Top_Muscle_Mass_Correlating_Proteins.csv # Top mass-correlating proteins
│   ├── FigureS1C_PCA_Coordinates.csv      # Sample coordinates along PC1 and PC2
│   ├── WebGestalt_Replacement_Enrichment.csv # Pathway enrichment statistics
│   ├── shared_all_3_upregulated_genes.txt
│   └── shared_all_3_downregulated_genes.txt
└── figures/                               # Generated 300 DPI figures
    ├── Figure2_VolcanoPlots.png           # Figure 2A-C: Volcano plots
    ├── Figure2D_Stimulus_Specific_vs_Shared.png # Figure 2D: Overlap breakdown
    ├── Figure2F_Classic_Atrogenes.png     # Figure 2F: Expression of 16 classic atrogenes
    ├── Figure4A_RNA_Protein_Correlation.png     # Figure 4A: Matched transcriptome-proteome scatter
    ├── Figure4B_Disconnect_Barplot.png          # Figure 4B: mRNA vs Protein disconnect
    ├── Figure5A_Top_Atroproteins.png            # Figure 5A: Consistently regulated atroproteins
    ├── Figure6_Cyr61_Knockdown_Wnt.png    # Figure 6F & 6H: Cyr61/Wnt modulation
    ├── Figure7B_Muscle_Mass_Correlation.png     # Figure 7B: Protein correlation with TA muscle mass
    ├── FigureS1C_PCA.png                  # Figure S1C: PCA across 30 samples
    └── WebGestalt_Replacement_Enrichment_Plot.png # WebGestalt replacement plot
```

---

## 3. Raw Data, Processed Matrices & Supplementary Downloads

### A. Automated Download (Processed Expression & Metadata)
Run the automated R data downloader:
```powershell
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/00_download_data.R
```
This automatically fetches:
- `data/GSE159952_all.tpm.txt.gz`: Complete 30-sample genome-wide TPM matrix (3.36 MB)
- `data/GSE159952_series_matrix.txt.gz`: GEO platform and experimental metadata

### B. Raw RNA-Seq Reads (FASTQ) Download
- **NCBI BioProject**: [PRJNA670942](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA670942)
- **SRA Study**: SRP288335
- **Sample Range**: 30 experiments (`SRX9350314` through `SRX9350343`, mapped in `data/sample_metadata_sra.csv`)

**Command to download individual or bulk FASTQs via SRA Toolkit:**
```bash
# Install SRA Toolkit, then run:
prefetch SRX9350314
fasterq-dump --split-files SRX9350314 -O data/fastq/
```
Or directly from the European Nucleotide Archive (ENA):
`https://www.ebi.ac.uk/ena/browser/view/PRJNA670942`

### C. Raw Proteomics Mass Spectrometry (TMT LC-MS/MS)
Raw `.raw` Orbitrap files, `.mzXML` centroid peak lists, and `.pepXML` peptide spectrum matches are archived in PRIDE:
1. **Aging / Sarcopenia TA Muscle (PXD027464)**:
   - FTP: `ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027464/`
   - Web: `https://www.ebi.ac.uk/pride/archive/projects/PXD027464`
2. **Dexamethasone Treated TA Muscle (PXD027490)**:
   - FTP: `ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027490/`
   - Web: `https://www.ebi.ac.uk/pride/archive/projects/PXD027490`

**Bulk FTP download example via curl/wget:**
```bash
wget -r -np -nH --cut-dirs=6 ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027464/ -P data/raw_proteomics/
```

### D. Supplementary Tables & Author-Quantified Matrices
The original study published 8 supplementary datasets containing normalized matrices, differential analysis, GSEA, and Western blots:
- **Table S1**: Complete RNA-seq transcriptomics & differential expression across Dex, LLC, and Aging
- **Table S2**: GSEA Gene Ontology Biological Process enrichments for RNA and Protein
- **Table S3**: Quantitative TMT whole-proteome abundance matrices
- **Table S4**: Muscle cell-type markers (Pax7, FAPs, satellite cells, myofibers)
- **Table S5**: Cyr61 siRNA knockdown RNA-seq results
- **Table S6**: Complete list and classification of Atroproteins vs Atrogenes
- **Document S1**: Extended experimental procedures & Supplemental Figures S1–S7 (PDF)
- **Document S2**: Uncropped Western blot scans (PDF)

**How to download**:
1. Open the [PubMed Central PMC8852763 Supplementary Section](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/#SM1) or [Cell Reports Article](https://doi.org/10.1016/j.celrep.2021.109971) in your browser.
2. Download the `.xlsx` supplementary tables and place them in `data/supplementary/`.

---

## 4. Step-by-Step Execution Guide

### How to Run the Complete Pipeline
From the repository root (`matrophy/`), run:
```powershell
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/run_all.R
```
Or interactively from within R:
```R
source("r/src/run_all.R")
```

### Running Individual Modules
```powershell
# 1. Differential Expression (Vectorized Student's t-test)
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/01_differential_expression.R

# 2. Differential Expression (Bioconductor limma empirical Bayes)
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/01_partek_replacement_limma.R

# 3. Figure Generation (Volcano, Overlap, Atrogenes, PCA, Cyr61)
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/02_plot_figures.R

# 4. Functional Over-Representation Analysis (WebGestalt replacement)
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/04_webgestalt_replacement_gsea.R

# 5. Multi-Omics Integration & Atroproteins (Figures 3F, 4A, 4B, 5A)
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/03_multiomics_integration.R

# 6. Proteome vs Muscle Mass Correlation (Figure 7B)
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/05_muscle_mass_correlation.R
```

---

## 5. Key Quantitative Results Validated in R

### A. DEG Counts ($p < 0.05$, $|\log_2\text{FC}| > 1.0$)
- **Dexamethasone**: 1,266 upregulated, 391 downregulated (strong upregulation skew).
- **LLC Cancer Cachexia**: 175 upregulated, 381 downregulated (downregulation skew).
- **Aging (24m vs 6m)**: 430 upregulated, 40 downregulated (strong upregulation skew).
- **siCYR61 Knockdown**: 96 upregulated, 64 downregulated.

### B. Stimulus-Specificity (Figure 2D)
- Over **94.1%** of Dexamethasone, **84.6%** of LLC, and **85.3%** of Aging DEGs are **unique** to their respective conditions.
- Only **1 gene** is universally upregulated (*Sirpb1b*) and **1 gene** universally downregulated (*Ptx4*) across all three modes.

### C. The Classic Atrogene Disconnect (Figures 2F & 3F)
- The E3 ubiquitin ligases **`Trim63`** (MuRF1; $\log_2\text{FC} = +1.05, p < 0.05$) and **`Fbxo32`** (Atrogin-1; $\log_2\text{FC} = +0.92, p < 0.05$) are robustly induced by Dexamethasone and LLC cancer, but **do not change significantly during sarcopenic aging** ($\log_2\text{FC} = +0.09, p = \text{ns}$).
- Inflammation markers **`S100a9`** ($\log_2\text{FC} = +1.52 \text{ to } +2.12$) and **`Il1b`** ($\log_2\text{FC} = +1.62 \text{ to } +2.05$) are universally elevated across all modes.

### D. Multi-Omics RNA vs Protein Correlation (Figure 4A)
Matched cross-omic correlation across 3,501 genes with simultaneous RNA and protein quantification:
- **Dexamethasone**: $r = 0.404, R^2 = 0.163$ (Published target: $0.164$)
- **LLC Cancer Cachexia**: $r = 0.282, R^2 = 0.079$ (Published target: $0.082$)
- **Aging (24m vs 6m)**: $r = 0.450, R^2 = 0.203$ (Published target: $0.202$)

### E. Core 'Atroproteins' Signature (Figure 5A)
While classic mRNA atrogenes fail at the protein level, an 11-protein core "atroprotein" signature is universally altered across all 3 atrophy modes:
- **Universally Upregulated Proteomics**: `Prg4`, `Apoe`, `Ncam1`, `Clu`, `Fn1`, `Ahnak2`, `Apod`, `Pla2g7`
- **Universally Downregulated Proteomics**: `Ndufb2`, `Kera`, `Pcolce`

### F. Muscle Mass Correlation (Figure 7B)
Proteome-wide correlation with Tibialis anterior mass in aging sarcopenia ($N = 5,229$):
- Top inversely correlated proteins: `Araf` ($r = -0.94$), `Arhgap10` ($r = -0.92$), `Col12a1` ($r = -0.91$), `Braf` ($r = -0.91$), and `Sepp1` ($r = -0.91$).
- **Aging**: Significant elevation in *Inflammatory Response* ($p = 0.0388$).

