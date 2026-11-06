# Reproduction: Integrated Genomic & Proteomic Analyses of Skeletal Muscle Atrophy

This repository contains the computational reproduction suite for the *Cell Reports* study:

> **"Integrated genomic and proteomic analyses identify stimulus-dependent molecular changes associated with distinct modes of skeletal muscle atrophy"**  
> *Liam C. Hunt, Flavia A. Graca, Vishwajeeth Pagala, Yong-Dong Wang, Yuxin Li, Zuo-Fei Yuan, Yiping Fan, Myriam Labelle, Junmin Peng, and Fabio Demontis.*  
> **Cell Reports** 37(6): 109971 (November 9, 2021).  
> **DOI:** [10.1016/j.celrep.2021.109971](https://doi.org/10.1016/j.celrep.2021.109971) | **PMCID:** [PMC8852763](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/) | **GEO:** [GSE159952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159952) | **PRIDE:** PXD027464 & PXD027490

---

## 1. Dedicated Reproduction Suites

The codebase is partitioned into two dedicated language-specific reproduction suites, each with its own comprehensive documentation:

- 📁 [**`r/README.md` — Complete R Pipeline**](r/README.md):
  - **100% R Implementation**: Replaces proprietary software (Partek, Spotfire, GraphPad Prism, WebGestalt) with open-source R & Bioconductor packages (`limma`, `stats::prcomp`, `ggplot2`, and hypergeometric over-representation testing).
  - High-speed vectorized statistical testing and publication-quality 300 DPI visualization.
- 📁 [**`python/README.md` — Complete Python Pipeline**](python/README.md):
  - NumPy, SciPy, Pandas, and Matplotlib implementation for independent cross-language validation.

```
matrophy/
├── README.md                              # This file: Main repository landing guide
├── test.R                                 # Quick-start script in repository root
├── data/                                  # Common RNA-seq datasets (GSE159952)
│   ├── GSE159952_all.tpm.txt.gz          # Complete TPM expression matrix (47,643 genes, 30 samples)
│   └── GSE159952_series_matrix.txt.gz    # GEO experiment metadata
│
├── r/                                     # Dedicated R Pipeline (See r/README.md)
│   ├── README.md                         # Detailed R documentation
│   ├── src/                              # R scripts (limma, t-test, plotting, GSEA)
│   ├── results/                          # DEG tables, PCA coordinates, enrichment results
│   └── figures/                          # 300 DPI figures reproduced in R
│
└── python/                                # Dedicated Python Pipeline (See python/README.md)
    ├── README.md                         # Detailed Python documentation
    ├── src/                              # Python scripts (differential expression & plotting)
    ├── results/                          # Python DEG tables and overlap summaries
    └── figures/                          # 300 DPI figures reproduced in Python
```

---

## 2. Software Translation: Replacing Closed Software with Open R Alternatives

A primary contribution of this reproduction is replacing the proprietary tools used by Hunt et al. with open-source R packages:

| Original Study Software | Commercial / Closed Role | Open-Source R Replacement | Implemented In |
| :--- | :--- | :--- | :--- |
| **Partek Genomics Suite** | Linear modeling & RNA-seq ANOVA/t-test | **Bioconductor `limma`** (`lmFit` + `eBayes`) | [`r/src/01_partek_replacement_limma.R`](r/src/01_partek_replacement_limma.R) |
| **Microsoft Excel** | Spreadsheet data management & t-tests | **Vectorized Base R** | [`r/src/01_differential_expression.R`](r/src/01_differential_expression.R) |
| **TIBCO Spotfire** | Principal Component Analysis (PCA) | **`stats::prcomp`** + 300 DPI Graphics | [`r/src/02_plot_figures.R`](r/src/02_plot_figures.R) (Figure S1C) |
| **GraphPad Prism** | Statistical plotting, bar charts & ANOVA | **Base R / `ggplot2`** | [`r/src/02_plot_figures.R`](r/src/02_plot_figures.R) (Figures 2, 2D, 2F, 6) |
| **WebGestalt** | Web-based Gene Ontology enrichment | **Hypergeometric Over-Representation Test** | [`r/src/04_webgestalt_replacement_gsea.R`](r/src/04_webgestalt_replacement_gsea.R) |
| **CLC Genomics Workbench** | Read mapping & feature read counting | **`Rsubread` / `featureCounts`** | Upstream raw FASTQ alignment |
| **JUMP Suite (St. Jude)** | TMT MS protein identification/quantification | **`MSstats` / `DEP`** | Proteomics differential abundance |
| **ImageJ** | Muscle cross-section myofiber CSA analysis | **`EBImage` / `magick`** | Automated myofiber segmentation |

---

## 3. Scientific Discoveries Validated Across Both Suites

### A. Stimulus-Specific Molecular Programs (>85–94% Unique)
Hunt et al. disproved the paradigm that diverse catabolic stressors share a common transcriptional program. Our reproduction confirms that the vast majority of differentially expressed genes are strictly stimulus-specific:

| Condition | Upregulated ($p<0.05, \log_2\text{FC}>1$) | Downregulated ($p<0.05, \log_2\text{FC}<-1$) | Total DEGs | % Stimulus-Specific |
| :--- | :---: | :---: | :---: | :---: |
| **Dexamethasone** | 1,266 | 391 | 1,657 | **94.1% Unique** |
| **LLC Cancer Cachexia** | 175 | 381 | 556 | **84.6% Unique** |
| **Aging (24m vs 6m)** | 430 | 40 | 470 | **85.3% Unique** |

Only **1 gene** is universally upregulated (*Sirpb1b*) and **1 gene** universally downregulated (*Ptx4*) across all three atrophic stressors.

### B. The Classic Atrogene Disconnect (Figure 2F)
- The widely cited E3 ubiquitin ligases **`Trim63`** (MuRF1) and **`Fbxo32`** (Atrogin-1) are induced in acute atrophic stressors (Dexamethasone and LLC Cancer), but **do not change during sarcopenic aging** ($\log_2\text{FC} = +0.09$, $p = \text{ns}$).
- Inflammation markers **`S100a9`** and **`Il1b`** are universally elevated across all three modes ($\log_2\text{FC} \approx +1.5 \text{ to } +2.1$).

### C. Cyr61 / Wnt Pathway Modulation (Figure 6)
- Knockdown of *Cyr61/CCN1* in aged muscle alters Wnt pathway components (*Wnt6*, *Wnt7b*, *Sfrp5*, *Rspo1*, *Kremen2*), explaining myofiber type switching.

### D. Multi-Omics Transcriptome-Proteome Disconnect (Figures 3, 4, 5)
- **Low Global Correlation**: RNA fold change explains only **8% to 20%** of protein changes ($R^2_{\text{Dex}} = 0.163$, $R^2_{\text{LLC}} = 0.079$, $R^2_{\text{Aging}} = 0.203$), precisely replicating the study's published $R^2$ values ($0.164, 0.082, 0.202$).
- **Core Atroproteins vs Atrogenes**: While classic atrogenes (*Fbxo32*, *Trim63*) fail to change at the protein level in aging, a true core panel of **atroproteins** (`Prg4`, `Apoe`, `Ncam1`, `Clu`, `Fn1`, `Ahnak2`, `Apod`, `Pla2g7` upregulated; `Ndufb2`, `Kera`, `Pcolce` downregulated) is conserved across all 3 conditions.
- **Muscle Mass Correlates (Figure 7B)**: Top proteins inversely correlated with muscle mass identified in aging sarcopenia (`Araf`, `Arhgap10`, `Col12a1`, `Braf`, `Sepp1`, all $r \le -0.91$).

---

## 4. Quick-Start Execution

### Run the Entire R Suite
```powershell
# Master pipeline via Rscript (runs all 6 modules)
& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/run_all.R

# Or from repository root in an interactive R console:
source("test.R")
```
*(Detailed options, module breakdowns, and limma details: see [**`r/README.md`**](r/README.md)).*

### Run the Entire Python Suite
```powershell
# Master pipeline via Python (runs all 4 modules)
python python/src/run_all.py
```
*(Detailed options and standalone scripts: see [**`python/README.md`**](python/README.md)).*
---

## 5. Raw Data, Processed Matrices & Supplementary Downloads

### A. Automated Download (Processed TPM Matrix & Metadata)
- **Via R**: `& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" r/src/00_download_data.R`
- **Via Python**: `python python/src/00_download_data.py`
- Downloads `data/GSE159952_all.tpm.txt.gz` (3.36 MB, 47,643 genes across 30 samples) and `data/GSE159952_series_matrix.txt.gz`.

### B. Raw RNA-Seq Reads (FASTQ) Download
- **NCBI BioProject**: [PRJNA670942](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA670942) | **SRA Study**: SRP288335
- 30 biological replicates mapped in [`data/sample_metadata_sra.csv`](data/sample_metadata_sra.csv) (`SRX9350314` – `SRX9350343`).
- Download using SRA Toolkit:
  ```bash
  prefetch SRX9350314
  fasterq-dump --split-files SRX9350314 -O data/fastq/
  ```

### C. Raw Proteomics Mass Spectrometry (TMT LC-MS/MS)
Raw `.raw` and `.pepXML` runs are archived in PRIDE:
- **Aging / Sarcopenia TA Muscle**: `PXD027464` ([FTP Download](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027464/))
- **Dexamethasone Treated TA Muscle**: `PXD027490` ([FTP Download](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027490/))

### D. Supplementary Tables & Author-Quantified Matrices
The study published 8 comprehensive supplementary datasets:
- **Table S1**: Complete RNA-seq transcriptomics & differential expression across Dex, LLC, and Aging
- **Table S2**: GSEA Gene Ontology Biological Process enrichments for RNA and Protein
- **Table S3**: Quantitative TMT whole-proteome abundance matrices
- **Table S4**: Muscle cell-type markers (Pax7, FAPs, satellite cells, myofibers)
- **Table S5**: Cyr61 siRNA knockdown RNA-seq results
- **Table S6**: Complete list and classification of Atroproteins vs Atrogenes
- **Document S1**: Extended experimental procedures & Supplemental Figures S1–S7 (PDF)
- **Document S2**: Uncropped Western blot scans (PDF)

Download link: [PubMed Central PMC8852763 Supplementary Section](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/#SM1) (place files in `data/supplementary/`).

