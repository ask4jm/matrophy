# Skeletal Muscle Atrophy Analysis — Python Reproduction Suite

This directory contains the **Python-based computational reproduction** of the landmark *Cell Reports* study:

> **"Integrated genomic and proteomic analyses identify stimulus-dependent molecular changes associated with distinct modes of skeletal muscle atrophy"**  
> *Liam C. Hunt, Flavia A. Graca, Vishwajeeth Pagala, Yong-Dong Wang, Yuxin Li, Zuo-Fei Yuan, Yiping Fan, Myriam Labelle, Junmin Peng, and Fabio Demontis.*  
> **Cell Reports** 37(6): 109971 (2021).  
> **DOI:** [10.1016/j.celrep.2021.109971](https://doi.org/10.1016/j.celrep.2021.109971) | **PMCID:** [PMC8852763](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/) | **GEO:** [GSE159952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159952)

---

## 1. Overview & Stack

The Python reproduction suite implements differential expression analysis, sample-level principal component analysis (PCA), overlap quantification, and publication-quality figure generation using standard scientific Python:
- **NumPy & Pandas**: Matrix manipulations, filtering, and table management
- **SciPy**: Two-tailed Student's t-test (`scipy.stats.ttest_ind`) and Singular Value Decomposition (`scipy.linalg.svd` for PCA)
- **Matplotlib**: 300 DPI publication figure generation

---

## 2. Directory Structure

```
python/
├── README.md                              # This file: Python pipeline documentation
├── src/
│   ├── 00_download_data.py               # Data Acquisition: Automated GEO download & repository guides
│   ├── run_all.py                        # Master Orchestrator: runs full pipeline end-to-end
│   ├── differential_expression.py        # Differential expression pipeline (Student's t-test)
│   ├── plot_reproduction_figures.py      # Publication-quality figure generation (300 DPI)
│   ├── 03_multiomics_integration.py      # RNA-Protein Correlation, Disconnect, Core Atroproteins
│   ├── 05_muscle_mass_correlation.py     # Sarcopenia Phenotype: Proteome vs Muscle Mass Correlation
│   └── extract_captions.py               # Extraction of paper captions and statistical parameters
├── results/                               # Generated tables and statistical summaries
│   ├── Dexamethasone_DEGs.csv            # Dex differential expression results
│   ├── LLC_Cancer_DEGs.csv               # LLC differential expression results
│   ├── Aging_DEGs.csv                    # Aging differential expression results
│   ├── siCYR61_DEGs.csv                  # Cyr61 knockdown differential expression results
│   ├── Figure2D_Overlap_Summary.csv      # Quantified stimulus-specific vs shared DEGs
│   ├── Figure2F_Classic_Atrogenes.csv    # Fold changes & p-values for classic atrogenes
│   ├── Figure6_Wnt_Pathways_siCYR61.csv  # Cyr61/Wnt target modulation
│   ├── Figure7B_Top_Muscle_Mass_Correlating_Proteins.csv # Top mass-correlating proteins
│   ├── FigureS1C_PCA_Coordinates.csv     # Sample PCA scores (PC1, PC2)
│   ├── shared_all_3_upregulated_genes.txt
│   └── shared_all_3_downregulated_genes.txt
└── figures/                               # Generated 300 DPI figures
    ├── Figure2_VolcanoPlots.png          # Figure 2A-C: Volcano plots
    ├── Figure2D_Stimulus_Specific_vs_Shared.png # Figure 2D: Overlap breakdown
    ├── Figure2F_Classic_Atrogenes.png    # Figure 2F: Expression of 16 classic atrogenes
    ├── Figure4A_RNA_Protein_Correlation.png     # Figure 4A: Matched transcriptome-proteome scatter
    ├── Figure4B_Disconnect_Barplot.png          # Figure 4B: mRNA vs Protein disconnect
    ├── Figure5A_Top_Atroproteins.png            # Figure 5A: Consistently regulated atroproteins
    ├── Figure6_Cyr61_Knockdown_Wnt.png   # Figure 6F & 6H: Cyr61/Wnt modulation
    ├── Figure7B_Muscle_Mass_Correlation.png     # Figure 7B: Protein correlation with TA muscle mass
    └── FigureS1C_PCA.png                 # Figure S1C: PCA across 30 samples
```

---

## 3. Raw Data, Processed Matrices & Supplementary Downloads

### A. Automated Python Download
```powershell
python python/src/00_download_data.py
```
This automatically retrieves:
- `data/GSE159952_all.tpm.txt.gz`: Complete 30-sample genome-wide TPM matrix (3.36 MB)
- `data/GSE159952_series_matrix.txt.gz`: GEO platform and experimental metadata

### B. Raw RNA-Seq Reads (FASTQ) Download
- **BioProject**: PRJNA670942 | **SRA Study**: SRP288335
- 30 biological replicates mapped in `data/sample_metadata_sra.csv` (`SRX9350314` – `SRX9350343`).
- Download using SRA Toolkit:
  ```bash
  fasterq-dump --split-files SRX9350314 -O data/fastq/
  ```

### C. Raw Proteomics Mass Spectrometry (TMT LC-MS/MS)
Raw `.raw` and `.pepXML` runs are available via PRIDE:
- `PXD027464`: Aging / Sarcopenia TA muscle
- `PXD027490`: Dexamethasone treated TA muscle
- Download URL: `ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027464/`

### D. Supplementary Spreadsheets (Cell Reports / PMC8852763)
Author-curated supplementary spreadsheets (Table S1 to Table S6):
- Visit: [PMC8852763 Supplementary Section](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/#SM1)
- Place downloaded Excel files into `data/supplementary/`.

---

## 4. How to Run

From the repository root (`matrophy/`):
```powershell
# Run the entire Python reproduction suite end-to-end
python python/src/run_all.py
```

Or execute individual steps:
```powershell
# 1. Download Data
python python/src/00_download_data.py

# 2. Run Differential Expression Pipeline
python python/src/differential_expression.py

# 3. Generate Transcriptomics Reproduction Figures (Figures 2, 2D, 2F, 6, S1C)
python python/src/plot_reproduction_figures.py

# 4. Multi-Omics Integration (Figures 4A, 4B, 5A)
python python/src/03_multiomics_integration.py

# 5. Sarcopenia Phenotype Muscle Mass Correlation (Figure 7B)
python python/src/05_muscle_mass_correlation.py
```

*(Note for Windows environments with QGIS Python: `& "C:\Program Files\QGIS 4.2.3\bin\python-qgis.bat" python/src/run_all.py`).*

---

## 4. Key Quantitative Results

### A. Differential Expression Counts ($p < 0.05$ and $|\log_2\text{FC}| > 1.0$)
- **Dexamethasone**: 1,266 upregulated, 391 downregulated (1,657 total DEGs).
- **LLC Cancer Cachexia**: 175 upregulated, 381 downregulated (556 total DEGs).
- **Aging (24m vs 6m)**: 430 upregulated, 40 downregulated (470 total DEGs).
- **siCYR61 Knockdown**: 96 upregulated, 64 downregulated (160 total DEGs).

### B. Stimulus-Specificity (Figure 2D)
- Over **94.1%** of Dexamethasone, **84.6%** of LLC, and **85.3%** of Aging DEGs are **unique** to their respective conditions.
- Only **1 gene** is universally upregulated (*Sirpb1b*) and **1 gene** universally downregulated (*Ptx4*) across all three modes.

### C. Classic Atrogene Profiles (Figure 2F)
- `Trim63` (MuRF1): Dex $\log_2\text{FC} = +1.05$ ($p < 0.05$), LLC $\log_2\text{FC} = +0.49$ ($p < 0.05$), Aging $\log_2\text{FC} = +0.09$ (ns).
- `Fbxo32` (Atrogin-1): Dex $\log_2\text{FC} = +0.92$ ($p < 0.05$), LLC $\log_2\text{FC} = +0.67$ ($p < 0.05$), Aging $\log_2\text{FC} = +0.34$ (ns).
- `S100a9` and `Il1b`: Universally induced across all 3 conditions ($\log_2\text{FC} \approx +1.5 \text{ to } +2.1$).

