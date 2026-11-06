# ==============================================================================
# Script: 00_download_data.R
# Purpose: Automated data downloader for raw and processed datasets
# Paper: Hunt et al. (Cell Reports 2021)
# ==============================================================================

cat("=================================================================\n")
cat("Data Acquisition & Repository Download Script (Hunt et al. 2021)\n")
cat("=================================================================\n\n")

# Determine repository root
repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
data_dir <- file.path(repo_root, "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. DOWNLOAD PROCESSED RNA-SEQ GENOME-WIDE TPM MATRIX (GEO: GSE159952)
# ------------------------------------------------------------------------------
geo_base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159952"

tpm_file <- file.path(data_dir, "GSE159952_all.tpm.txt.gz")
tpm_url  <- paste0(geo_base, "/suppl/GSE159952_all.tpm.txt.gz")

if (!file.exists(tpm_file)) {
  cat(sprintf("Downloading processed TPM matrix from GEO:\n  %s\n", tpm_url))
  download.file(tpm_url, destfile = tpm_file, mode = "wb", quiet = FALSE)
  cat("Download complete!\n\n")
} else {
  cat(sprintf("GSE159952_all.tpm.txt.gz already exists (%0.2f MB).\n\n", file.size(tpm_file) / 1e6))
}

matrix_file <- file.path(data_dir, "GSE159952_series_matrix.txt.gz")
matrix_url  <- paste0(geo_base, "/matrix/GSE159952_series_matrix.txt.gz")

if (!file.exists(matrix_file)) {
  cat(sprintf("Downloading GEO Series Matrix metadata:\n  %s\n", matrix_url))
  download.file(matrix_url, destfile = matrix_file, mode = "wb", quiet = FALSE)
  cat("Download complete!\n\n")
} else {
  cat(sprintf("GSE159952_series_matrix.txt.gz already exists (%0.2f KB).\n\n", file.size(matrix_file) / 1e3))
}

# ------------------------------------------------------------------------------
# 2. INSTRUCTIONS FOR DOWNLOADING RAW RNA-SEQ READS (NCBI SRA / ENA)
# ------------------------------------------------------------------------------
cat("=================================================================\n")
cat("RAW RNA-SEQ READS (FASTQ) DOWNLOAD INSTRUCTIONS\n")
cat("BioProject: PRJNA670942 | SRA Study: SRP288335\n")
cat("=================================================================\n")
cat("To download raw FASTQ sequence files for all 30 samples:\n\n")
cat("Option A: Using SRA-Toolkit (fasterq-dump / prefetch):\n")
cat("  1. Install SRA Toolkit: https://github.com/ncbi/sra-tools\n")
cat("  2. Download by SRX experiment (SRX9350314 to SRX9350343):\n")
cat("     prefetch SRX9350314\n")
cat("     fasterq-dump --split-files SRX9350314 -O data/fastq/\n\n")

cat("Option B: Direct download via European Nucleotide Archive (ENA):\n")
cat("  https://www.ebi.ac.uk/ena/browser/view/PRJNA670942\n\n")

# ------------------------------------------------------------------------------
# 3. INSTRUCTIONS FOR DOWNLOADING RAW PROTEOMICS DATA (PRIDE)
# ------------------------------------------------------------------------------
cat("=================================================================\n")
cat("RAW MASS SPECTROMETRY PROTEOMICS (TMT LC-MS/MS) DOWNLOAD\n")
cat("ProteomeXchange / PRIDE Accessions: PXD027464 & PXD027490\n")
cat("=================================================================\n")
cat("Raw mass spectrometry files (.raw), peak lists (.mzXML), and peptide\n")
cat("identifications (.pepXML) are hosted on EMBL-EBI PRIDE Archive:\n\n")
cat("1. Aging / Sarcopenia Muscle (PXD027464):\n")
cat("   FTP:   ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027464/\n")
cat("   HTTP:  https://www.ebi.ac.uk/pride/archive/projects/PXD027464\n\n")
cat("2. Dexamethasone Treated Muscle (PXD027490):\n")
cat("   FTP:   ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027490/\n")
cat("   HTTP:  https://www.ebi.ac.uk/pride/archive/projects/PXD027490\n\n")
cat("Command-line download example (via curl or wget):\n")
cat("  wget -r -np -nH --cut-dirs=6 ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027464/ -P data/raw_proteomics/\n\n")

# ------------------------------------------------------------------------------
# 4. SUPPLEMENTARY TABLES (CELL REPORTS / PMC)
# ------------------------------------------------------------------------------
cat("=================================================================\n")
cat("SUPPLEMENTARY TABLES & DATASETS (Cell Reports / PMC8852763)\n")
cat("=================================================================\n")
cat("The publication provides author-quantified supplementary spreadsheets:\n")
cat("  - Table S1: Complete RNA-seq differential expression results across Dex, LLC, Aging\n")
cat("  - Table S2: Gene Set Enrichment Analysis (GSEA) results for RNA and Protein\n")
cat("  - Table S3: Quantitative TMT mass spectrometry proteomics results\n")
cat("  - Table S4: Muscle single-cell & myofiber marker abundance\n")
cat("  - Table S5: Cyr61 siRNA knockdown RNA-seq results\n")
cat("  - Table S6: Categorization of Atroproteins & Atrogenes\n")
cat("  - Supplemental Document 1: Extended methods & Supplemental Figures S1-S7 (PDF)\n\n")
cat("To download these files:\n")
cat("  Open in web browser: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/#SM1\n")
cat("  Or via Cell Reports: https://doi.org/10.1016/j.celrep.2021.109971\n")
cat("  Save the downloaded .xlsx files into: data/supplementary/\n")
cat("=================================================================\n")

