"""
Script: 00_download_data.py
Purpose: Automated dataset acquisition for raw, processed, and supplementary data
Paper: Hunt et al. (Cell Reports 2021)
"""

import os
import urllib.request
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
data_dir = os.path.join(repo_root, 'data')
os.makedirs(data_dir, exist_ok=True)

print("=" * 65)
print("Data Acquisition & Repository Download Script (Hunt et al. 2021)")
print("=" * 65)

# 1. Download Processed GEO Files
geo_base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159952"
files_to_download = {
    "GSE159952_all.tpm.txt.gz": f"{geo_base}/suppl/GSE159952_all.tpm.txt.gz",
    "GSE159952_series_matrix.txt.gz": f"{geo_base}/matrix/GSE159952_series_matrix.txt.gz"
}

for fname, url in files_to_download.items():
    dest = os.path.join(data_dir, fname)
    if not os.path.exists(dest):
        print(f"Downloading {fname} from {url}...")
        urllib.request.urlretrieve(url, dest)
        print(f"Saved to {dest}")
    else:
        size_mb = os.path.getsize(dest) / 1e6
        print(f"{fname} already exists ({size_mb:.2f} MB).")

print("\n" + "=" * 65)
print("RAW RNA-SEQ SEQUENCING READS (NCBI SRA / ENA)")
print("BioProject: PRJNA670942 | SRA Study: SRP288335")
print("=" * 65)
print("Samples can be downloaded via NCBI SRA Toolkit:")
print("  fasterq-dump --split-files <SRX_ACCESSION> -O data/fastq/")
print("  (See data/sample_metadata_sra.csv for mapping of all 30 samples to SRX accessions)\n")

print("=" * 65)
print("RAW TMT PROTEOMICS MASS SPECTROMETRY (PRIDE)")
print("PXD027464 (Aging Muscle) | PXD027490 (Dexamethasone Muscle)")
print("=" * 65)
print("Download raw .raw and .pepXML files via FTP/HTTPS:")
print("  ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027464/")
print("  ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD027490/\n")

print("=" * 65)
print("SUPPLEMENTARY TABLES & DATASETS (Cell Reports / PMC8852763)")
print("=" * 65)
print("Download Table S1-S6 spreadsheets directly from:")
print("  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8852763/#SM1")
print("  Save spreadsheets into: data/supplementary/")
print("=" * 65)

