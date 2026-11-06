import gzip
import pandas as pd

titles, gsms, sras = [], [], []
with gzip.open('data/GSE159952_series_matrix.txt.gz', 'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            titles = [x.strip(' \t\n\r"').replace('Tibialis anterior ', '') for x in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_geo_accession'):
            gsms = [x.strip(' \t\n\r"') for x in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_relation') and 'sra?term=' in line:
            sras = [x.strip(' \t\n\r"').split('term=')[-1] for x in line.strip().split('\t')[1:]]

df = pd.DataFrame({'Sample_Name': titles, 'GEO_Accession': gsms, 'SRA_Accession': sras})
df.to_csv('data/sample_metadata_sra.csv', index=False)
print(f"Saved data/sample_metadata_sra.csv with {len(df)} samples")
print(df.to_string())

