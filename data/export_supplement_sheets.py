import os
import pandas as pd

supp_dir = 'data/supplementary'

supp_files = {
    'Table_S1_RNAseq.xlsx': 'NIHMS1756261-supplement-2.xlsx',
    'Table_S2_GSEA.xlsx': 'NIHMS1756261-supplement-3.xlsx',
    'Table_S3_Proteomics.xlsx': 'NIHMS1756261-supplement-4.xlsx',
    'Table_S4_CellMarkers.xlsx': 'NIHMS1756261-supplement-5.xlsx',
    'Table_S5_Cyr61siRNA.xlsx': 'NIHMS1756261-supplement-6.xlsx',
    'Table_S6_Atroproteins.xlsx': 'NIHMS1756261-supplement-7.xlsx',
}

for prefix, fname in supp_files.items():
    fpath = os.path.join(supp_dir, fname)
    if os.path.exists(fpath):
        print(f"Reading {fname}...")
        xl = pd.ExcelFile(fpath)
        for sname in xl.sheet_names:
            clean_sname = sname.replace(' ', '_').replace('-', '_')
            out_csv = os.path.join(supp_dir, f"{prefix.replace('.xlsx', '')}_{clean_sname}.csv")
            try:
                df = xl.parse(sname)
                df.to_csv(out_csv, index=False)
                print(f"  Exported sheet '{sname}' ({df.shape[0]} rows, {df.shape[1]} cols) -> {out_csv}")
            except Exception as e:
                print(f"  Error exporting '{sname}': {e}")

print("\nAll supplementary sheets exported to CSV successfully!")

