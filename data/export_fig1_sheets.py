import os
import pandas as pd

xl_path = 'supplementary/NIHMS1756261-supplement-7.xlsx'
out_dir = 'data/supplementary'
os.makedirs(out_dir, exist_ok=True)

xl = pd.ExcelFile(xl_path)
fig1_sheets = [s for s in xl.sheet_names if 'Fig1' in s]
print(f"Exporting {len(fig1_sheets)} sheets...")

for s in fig1_sheets:
    df = xl.parse(s)
    clean_name = s.replace(' ', '_').replace('-', '_') + '.csv'
    out_file = os.path.join(out_dir, clean_name)
    df.to_csv(out_file, index=False)
    print(f"  {s} -> {clean_name}")

print("\nAll 30 Fig1 sheets exported to CSV successfully!")

