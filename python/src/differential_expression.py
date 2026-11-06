import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
data_path = os.path.join(repo_root, 'data', 'GSE159952_all.tpm.txt.gz')
results_dir = os.path.join(repo_root, 'python', 'results')
os.makedirs(results_dir, exist_ok=True)

print(f"Loading GSE159952 TPM data from {data_path}...")
df = pd.read_csv(data_path, sep='\t')
print(f"Data loaded: {df.shape[0]} genes, {df.shape[1]} columns")

# Filter out genes with zero expression across all samples
expr_cols = [c for c in df.columns if c != 'ID']
df_filtered = df[df[expr_cols].sum(axis=1) > 0].copy()
print(f"Genes with non-zero expression: {df_filtered.shape[0]}")

# Define experimental comparisons: (name, test_cols, ctrl_cols)
comparisons = {
    'Dexamethasone': (
        ['dexamethasone_1', 'dexamethasone_2', 'dexamethasone_3', 'dexamethasone_4'],
        ['Control_dexamethasone_1', 'Control_dexamethasone_2', 'Control_dexamethasone_3', 'Control_dexamethasone_4']
    ),
    'LLC_Cancer': (
        ['LLC_1', 'LLC_2', 'LLC_3', 'LLC_4'],
        ['Control_LLC_1', 'Control_LLC_2', 'Control_LLC_3', 'Control_LLC_4']
    ),
    'Aging': (
        ['24month_1', '24month_2', '24month_3'],
        ['6month_1', '6month_2', '6month_3']
    ),
    'siCYR61': (
        ['SICYR61_1', 'SICYR61_2', 'SICYR61_3', 'SICYR61_4'],
        ['SINT_1', 'SINT_2', 'SINT_3', 'SINT_4']
    )
}

results = {}

for name, (test_cols, ctrl_cols) in comparisons.items():
    print(f"\nAnalyzing {name} ({len(test_cols)} test vs {len(ctrl_cols)} ctrl)...")
    test_vals = df_filtered[test_cols].values
    ctrl_vals = df_filtered[ctrl_cols].values
    
    mean_test = np.mean(test_vals, axis=1)
    mean_ctrl = np.mean(ctrl_vals, axis=1)
    
    log2fc = np.log2(mean_test + 1e-4) - np.log2(mean_ctrl + 1e-4)
    
    t_stat, p_val = stats.ttest_ind(test_vals, ctrl_vals, axis=1, equal_var=True, nan_policy='omit')
    p_val = np.nan_to_num(p_val, nan=1.0)
    
    res_df = pd.DataFrame({
        'gene': df_filtered['ID'].values,
        'mean_test': mean_test,
        'mean_ctrl': mean_ctrl,
        'log2FC': log2fc,
        'p_value': p_val,
        'minus_log10_p': -np.log10(np.clip(p_val, 1e-300, 1.0))
    })
    
    res_df['is_up'] = (res_df['p_value'] < 0.05) & (res_df['log2FC'] > 1.0)
    res_df['is_down'] = (res_df['p_value'] < 0.05) & (res_df['log2FC'] < -1.0)
    res_df['is_sig'] = res_df['is_up'] | res_df['is_down']
    
    n_up = res_df['is_up'].sum()
    n_down = res_df['is_down'].sum()
    print(f"  Significantly Upregulated (p < 0.05, log2FC > 1): {n_up}")
    print(f"  Significantly Downregulated (p < 0.05, log2FC < -1): {n_down}")
    print(f"  Total Significant DEGs: {n_up + n_down}")
    
    results[name] = res_df
    res_df.to_csv(os.path.join(results_dir, f"{name}_DEGs.csv"), index=False)

print("\nDifferential expression analysis completed successfully!")
