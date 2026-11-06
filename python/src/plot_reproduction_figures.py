import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
results_dir = os.path.join(repo_root, 'python', 'results')
figures_dir = os.path.join(repo_root, 'python', 'figures')
data_path = os.path.join(repo_root, 'data', 'GSE159952_all.tpm.txt.gz')

os.makedirs(figures_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# 1. LOAD DIFFERENTIAL EXPRESSION RESULTS
dex_df = pd.read_csv(os.path.join(results_dir, 'Dexamethasone_DEGs.csv'))
llc_df = pd.read_csv(os.path.join(results_dir, 'LLC_Cancer_DEGs.csv'))
age_df = pd.read_csv(os.path.join(results_dir, 'Aging_DEGs.csv'))
cyr_df = pd.read_csv(os.path.join(results_dir, 'siCYR61_DEGs.csv'))

# 2. FIGURE 2A-C: VOLCANO PLOTS
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), dpi=300)

datasets = [
    ('Dexamethasone vs Control', dex_df, axes[0]),
    ('LLC Cancer Cachexia vs Control', llc_df, axes[1]),
    ('Aging (24m vs 6m)', age_df, axes[2])
]

highlight_genes = ['Trim63', 'Fbxo32', 'Cyr61', 'Psmd11', 'Foxo1', 'Foxo3', 'S100a9', 'Il1b', 'Cdkn1a', 'Ctsl', 'Apod', 'Apoe']

for title, df, ax in datasets:
    capped_p = np.clip(df['minus_log10_p'], 0, 15)
    
    ax.scatter(df.loc[~df['is_sig'], 'log2FC'], capped_p[~df['is_sig']],
               c='lightgray', s=8, alpha=0.5, label='Not Significant')
    
    ax.scatter(df.loc[df['is_up'], 'log2FC'], capped_p[df['is_up']],
               c='#d95f02', s=12, alpha=0.7, label=f'Up ({df["is_up"].sum()})')
    
    ax.scatter(df.loc[df['is_down'], 'log2FC'], capped_p[df['is_down']],
               c='#1f78b4', s=12, alpha=0.7, label=f'Down ({df["is_down"].sum()})')
    
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(1.0, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(-1.0, color='gray', linestyle='--', linewidth=0.8)
    
    for g in highlight_genes:
        row = df[df['gene'] == g]
        if not row.empty:
            r = row.iloc[0]
            if r['is_sig']:
                ax.scatter(r['log2FC'], min(r['minus_log10_p'], 15), c='black', s=35, zorder=5)
                ax.annotate(g, (r['log2FC'], min(r['minus_log10_p'], 15)),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize=8, fontweight='bold',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none'))

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('log2(Fold Change)', fontsize=11)
    ax.set_ylabel('-log10(p-value)', fontsize=11)
    ax.set_xlim(-8, 8)
    ax.set_ylim(-0.5, 16)
    ax.legend(loc='upper left', frameon=True, fontsize=8)

plt.tight_layout()
fig2_path = os.path.join(figures_dir, 'Figure2_VolcanoPlots.png')
plt.savefig(fig2_path)
plt.close()
print(f"Figure 2 Volcano Plots saved to {fig2_path}")

# 3. FIGURE 2D: OVERLAP & BARPLOT
up_dex = set(dex_df.loc[dex_df['is_up'], 'gene'])
up_llc = set(llc_df.loc[llc_df['is_up'], 'gene'])
up_age = set(age_df.loc[age_df['is_up'], 'gene'])

down_dex = set(dex_df.loc[dex_df['is_down'], 'gene'])
down_llc = set(llc_df.loc[llc_df['is_down'], 'gene'])
down_age = set(age_df.loc[age_df['is_down'], 'gene'])

all_3_up = up_dex & up_llc & up_age
dex_llc_up = (up_dex & up_llc) - all_3_up
dex_age_up = (up_dex & up_age) - all_3_up
llc_age_up = (up_llc & up_age) - all_3_up
unique_dex_up = up_dex - up_llc - up_age
unique_llc_up = up_llc - up_dex - up_age
unique_age_up = up_age - up_dex - up_llc

all_3_down = down_dex & down_llc & down_age
dex_llc_down = (down_dex & down_llc) - all_3_down
dex_age_down = (down_dex & down_age) - all_3_down
llc_age_down = (down_llc & down_age) - all_3_down
unique_dex_down = down_dex - down_llc - down_age
unique_llc_down = down_llc - down_dex - down_age
unique_age_down = down_age - down_dex - down_llc

overlap_summary = pd.DataFrame([
    {'Category': 'Upregulated (p<0.05, log2FC>1)',
     'Total Dex': len(up_dex), 'Total LLC': len(up_llc), 'Total Aging': len(up_age),
     'Unique Dex': len(unique_dex_up), 'Unique LLC': len(unique_llc_up), 'Unique Aging': len(unique_age_up),
     'Dex & LLC': len(dex_llc_up), 'Dex & Aging': len(dex_age_up), 'LLC & Aging': len(llc_age_up),
     'Shared by All 3': len(all_3_up)},
    {'Category': 'Downregulated (p<0.05, log2FC<-1)',
     'Total Dex': len(down_dex), 'Total LLC': len(down_llc), 'Total Aging': len(down_age),
     'Unique Dex': len(unique_dex_down), 'Unique LLC': len(unique_llc_down), 'Unique Aging': len(unique_age_down),
     'Dex & LLC': len(dex_llc_down), 'Dex & Aging': len(dex_age_down), 'LLC & Aging': len(llc_age_down),
     'Shared by All 3': len(all_3_down)}
])
overlap_summary.to_csv(os.path.join(results_dir, 'Figure2D_Overlap_Summary.csv'), index=False)

fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
categories = ['Dex Upregulated', 'LLC Upregulated', 'Aging Upregulated',
              'Dex Downregulated', 'LLC Downregulated', 'Aging Downregulated']
unique_vals = [len(unique_dex_up), len(unique_llc_up), len(unique_age_up),
               len(unique_dex_down), len(unique_llc_down), len(unique_age_down)]
shared_2_vals = [len(dex_llc_up)+len(dex_age_up), len(dex_llc_up)+len(llc_age_up), len(dex_age_up)+len(llc_age_up),
                 len(dex_llc_down)+len(dex_age_down), len(dex_llc_down)+len(llc_age_down), len(dex_age_down)+len(llc_age_down)]
shared_3_vals = [len(all_3_up), len(all_3_up), len(all_3_up),
                 len(all_3_down), len(all_3_down), len(all_3_down)]

x = np.arange(len(categories))
width = 0.55

ax.bar(x, unique_vals, width, label='Stimulus-Specific (Unique)', color='#4daf4a')
ax.bar(x, shared_2_vals, width, bottom=unique_vals, label='Shared by 2 Stimuli', color='#377eb8')
bottom_3 = [u + s2 for u, s2 in zip(unique_vals, shared_2_vals)]
ax.bar(x, shared_3_vals, width, bottom=bottom_3, label='Shared by All 3 Stimuli', color='#e41a1c')

ax.set_ylabel('Number of DEGs (p < 0.05, |log2FC| > 1)', fontsize=11)
ax.set_title('Figure 2D Reproduction: Stimulus-Specific vs Shared Transcriptional Programs', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=25, ha='right', fontsize=9)
ax.legend(loc='upper right', frameon=True)

for i in range(len(categories)):
    total = unique_vals[i] + shared_2_vals[i] + shared_3_vals[i]
    if total > 0:
        pct_unique = (unique_vals[i] / total) * 100
        ax.annotate(f"{pct_unique:.1f}%\nUnique", (x[i], total + 20),
                    ha='center', fontsize=8, fontweight='bold', color='#2b6329')

plt.tight_layout()
fig2d_path = os.path.join(figures_dir, 'Figure2D_Stimulus_Specific_vs_Shared.png')
plt.savefig(fig2d_path)
plt.close()

# 4. FIGURE 2F: CLASSIC ATROGENES
classic_atrogenes = [
    'Trim63', 'Fbxo32', 'Psmd11', 'Fbxo30', 'Fbxo21', 'Foxo1', 'Foxo3',
    'Ctsl', 'Cdkn1a', 'Gadd45a', 'Mstn', 'Casp3', 'Ubr5', 'S100a9', 'Il1b', 'Cyr61'
]

atrogene_records = []
for gene in classic_atrogenes:
    dex_row = dex_df[dex_df['gene'] == gene]
    llc_row = llc_df[llc_df['gene'] == gene]
    age_row = age_df[age_df['gene'] == gene]
    
    dex_fc = dex_row['log2FC'].values[0] if not dex_row.empty else 0.0
    dex_p = dex_row['p_value'].values[0] if not dex_row.empty else 1.0
    llc_fc = llc_row['log2FC'].values[0] if not llc_row.empty else 0.0
    llc_p = llc_row['p_value'].values[0] if not llc_row.empty else 1.0
    age_fc = age_row['log2FC'].values[0] if not age_row.empty else 0.0
    age_p = age_row['p_value'].values[0] if not age_row.empty else 1.0
    
    atrogene_records.append({
        'Gene': gene,
        'Dex_log2FC': dex_fc, 'Dex_pval': dex_p,
        'LLC_log2FC': llc_fc, 'LLC_pval': llc_p,
        'Aging_log2FC': age_fc, 'Aging_pval': age_p
    })

atrogene_res_df = pd.DataFrame(atrogene_records)
atrogene_res_df.to_csv(os.path.join(results_dir, 'Figure2F_Classic_Atrogenes.csv'), index=False)

fig, ax = plt.subplots(figsize=(12, 6), dpi=300)
x_pos = np.arange(len(classic_atrogenes))
bar_width = 0.25

ax.bar(x_pos - bar_width, atrogene_res_df['Dex_log2FC'], bar_width, label='Dexamethasone', color='#e41a1c', alpha=0.85)
ax.bar(x_pos, atrogene_res_df['LLC_log2FC'], bar_width, label='LLC Cancer Cachexia', color='#377eb8', alpha=0.85)
ax.bar(x_pos + bar_width, atrogene_res_df['Aging_log2FC'], bar_width, label='Aging (24m vs 6m)', color='#4daf4a', alpha=0.85)

ax.axhline(0, color='black', linewidth=0.8)
ax.axhline(1.0, color='gray', linestyle=':', linewidth=0.7)
ax.axhline(-1.0, color='gray', linestyle=':', linewidth=0.7)

ax.set_ylabel('log2(Fold Change)', fontsize=11)
ax.set_title('Figure 2F Reproduction: Expression of Classic Atrogenes Across Distinct Atrophy Stimuli', fontsize=12, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(classic_atrogenes, rotation=35, ha='right', fontsize=9, fontweight='medium')
ax.legend(loc='upper right', frameon=True)

for i, r in atrogene_res_df.iterrows():
    if r['Dex_pval'] < 0.05:
        y = r['Dex_log2FC'] + (0.1 if r['Dex_log2FC'] >= 0 else -0.3)
        ax.text(x_pos[i] - bar_width, y, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
    if r['LLC_pval'] < 0.05:
        y = r['LLC_log2FC'] + (0.1 if r['LLC_log2FC'] >= 0 else -0.3)
        ax.text(x_pos[i], y, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
    if r['Aging_pval'] < 0.05:
        y = r['Aging_log2FC'] + (0.1 if r['Aging_log2FC'] >= 0 else -0.3)
        ax.text(x_pos[i] + bar_width, y, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
fig2f_path = os.path.join(figures_dir, 'Figure2F_Classic_Atrogenes.png')
plt.savefig(fig2f_path)
plt.close()

# 5. FIGURE S1C: PCA
from scipy.linalg import svd
tpm_data = pd.read_csv(data_path, sep='\t')
sample_cols = [c for c in tpm_data.columns if c != 'ID']
log_expr = np.log2(tpm_data[sample_cols].values + 1.0)
gene_vars = np.var(log_expr, axis=1)
top_idx = np.argsort(gene_vars)[-5000:]
filtered_matrix = log_expr[top_idx, :]
mean_centered = filtered_matrix - np.mean(filtered_matrix, axis=1, keepdims=True)

U, S, Vt = svd(mean_centered, full_matrices=False)
pc1 = Vt[0, :] * S[0]
pc2 = Vt[1, :] * S[1]
var_explained = (S ** 2) / np.sum(S ** 2) * 100

pca_df = pd.DataFrame({'Sample': sample_cols, 'PC1': pc1, 'PC2': pc2})

def assign_group(s):
    if '24month' in s: return 'Aging (24-Month)'
    elif '6month' in s: return 'Young (6-Month)'
    elif s.startswith('dexamethasone'): return 'Dexamethasone'
    elif 'Control_dexamethasone' in s: return 'Control (Dex)'
    elif s.startswith('LLC'): return 'LLC Cachexia'
    elif 'Control_LLC' in s: return 'Control (LLC)'
    elif 'SICYR61' in s: return 'siCYR61'
    elif 'SINT' in s: return 'siNT (Control)'
    return 'Other'

pca_df['Group'] = pca_df['Sample'].apply(assign_group)
pca_df.to_csv(os.path.join(results_dir, 'FigureS1C_PCA_Coordinates.csv'), index=False)

fig, ax = plt.subplots(figsize=(9, 7), dpi=300)
colors = {
    'Aging (24-Month)': '#e41a1c', 'Young (6-Month)': '#377eb8',
    'Dexamethasone': '#ff7f00', 'Control (Dex)': '#fdbf6f',
    'LLC Cachexia': '#984ea3', 'Control (LLC)': '#cab2d6',
    'siCYR61': '#4daf4a', 'siNT (Control)': '#b2df8a'
}
markers = {
    'Aging (24-Month)': 'o', 'Young (6-Month)': 'o',
    'Dexamethasone': 's', 'Control (Dex)': 's',
    'LLC Cachexia': '^', 'Control (LLC)': '^',
    'siCYR61': 'D', 'siNT (Control)': 'D'
}

for group, grp_df in pca_df.groupby('Group'):
    ax.scatter(grp_df['PC1'], grp_df['PC2'],
               c=colors.get(group, 'black'), marker=markers.get(group, 'o'),
               s=90, alpha=0.85, edgecolors='black', linewidth=0.8, label=group)

ax.set_xlabel(f'PC1 ({var_explained[0]:.1f}% Variance)', fontsize=11)
ax.set_ylabel(f'PC2 ({var_explained[1]:.1f}% Variance)', fontsize=11)
ax.set_title('Figure S1C Reproduction: PCA of Skeletal Muscle RNA-seq Profiles', fontsize=12, fontweight='bold')
ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left', frameon=True)
ax.grid(True, linestyle=':', alpha=0.5)

plt.tight_layout()
figs1c_path = os.path.join(figures_dir, 'FigureS1C_PCA.png')
plt.savefig(figs1c_path)
plt.close()

# 6. FIGURE 6: CYR61 / WNT PATHWAY
wnt_targets = ['Wnt6', 'Wnt7b', 'Sfrp5', 'Rspo1', 'Kremen2', 'Wnt4', 'Fzd1', 'Fzd4', 'Dkk3', 'Axin2', 'Ctnnb1']
wnt_records = []
for gene in wnt_targets:
    row = cyr_df[cyr_df['gene'] == gene]
    if not row.empty:
        r = row.iloc[0]
        wnt_records.append({
            'Gene': gene, 'log2FC_siCYR61': r['log2FC'],
            'pval_siCYR61': r['p_value'], 'minus_log10_p': r['minus_log10_p'], 'is_sig': r['is_sig']
        })

wnt_res_df = pd.DataFrame(wnt_records)
wnt_res_df.to_csv(os.path.join(results_dir, 'Figure6_Wnt_Pathways_siCYR61.csv'), index=False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
capped_p_cyr = np.clip(cyr_df['minus_log10_p'], 0, 10)
ax1.scatter(cyr_df.loc[~cyr_df['is_sig'], 'log2FC'], capped_p_cyr[~cyr_df['is_sig']], c='lightgray', s=8, alpha=0.5, label='Not Sig')
ax1.scatter(cyr_df.loc[cyr_df['is_up'], 'log2FC'], capped_p_cyr[cyr_df['is_up']], c='#d95f02', s=15, alpha=0.7, label=f'Up ({cyr_df["is_up"].sum()})')
ax1.scatter(cyr_df.loc[cyr_df['is_down'], 'log2FC'], capped_p_cyr[cyr_df['is_down']], c='#1f78b4', s=15, alpha=0.7, label=f'Down ({cyr_df["is_down"].sum()})')

ax1.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8)
ax1.axvline(1.0, color='gray', linestyle='--', linewidth=0.8)
ax1.axvline(-1.0, color='gray', linestyle='--', linewidth=0.8)

for g in ['Cyr61', 'Wnt6', 'Wnt7b', 'Sfrp5', 'Rspo1', 'Kremen2']:
    row = cyr_df[cyr_df['gene'] == g]
    if not row.empty:
        r = row.iloc[0]
        ax1.scatter(r['log2FC'], min(r['minus_log10_p'], 10), c='black', s=40, zorder=5)
        ax1.annotate(g, (r['log2FC'], min(r['minus_log10_p'], 10)),
                     xytext=(5, 5), textcoords='offset points', fontsize=9, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none'))

ax1.set_title('Figure 6F: siCYR61 vs siNT in Aged Muscle', fontsize=11, fontweight='bold')
ax1.set_xlabel('log2(Fold Change)', fontsize=10)
ax1.set_ylabel('-log10(p-value)', fontsize=10)
ax1.set_xlim(-6, 6)
ax1.legend(loc='upper right', fontsize=8)

x_wnt = np.arange(len(wnt_res_df))
ax2.bar(x_wnt, wnt_res_df['log2FC_siCYR61'], color='#7570b3', alpha=0.85, width=0.5)
ax2.axhline(0, color='black', linewidth=0.8)
ax2.axhline(1.0, color='gray', linestyle=':', linewidth=0.7)
ax2.axhline(-1.0, color='gray', linestyle=':', linewidth=0.7)

ax2.set_xticks(x_wnt)
ax2.set_xticklabels(wnt_res_df['Gene'], rotation=35, ha='right', fontsize=9)
ax2.set_ylabel('log2(Fold Change: siCYR61 / siNT)', fontsize=10)
ax2.set_title('Figure 6H: Modulation of Wnt Pathway Components', fontsize=11, fontweight='bold')

for i, r in wnt_res_df.iterrows():
    if r['pval_siCYR61'] < 0.05:
        y = r['log2FC_siCYR61'] + (0.05 if r['log2FC_siCYR61'] >= 0 else -0.15)
        ax2.text(x_wnt[i], y, '*', ha='center', va='bottom', fontsize=11, fontweight='bold')

plt.tight_layout()
fig6_path = os.path.join(figures_dir, 'Figure6_Cyr61_Knockdown_Wnt.png')
plt.savefig(fig6_path)
plt.close()

print("All Python reproduction figures generated successfully!")
