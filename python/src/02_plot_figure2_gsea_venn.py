"""
Script: 02_plot_figure2_gsea_venn.py
Purpose: Reproduce Figure 2 GSEA NES Barplots, 3-Way Venn Diagrams, and 4-Panel Atrogenes in Python
Paper: Hunt et al. (Cell Reports 2021)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
supp_dir = os.path.join(repo_root, "data", "supplementary")
fig_dir = os.path.join(repo_root, "python", "figures")
os.makedirs(fig_dir, exist_ok=True)

# ------------------------------------------------------------------------------
# 1. GSEA NES HORIZONTAL BARPLOTS (FIGURE 2A, 2B, 2C)
# ------------------------------------------------------------------------------
def plot_gsea_condition(ax, csv_path, title):
    df = pd.read_csv(csv_path, skiprows=4)
    df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
    df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
    df = df.dropna(subset=['NES', 'Description'])

    up = df[df['NES'] > 0].sort_values('NES')
    down = df[df['NES'] < 0].sort_values('NES')
    plot_df = pd.concat([down, up]).reset_index(drop=True)

    y_pos = np.arange(len(plot_df))
    colors = []
    for _, r in plot_df.iterrows():
        if r['NES'] > 0:
            colors.append('#2b5c8f' if r['FDR'] <= 0.05 else '#92b6d8')
        else:
            colors.append('#e66101' if r['FDR'] <= 0.05 else '#fdb863')

    ax.barh(y_pos, plot_df['NES'], color=colors, edgecolor='white', height=0.7)
    ax.axvline(0, color='black', linewidth=1)
    ax.axvline(-1, color='gray', linestyle=':', linewidth=0.7)
    ax.axvline(1, color='gray', linestyle=':', linewidth=0.7)

    for i, r in plot_df.iterrows():
        desc = r['Description']
        if len(desc) > 42:
            desc = desc[:40] + '...'
        val = r['NES']
        if val < 0:
            ax.text(0.05, i, desc, va='center', ha='left', fontsize=7.5, fontweight='bold')
        else:
            ax.text(-0.05, i, desc, va='center', ha='right', fontsize=7.5, fontweight='bold')

    ax.set_yticks([])
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=9, fontweight='bold')
    ax.set_title(f'GSEA: {title}', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig, axes = plt.subplots(3, 1, figsize=(9, 12), dpi=300)
plot_gsea_condition(axes[0], os.path.join(supp_dir, 'Table_S2_GSEA_GSEA_RNA_Dexa.csv'), 'Dexamethasone')
plot_gsea_condition(axes[1], os.path.join(supp_dir, 'Table_S2_GSEA_GSEA_RNA_LLC.csv'), 'LLC Cancer Cachexia')
plot_gsea_condition(axes[2], os.path.join(supp_dir, 'Table_S2_GSEA_GSEA_RNA_Aging.csv'), 'Aging (Sarcopenia)')

plt.tight_layout()
gsea_out = os.path.join(fig_dir, 'Figure2_GSEA_NES_Bars.png')
plt.savefig(gsea_out)
plt.close()
print(f"Saved: {gsea_out}")

# ------------------------------------------------------------------------------
# 2. VENN DIAGRAMS (FIGURE 2D)
# ------------------------------------------------------------------------------
# Matplotlib fallback 3-circle Venn using patches if matplotlib_venn not installed
from matplotlib.patches import Circle

def draw_simple_venn(ax, title, counts):
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-3.5, 3.8)
    ax.axis('off')
    ax.set_title(title, fontsize=13, fontweight='bold')

    c_aging = Circle((0, 1.2), 1.45, facecolor='#619cff', alpha=0.4, edgecolor='#619cff', linewidth=2)
    c_dex = Circle((-1.2, -0.7), 1.65, facecolor='#00ba38', alpha=0.4, edgecolor='#00ba38', linewidth=2)
    c_llc = Circle((1.2, -0.7), 1.45, facecolor='#f8766d', alpha=0.4, edgecolor='#f8766d', linewidth=2)

    ax.add_patch(c_aging)
    ax.add_patch(c_dex)
    ax.add_patch(c_llc)

    # Labels
    ax.text(0, 1.8, f"Aging\n{counts['Aging']}", ha='center', va='center', fontweight='bold', fontsize=10)
    ax.text(-1.6, -0.8, f"Dexamethasone\n{counts['Dex']}", ha='center', va='center', fontweight='bold', fontsize=10)
    ax.text(1.5, -0.8, f"LLC\n{counts['LLC']}", ha='center', va='center', fontweight='bold', fontsize=10)

    # Overlaps
    ax.text(-0.7, 0.45, str(counts['Dex_Aging']), ha='center', va='center', fontweight='bold', fontsize=10)
    ax.text(0.7, 0.45, str(counts['LLC_Aging']), ha='center', va='center', fontweight='bold', fontsize=10)
    ax.text(0, -1.2, str(counts['Dex_LLC']), ha='center', va='center', fontweight='bold', fontsize=10)
    ax.text(0, -0.2, str(counts['All_3']), ha='center', va='center', fontweight='bold', fontsize=11)

fig, axes = plt.subplots(1, 2, figsize=(11, 5.5), dpi=300)
draw_simple_venn(axes[0], "UP-REGULATED", {'Aging': 263, 'Dex': 942, 'LLC': 124, 'Dex_Aging': 34, 'LLC_Aging': 5, 'Dex_LLC': 30, 'All_3': 1})
draw_simple_venn(axes[1], "DOWN-REGULATED", {'Aging': 9, 'Dex': 203, 'LLC': 131, 'Dex_Aging': 0, 'LLC_Aging': 0, 'Dex_LLC': 16, 'All_3': 1})

plt.tight_layout()
venn_out = os.path.join(fig_dir, 'Figure2D_Venn_Diagrams.png')
plt.savefig(venn_out)
plt.close()
print(f"Saved: {venn_out}")

# ------------------------------------------------------------------------------
# 3. 4-PANEL INDIVIDUAL ATROGENES (FIGURE 2F)
# ------------------------------------------------------------------------------
raw_2f = pd.read_csv(os.path.join(supp_dir, 'Table_S6_Fig2F_Atrogenes_Individual.csv'))

row_starts = {'Fbxo32': 0, 'Trim63': 5, 'Ctsl': 10, 'Psmd11': 15}

def get_gene_data(gene):
    idx = row_starts[gene]
    aging_row = raw_2f.iloc[idx, 1:10].dropna().values.astype(float)
    dex_row = raw_2f.iloc[idx+1, 1:10].dropna().values.astype(float)
    llc_row = raw_2f.iloc[idx+2, 1:10].dropna().values.astype(float)

    return {
        'Aging': (aging_row[:3], aging_row[3:]),
        'Dexamethasone': (dex_row[:4], dex_row[4:]),
        'LLC': (llc_row[:3], llc_row[3:])
    }

def draw_atrogene(ax, gene_data, gene_name, y_max=4.0):
    conditions = ['Aging', 'Dexamethasone', 'LLC']
    x = np.arange(len(conditions))
    width = 0.35

    for i, cond in enumerate(conditions):
        c, t = gene_data[cond]
        mc, mt = np.mean(c), np.mean(t)
        sec = np.std(c, ddof=1) / np.sqrt(len(c))
        set_ = np.std(t, ddof=1) / np.sqrt(len(t))

        p = ttest_ind(c, t).pvalue
        p_str = "P<=0.001" if p < 0.001 else f"P={p:.3f}"

        ax.bar(x[i] - width/2, mc, width=width*0.85, color='white', edgecolor='black', linewidth=1.2, yerr=sec, capsize=3)
        ax.bar(x[i] + width/2, mt, width=width*0.85, color='white', edgecolor='#e41a1c', linewidth=1.2, yerr=set_, capsize=3, ecolor='#e41a1c')

        np.random.seed(42 + i)
        ax.scatter(np.random.normal(x[i] - width/2, 0.04, len(c)), c, color='grey', edgecolors='black', s=20, alpha=0.7)
        ax.scatter(np.random.normal(x[i] + width/2, 0.04, len(t)), t, color='#e41a1c', edgecolors='#e41a1c', s=20, alpha=0.8)

        ax.text(x[i], y_max * 0.92, p_str, ha='center', va='center', fontsize=8.5, fontweight='bold')

    ax.axhline(1, color='gray', linestyle='--', linewidth=0.8)
    ax.set_ylim(0, y_max)
    ax.set_xticks(x)
    ax.set_xticklabels(conditions, fontweight='bold', fontsize=9)
    ax.set_ylabel('relative mRNA level', fontweight='bold', fontsize=9)
    ax.set_title(f"${gene_name}$ mRNA", fontweight='bold', fontsize=11, style='italic')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=300)
draw_atrogene(axes[0, 0], get_gene_data('Fbxo32'), 'Fbxo32', y_max=4.0)
axes[0, 0].legend(['Control', 'Atrophy'], loc='upper left', frameon=False, fontsize=8.5)

draw_atrogene(axes[0, 1], get_gene_data('Trim63'), 'Trim63', y_max=4.0)
draw_atrogene(axes[1, 0], get_gene_data('Ctsl'), 'Ctsl', y_max=4.0)
draw_atrogene(axes[1, 1], get_gene_data('Psmd11'), 'Psmd11', y_max=2.0)

plt.tight_layout()
atrogenes_out = os.path.join(fig_dir, 'Figure2F_Individual_Atrogenes.png')
plt.savefig(atrogenes_out)
plt.close()
print(f"Saved: {atrogenes_out}")

print("\nAll Figure 2 GSEA, Venn, and Atrogene charts generated successfully!")
