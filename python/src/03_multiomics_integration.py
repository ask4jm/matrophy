"""
Script: 03_multiomics_integration.py
Purpose: Reproduce Figure 4A (RNA-Protein Correlation), Figure 4B (Transcriptome-Proteome Disconnect),
         and Figure 5A (Core Atroproteins) in Python using Matplotlib/Seaborn.
Paper: Hunt et al. (Cell Reports 2021)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    print("=================================================================")
    print("Python Step 3: Multi-Omics Integration & Atroprotein Analysis")
    print("=================================================================")

    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    supp_dir = os.path.join(repo_root, "data", "supplementary")
    results_dir = os.path.join(repo_root, "python", "results")
    figures_dir = os.path.join(repo_root, "python", "figures")

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    cross_csv = os.path.join(supp_dir, "Table_S3_Proteomics_Cross_Comparisons.csv")
    if not os.path.exists(cross_csv):
        raise FileNotFoundError(f"Missing file: {cross_csv}")

    print(f"Loading matched multi-omics matrix: {cross_csv}")
    df = pd.read_csv(cross_csv)
    print(f"Loaded {len(df)} matched genes with both mRNA and Protein measurements.")

    # --------------------------------------------------------------------------
    # 1. FIGURE 4A: RNA vs PROTEIN CORRELATION SCATTER PLOTS
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), sharey=True)

    conditions = [
        ("Dexa", "Dexa RNA", "Dexa Protein", "#e41a1c", "Dexamethasone (Dex)", 0.164),
        ("LLC", "LLC RNA", "LLC Protein", "#377eb8", "LLC Cancer Cachexia", 0.082),
        ("Aging", "Aging RNA", "Aging Protein", "#4daf4a", "Aging (24m vs 6m)", 0.202)
    ]

    for ax, (code, rna_col, prot_col, color, title, target_r2) in zip(axes, conditions):
        sub = df[[rna_col, prot_col]].dropna()
        x = sub[rna_col].values
        y = sub[prot_col].values

        r = np.corrcoef(x, y)[0, 1]
        r2 = r ** 2

        ax.scatter(x, y, alpha=0.25, s=15, color=color, edgecolors="none")

        # Linear trendline
        m, b = np.polyfit(x, y, 1)
        x_vals = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_vals, m * x_vals + b, color="black", lw=1.5, linestyle="--")

        ax.axhline(0, color="gray", lw=0.8, linestyle=":")
        ax.axvline(0, color="gray", lw=0.8, linestyle=":")

        ax.set_title(f"{title}\nr = {r:.3f}, R² = {r2:.3f} (Published: {target_r2})", fontsize=12, fontweight="bold")
        ax.set_xlabel("mRNA log2(Fold Change)", fontsize=11, fontweight="bold")
        if code == "Dexa":
            ax.set_ylabel("Protein log2(Fold Change)", fontsize=11, fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.5)

    plt.suptitle("Figure 4A Reproduction: Transcriptome vs Proteome Correlation across Muscle Atrophy Modes", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig4a_path = os.path.join(figures_dir, "Figure4A_RNA_Protein_Correlation.png")
    plt.savefig(fig4a_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig4a_path}")

    # --------------------------------------------------------------------------
    # 2. FIGURE 4B: TRANSCRIPTOME-PROTEOME DISCONNECT BAR PLOT
    # --------------------------------------------------------------------------
    def bin_disconnect(rna_series, prot_series):
        mask_down = rna_series < -0.5
        mask_mid = (rna_series >= -0.5) & (rna_series <= 0.5)
        mask_up = rna_series > 0.5
        return [
            prot_series[mask_down].mean(),
            prot_series[mask_mid].mean(),
            prot_series[mask_up].mean()
        ]

    disc_dex = bin_disconnect(df["Dexa RNA"], df["Dexa Protein"])
    disc_llc = bin_disconnect(df["LLC RNA"], df["LLC Protein"])
    disc_age = bin_disconnect(df["Aging RNA"], df["Aging Protein"])

    fig, ax = plt.subplots(figsize=(9, 6))
    x_pos = np.arange(3)
    width = 0.25

    ax.bar(x_pos - width, disc_dex, width, label="Dexamethasone", color="#e41a1c", edgecolor="black")
    ax.bar(x_pos, disc_llc, width, label="LLC Cancer", color="#377eb8", edgecolor="black")
    ax.bar(x_pos + width, disc_age, width, label="Aging", color="#4daf4a", edgecolor="black")

    ax.axhline(0, color="black", lw=1)
    ax.axhline(0.1, color="gray", linestyle=":", lw=0.8)
    ax.axhline(-0.1, color="gray", linestyle=":", lw=0.8)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(["mRNA Down (< -0.5)", "mRNA Unchanged (-0.5 to 0.5)", "mRNA Up (> 0.5)"], fontsize=11, fontweight="bold")
    ax.set_ylabel("Mean Protein log2(Fold Change)", fontsize=12, fontweight="bold")
    ax.set_title("Figure 4B Reproduction: Protein Abundance by mRNA Change Bins", fontsize=13, fontweight="bold")
    ax.legend(frameon=True, facecolor="white", edgecolor="gray")
    ax.grid(True, linestyle=":", alpha=0.5, axis="y")

    fig4b_path = os.path.join(figures_dir, "Figure4B_Disconnect_Barplot.png")
    plt.tight_layout()
    plt.savefig(fig4b_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig4b_path}")

    # --------------------------------------------------------------------------
    # 3. FIGURE 5A: CORE ATROPROTEINS ACROSS ALL MODES
    # --------------------------------------------------------------------------
    display_order = ["Prg4", "Apoe", "Ncam1", "Clu", "Fn1", "Ahnak2", "Apod", "Pla2g7", "Ndufb2", "Kera", "Pcolce"]
    sub_df = df[df["GN"].isin(display_order)].copy()
    sub_df["sort_key"] = sub_df["GN"].map({g: i for i, g in enumerate(display_order)})
    sub_df = sub_df.sort_values("sort_key")

    fig, ax = plt.subplots(figsize=(14, 6))
    x_idx = np.arange(len(sub_df))
    w = 0.25

    ax.bar(x_idx - w, sub_df["Aging Protein"], w, label="Aging (24m vs 6m)", color="#4daf4a", edgecolor="black")
    ax.bar(x_idx, sub_df["Dexa Protein"], w, label="Dexamethasone", color="#e41a1c", edgecolor="black")
    ax.bar(x_idx + w, sub_df["LLC Protein"], w, label="LLC Cancer Cachexia", color="#377eb8", edgecolor="black")

    ax.axhline(0, color="black", lw=1)
    ax.axhline(0.3, color="gray", linestyle="--", lw=0.8)
    ax.axhline(-0.3, color="gray", linestyle="--", lw=0.8)

    ax.set_xticks(x_idx)
    ax.set_xticklabels(sub_df["GN"], rotation=45, ha="right", fontsize=11, fontweight="bold")
    ax.set_ylabel("Protein log2(Fold Change)", fontsize=12, fontweight="bold")
    ax.set_title("Figure 5A Reproduction: Consistently Regulated 'Atroproteins' Across 3 Modes", fontsize=14, fontweight="bold")
    ax.legend(frameon=True, facecolor="white", edgecolor="gray")
    ax.grid(True, linestyle=":", alpha=0.5, axis="y")

    fig5a_path = os.path.join(figures_dir, "Figure5A_Top_Atroproteins.png")
    plt.tight_layout()
    plt.savefig(fig5a_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig5a_path}")

    print("03_multiomics_integration.py completed successfully!")

if __name__ == "__main__":
    main()

