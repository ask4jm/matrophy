"""
Script: 05_muscle_mass_correlation.py
Purpose: Reproduce Figure 7B — Correlation of Protein Abundance with Muscle Mass in Aging Atrophy
Paper: Hunt et al. (Cell Reports 2021)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("=================================================================")
    print("Python Step 4: Proteome Correlation with Muscle Mass (Figure 7B)")
    print("=================================================================")

    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    supp_dir = os.path.join(repo_root, "data", "supplementary")
    results_dir = os.path.join(repo_root, "python", "results")
    figures_dir = os.path.join(repo_root, "python", "figures")

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    corr_csv = os.path.join(supp_dir, "Table_S3_Proteomics_Aging_Muscle_Correlation.csv")
    if not os.path.exists(corr_csv):
        raise FileNotFoundError(f"Missing file: {corr_csv}")

    print(f"Loading Aging-Muscle correlation dataset: {corr_csv}")
    corr_df = pd.read_csv(corr_csv)
    print(f"Loaded {len(corr_df)} proteins with muscle mass correlation metrics.")

    # Sort to find top negative and top positive
    top_neg = corr_df.sort_values("ave correlation").head(10)
    top_pos = corr_df.sort_values("ave correlation", ascending=False).head(10)

    print("\nTop 5 Negative Correlating Proteins (Inversely related to muscle mass):")
    print(top_neg[["GN", "raw TA mass correlation", "TA/tibia correlation", "ave correlation"]].head(5))

    out_csv = os.path.join(results_dir, "Figure7B_Top_Muscle_Mass_Correlating_Proteins.csv")
    pd.concat([top_neg, top_pos]).to_csv(out_csv, index=False)

    # Plot Figure 7B
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: Top 10 Inversely Correlated
    axes[0].bar(top_neg["GN"], top_neg["ave correlation"], color="#e41a1c", edgecolor="black")
    axes[0].axhline(0, color="black", lw=1)
    axes[0].axhline(-0.9, color="gray", linestyle="--", lw=0.8)
    axes[0].set_ylabel("Pearson Correlation with Muscle Mass (r)", fontsize=11, fontweight="bold")
    axes[0].set_title("Top Proteins Inversely Correlated\nwith Muscle Mass (Hunt et al. Fig 7B)", fontsize=12, fontweight="bold")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].set_ylim(-1.0, 0)
    axes[0].grid(True, linestyle=":", alpha=0.5, axis="y")

    # Panel 2: Proteome-wide correlation distribution
    axes[1].hist(corr_df["ave correlation"].dropna(), bins=30, color="#377eb8", edgecolor="white")
    mean_r = corr_df["ave correlation"].mean()
    axes[1].axvline(0, color="black", lw=1)
    axes[1].axvline(mean_r, color="gold", linestyle="--", lw=2, label=f"Mean r = {mean_r:.2f}")
    axes[1].set_xlabel("Average Correlation with Muscle Mass", fontsize=11, fontweight="bold")
    axes[1].set_ylabel("Number of Proteins", fontsize=11, fontweight="bold")
    axes[1].set_title(f"Proteome-Wide Muscle Mass\nCorrelation Distribution (N = {len(corr_df)})", fontsize=12, fontweight="bold")
    axes[1].legend()
    axes[1].grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig7b_path = os.path.join(figures_dir, "Figure7B_Muscle_Mass_Correlation.png")
    plt.savefig(fig7b_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig7b_path}")

    print("05_muscle_mass_correlation.py completed successfully!")

if __name__ == "__main__":
    main()

