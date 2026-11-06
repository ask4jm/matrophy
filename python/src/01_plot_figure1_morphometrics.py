"""
Script: 01_plot_figure1_morphometrics.py
Purpose: Reproduce Figure 1A, 1B, and 1C (Physiological & Histological Phenotyping) in Python
Paper: Hunt et al. (Cell Reports 2021)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def draw_single_bar(ax, ctrl_vals, trt_vals, title, ylab, trt_name, trt_color):
    c = pd.Series(ctrl_vals).dropna().values.astype(float)
    t = pd.Series(trt_vals).dropna().values.astype(float)

    m_c, m_t = np.mean(c), np.mean(t)
    se_c = np.std(c, ddof=1) / np.sqrt(len(c)) if len(c) > 1 else 0
    se_t = np.std(t, ddof=1) / np.sqrt(len(t)) if len(t) > 1 else 0

    p_val = ttest_ind(c, t).pvalue
    p_str = "P<=0.001" if p_val < 0.001 else f"P={p_val:.3f}"
    pct = ((m_t - m_c) / m_c) * 100
    pct_str = f"{pct:+.0f}%"

    x_pos = [0, 1]
    ax.bar(0, m_c, width=0.6, color="white", edgecolor="black", linewidth=1.2, yerr=se_c, capsize=4)
    ax.bar(1, m_t, width=0.6, color="white", edgecolor=trt_color, linewidth=1.2, yerr=se_t, capsize=4, ecolor=trt_color)

    # Jittered data points
    np.random.seed(42)
    jit_c = np.random.normal(0, 0.06, size=len(c))
    jit_t = np.random.normal(1, 0.06, size=len(t))

    ax.scatter(jit_c, c, color="grey", edgecolors="black", s=20, alpha=0.7, zorder=3)
    ax.scatter(jit_t, t, color=trt_color, edgecolors=trt_color, s=20, alpha=0.8, zorder=3)

    y_max = max(np.max(c), np.max(t), m_c + se_c, m_t + se_t) * 1.25
    ax.set_ylim(0, y_max)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(["Ctrl", trt_name], fontweight="bold")
    ax.set_title(title, fontweight="bold", fontsize=10)
    ax.set_ylabel(ylab, fontsize=9, fontweight="bold")

    # Annotate p-value at top
    ax.text(0.5, y_max * 0.90, p_str, ha="center", va="center", fontsize=8.5, fontweight="bold")
    # Annotate percent change at bottom right
    ax.text(1.25, y_max * 0.08, pct_str, ha="right", va="bottom", fontsize=8.5, fontweight="bold")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def parse_fiber_csv(filepath):
    df = pd.read_csv(filepath)
    # 3 rows: 2A, 2X, 2B
    cols = df.columns
    trt_idx = [i for i, c in enumerate(cols) if any(k in c.lower() for k in ['dexa', 'llc', 'old'])]
    idx = trt_idx[0] if trt_idx else len(cols) // 2

    res = {}
    for r in range(3):
        fname = str(df.iloc[r, 0])
        if pd.isna(fname) or fname.strip() == '':
            fname = ['2A', '2X', '2B'][r]
        c_vals = pd.to_numeric(df.iloc[r, 1:idx], errors='coerce').dropna().values
        t_vals = pd.to_numeric(df.iloc[r, idx:], errors='coerce').dropna().values
        res[fname] = (c_vals, t_vals)
    return res

def draw_fiber_grouped(ax, fiber_dict, title, ylab, trt_name, trt_color, ylim_100=False):
    fiber_names = list(fiber_dict.keys())
    x = np.arange(len(fiber_names))
    width = 0.35

    p_vals = []
    max_val = 0

    for i, fn in enumerate(fiber_names):
        c, t = fiber_dict[fn]
        m_c, m_t = np.mean(c), np.mean(t)
        se_c = np.std(c, ddof=1) / np.sqrt(len(c)) if len(c) > 1 else 0
        se_t = np.std(t, ddof=1) / np.sqrt(len(t)) if len(t) > 1 else 0

        p = ttest_ind(c, t).pvalue
        p_vals.append(p)

        ax.bar(x[i] - width/2, m_c, width=width*0.9, color="white", edgecolor="black", linewidth=1.2, yerr=se_c, capsize=3)
        ax.bar(x[i] + width/2, m_t, width=width*0.9, color="white", edgecolor=trt_color, linewidth=1.2, yerr=se_t, capsize=3, ecolor=trt_color)

        np.random.seed(42 + i)
        ax.scatter(np.random.normal(x[i] - width/2, 0.04, len(c)), c, color="grey", edgecolors="black", s=15, alpha=0.7)
        ax.scatter(np.random.normal(x[i] + width/2, 0.04, len(t)), t, color=trt_color, edgecolors=trt_color, s=15, alpha=0.8)

        local_max = max(np.max(c) if len(c)>0 else 0, np.max(t) if len(t)>0 else 0, m_c + se_c, m_t + se_t)
        if local_max > max_val:
            max_val = local_max

    y_max = 100 if ylim_100 else max_val * 1.25
    ax.set_ylim(0, y_max)
    ax.set_xticks(x)
    ax.set_xticklabels(fiber_names, fontweight="bold")
    ax.set_title(title, fontweight="bold", fontsize=10)
    ax.set_ylabel(ylab, fontsize=9, fontweight="bold")

    p_str = f"P=  {p_vals[0]:.3f}  {p_vals[1]:.3f}  {p_vals[2]:.3f}"
    ax.set_xlabel(p_str, fontsize=8, fontweight="bold")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_condition(cond_id, title_main, trt_name, trt_color, supp_dir, out_path):
    files_map = {
        "Dexa": {
            "body": "Fig1A_dexabody.csv", "tibia": "Fig1A_dexatibia.csv",
            "ta": "Fig1A_dexaTA.csv", "ta_body": "Fig1A_dexaTA_body.csv", "ta_tibia": "Fig1A_dexaTA_tibia.csv",
            "all_fnum": "Fig1A_dexa_all_fiber_number.csv", "all_fsize": "Fig1A_dexa_all_fiber_size.csv",
            "fnum": "Fig1A_dexa_fiber_number.csv", "ftype": "Fig1A_dexa_fiber_type.csv", "fsize": "Fig1A_dexa_fiber_size.csv"
        },
        "LLC": {
            "body": "Fig1B_LLC_body.csv", "tibia": "Fig1B_LLC_tibia.csv",
            "ta": "Fig1B_LLC_TA.csv", "ta_body": "Fig1B_LLC_TA_body.csv", "ta_tibia": "Fig1B_LLC_TA_tibia.csv",
            "all_fnum": "Fig1B_LLC_total_fiber.csv", "all_fsize": "Fig1B_LLC_total_size.csv",
            "fnum": "Fig1B_LLC_fiber_number.csv", "ftype": "Fig1B_LLC_fiber_type.csv", "fsize": "Fig1B_LLC_fiber_size.csv"
        },
        "Aging": {
            "body": "Fig1C_old_body.csv", "tibia": "Fig1C_old_tibia.csv",
            "ta": "Fig1C_old_TA.csv", "ta_body": "Fig1C_old_TA_body.csv", "ta_tibia": "Fig1C_old_TA_tibia.csv",
            "all_fnum": "Fig1C_old_all_fiber_number.csv", "all_fsize": "Fig1C_old_all_fiber_size.csv",
            "fnum": "Fig1C_old_fiber_number.csv", "ftype": "Fig1C_old_fiber_type.csv", "fsize": "Fig1C_old_fiber_size.csv"
        }
    }[cond_id]

    fig, axes = plt.subplots(2, 5, figsize=(18, 7.5), dpi=300)

    # Row 1: Morphometrics
    df_body = pd.read_csv(os.path.join(supp_dir, files_map['body']))
    draw_single_bar(axes[0, 0], df_body.iloc[:, 0], df_body.iloc[:, 1], "Body mass", "body mass (g)", trt_name, trt_color)

    df_tibia = pd.read_csv(os.path.join(supp_dir, files_map['tibia']))
    draw_single_bar(axes[0, 1], df_tibia.iloc[:, 0], df_tibia.iloc[:, 1], "Tibia length", "Tibia length (mm)", trt_name, trt_color)

    df_ta = pd.read_csv(os.path.join(supp_dir, files_map['ta']))
    draw_single_bar(axes[0, 2], df_ta.iloc[:, 0], df_ta.iloc[:, 1], "TA mass", "TA mass (mg)", trt_name, trt_color)

    df_tabody = pd.read_csv(os.path.join(supp_dir, files_map['ta_body']))
    draw_single_bar(axes[0, 3], df_tabody.iloc[:, 0], df_tabody.iloc[:, 1], "TA/body mass", "TA mass/body mass (mg/g)", trt_name, trt_color)

    df_tatibia = pd.read_csv(os.path.join(supp_dir, files_map['ta_tibia']))
    draw_single_bar(axes[0, 4], df_tatibia.iloc[:, 0], df_tatibia.iloc[:, 1], "TA/tibia length", "TA mass/tibia length (mg/mm)", trt_name, trt_color)

    # Row 2: Myofibers
    df_fnum = pd.read_csv(os.path.join(supp_dir, files_map['all_fnum']))
    draw_single_bar(axes[1, 0], df_fnum.iloc[:, 0], df_fnum.iloc[:, 1], "All myofiber number", "total myofiber number", trt_name, trt_color)

    df_fsize = pd.read_csv(os.path.join(supp_dir, files_map['all_fsize']))
    draw_single_bar(axes[1, 1], df_fsize.iloc[:, 0], df_fsize.iloc[:, 1], "All myofiber size", "Feret's minimal diameter (µm)", trt_name, trt_color)

    f_num = parse_fiber_csv(os.path.join(supp_dir, files_map['fnum']))
    draw_fiber_grouped(axes[1, 2], f_num, "Myofiber number", "number of myofibers", trt_name, trt_color)

    f_prop = parse_fiber_csv(os.path.join(supp_dir, files_map['ftype']))
    draw_fiber_grouped(axes[1, 3], f_prop, "Myofiber proportion", "% of myofibers", trt_name, trt_color, ylim_100=True)

    f_size = parse_fiber_csv(os.path.join(supp_dir, files_map['fsize']))
    draw_fiber_grouped(axes[1, 4], f_size, "Myofiber size", "Feret's minimal diameter (µm)", trt_name, trt_color)

    plt.suptitle(title_main, fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")

def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    supp_dir = os.path.join(repo_root, "data", "supplementary")
    fig_dir = os.path.join(repo_root, "python", "figures")
    os.makedirs(fig_dir, exist_ok=True)

    plot_condition("Dexa", "Figure 1A: Dexamethasone Atrophy Phenotype", "Dexa", "#e41a1c", supp_dir,
                   os.path.join(fig_dir, "Figure1A_Dexamethasone_Morphometrics.png"))

    plot_condition("LLC", "Figure 1B: LLC Cancer Cachexia Phenotype", "LLC", "#377eb8", supp_dir,
                   os.path.join(fig_dir, "Figure1B_LLC_Morphometrics.png"))

    plot_condition("Aging", "Figure 1C: Sarcopenic Aging Phenotype", "Old", "#4daf4a", supp_dir,
                   os.path.join(fig_dir, "Figure1C_Aging_Morphometrics.png"))

if __name__ == "__main__":
    main()

