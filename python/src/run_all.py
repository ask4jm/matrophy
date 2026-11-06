"""
Master Orchestrator: run_all.py
Runs the complete Python reproduction suite for Hunt et al. (Cell Reports 2021).
Generates all differential expression, publication figures, multi-omics integration, and muscle mass correlations.
"""

import os
import sys
import subprocess

def run_step(script_name):
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), script_name)
    print(f"\n>>> Running {script_name}...")
    res = subprocess.run([sys.executable, script_path], check=True)
    if res.returncode != 0:
        print(f"Error running {script_name}!")
        sys.exit(res.returncode)

def main():
    print("=================================================================")
    print("100% Python Reproduction Suite (Hunt et al., Cell Reports 2021)")
    print("=================================================================")

    run_step("01_plot_figure1_morphometrics.py")
    run_step("differential_expression.py")
    run_step("plot_reproduction_figures.py")
    run_step("02_plot_figure2_gsea_venn.py")
    run_step("03_multiomics_integration.py")
    run_step("05_muscle_mass_correlation.py")

    print("\n=================================================================")
    print("All Python reproduction modules executed successfully!")
    print("Output figures: python/figures/")
    print("Output results: python/results/")
    print("=================================================================")

if __name__ == "__main__":
    main()

