# ==============================================================================
# Script: run_all.R
# Purpose: Master execution script for 100% R-based Reproduction of Hunt et al. 2021
# ==============================================================================

cat("=================================================================\n")
cat("100% R Reproduction Pipeline (Hunt et al., Cell Reports 2021)\n")
cat("Replacing Proprietary Tools (Partek, Spotfire, Prism, WebGestalt)\n")
cat("=================================================================\n\n")

# Robust resolution of script directory across interactive and CLI execution
find_script_dir <- function() {
  # 1. If running via Rscript command line
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  # 2. If sourced in R session, check call stack
  for (i in rev(seq_len(sys.nframe()))) {
    frame <- sys.frame(i)
    if (!is.null(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile)))
    }
  }
  # 3. Check relative to current working directory
  if (dir.exists("r/src")) {
    return(normalizePath("r/src"))
  }
  if (file.exists("01_differential_expression.R")) {
    return(getwd())
  }
  # 4. Fallback to default project location
  return("C:/ME/matrophy/r/src")
}

script_dir <- find_script_dir()
cat(sprintf("Using R scripts directory: %s\n\n", script_dir))

cat("Module 0: Morphometrics & Histological Phenotyping (Figures 1A, 1B, 1C)...\n")
source(file.path(script_dir, "01_plot_figure1_morphometrics.R"))

cat("\nModule 1: Vectorized Student's t-test Differential Expression...\n")
source(file.path(script_dir, "01_differential_expression.R"))

cat("\nModule 2: Partek Replacement with Bioconductor limma...\n")
source(file.path(script_dir, "01_partek_replacement_limma.R"))

cat("\nModule 3: Spotfire & Prism Replacement: Publication Figures...\n")
source(file.path(script_dir, "02_plot_figures.R"))
source(file.path(script_dir, "02_plot_figure2_gsea_venn.R"))

cat("\nModule 4: WebGestalt Replacement: Functional Over-Representation Analysis...\n")
source(file.path(script_dir, "04_webgestalt_replacement_gsea.R"))

cat("\nModule 5: Multi-Omics Integration: RNA-Protein Correlation & Disconnect (Figures 3, 4, 5)...\n")
source(file.path(script_dir, "03_multiomics_integration.R"))

cat("\nModule 6: Muscle Mass Correlation: Sarcopenia Phenotype Markers (Figure 7B)...\n")
source(file.path(script_dir, "05_muscle_mass_correlation.R"))

cat("\n=================================================================\n")
cat("Complete 100% R Pipeline Execution Finished Successfully!\n")
cat("All Results saved in: r/results/\n")
cat("All Figures saved in: r/figures/\n")
cat("=================================================================\n")
