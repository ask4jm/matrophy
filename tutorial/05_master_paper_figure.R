# ==============================================================================
# TUTORIAL 05: Master Composite Publication Figure
# Topic: Combining Volcanoes, Barplots, and Multi-Omics into a Full Multi-Panel
# Reproduces: Cell Reports Publication-Grade Master Composite (Figures 2 & 4)
# ==============================================================================

repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
custom_lib <- file.path(repo_root, "r", "r_libs")
if (dir.exists(custom_lib)) .libPaths(c(custom_lib, .libPaths()))

library(ggplot2)
library(patchwork)

tutorial_dir <- file.path(repo_root, "tutorial")
fig_dir <- file.path(tutorial_dir, "figures")

cat("=================================================================\n")
cat("Generating Master Publication Composite Figure\n")
cat("=================================================================\n")

# Sourcing components from earlier tutorials
cat("Loading components from Tutorial 02, 03, and 04...\n")
source(file.path(tutorial_dir, "02_volcano_plots_step_by_step.R"))
source(file.path(tutorial_dir, "03_barplots_and_atrogenes.R"))
source(file.path(tutorial_dir, "04_multiomics_correlation.R"))

# ------------------------------------------------------------------------------
# STEP 1: ASSEMBLE PANELS
# ------------------------------------------------------------------------------
# Top Row: 3 Volcano plots side-by-side
row_volcano <- (p_dexa | p_llc | p_aging)

# Middle Row: Stimulus-Specificity (Left) & Classic Atrogenes (Right)
row_bars <- (p_overlap | p_atrogenes) + plot_layout(widths = c(1, 1.4))

# Bottom Row: Multi-Omics Correlation Scatter Plots
row_multi <- p_multi

# Assemble all rows using patchwork's '/' operator
master_figure <- (row_volcano / row_bars / row_multi) +
  plot_layout(heights = c(1.0, 1.1, 1.0)) +
  plot_annotation(
    title = "Integrated Transcriptomic & Proteomic Profiling of Skeletal Muscle Atrophy",
    subtitle = "Replication of Liam C. Hunt et al., Cell Reports 37(6): 109971 (2021)",
    tag_levels = 'A', # Automatically labels panels A, B, C, D... for publication!
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

print(master_figure)
# ------------------------------------------------------------------------------
# STEP 2: SAVE FINAL HIGH-RES PUBLICATION FIGURES
# ------------------------------------------------------------------------------
# 1. 300 DPI PNG (for presentations, documents, reviews)
out_png <- file.path(fig_dir, "Tutorial_Master_Paper_Ready_Figure.png")
cat(sprintf("\nSaving 300 DPI Master PNG: %s\n", out_png))
ggsave(
  filename = out_png,
  plot = master_figure,
  width = 17,
  height = 16,
  dpi = 300
)

# 2. Vector PDF (lossless, editable in Adobe Illustrator or Inkscape)
out_pdf <- file.path(fig_dir, "Tutorial_Master_Paper_Ready_Figure.pdf")
cat(sprintf("Saving Lossless Vector PDF: %s\n", out_pdf))
ggsave(
  filename = out_pdf,
  plot = master_figure,
  width = 17,
  height = 16,
  device = "pdf"
)

cat("\nMaster Paper-Ready Figure generated successfully in tutorial/figures/!\n")
cat("=================================================================\n")

