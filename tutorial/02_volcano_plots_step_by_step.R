# ==============================================================================
# TUTORIAL 02: Volcano Plots Step-by-Step with ggplot2
# Topic: The Grammar of Graphics, Layers, Thresholds, Landmark Labels & Patchwork
# Reproduces: Figure 2A, 2B, 2C (Hunt et al., Cell Reports 2021)
# ==============================================================================

# Ensure local packages (ggplot2, patchwork, scales) are loaded
repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
custom_lib <- file.path(repo_root, "r", "r_libs")
if (dir.exists(custom_lib)) .libPaths(c(custom_lib, .libPaths()))

library(ggplot2)
library(patchwork)

# Source data cleaning from Tutorial 01
tutorial_dir <- file.path(repo_root, "tutorial")
fig_dir <- file.path(tutorial_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

source(file.path(tutorial_dir, "01_r_basics_and_data_loading.R"))

cat("\n=================================================================\n")
cat("Starting Volcano Plot Tutorial\n")
cat("=================================================================\n")

# ------------------------------------------------------------------------------
# STEP 1: WHAT IS A VOLCANO PLOT?
# ------------------------------------------------------------------------------
# A volcano plot combines two statistical measures for every single gene:
# - X-axis: Magnitude of biological effect -> log2(Fold Change).
#           Positive = induced/upregulated; Negative = repressed/downregulated.
# - Y-axis: Statistical significance -> -log10(P-value).
#           Higher up = more statistically significant (p < 0.05 -> -log10 > 1.301).

# ------------------------------------------------------------------------------
# STEP 2: THE "GRAMMAR OF GRAPHICS" IN GGPLOT2
# ------------------------------------------------------------------------------
# Every ggplot is built from 3 core components:
# 1. DATA: The data frame containing our columns (e.g. clean_dexa).
# 2. AES (Aesthetics): Mapping columns to visual dimensions (x = log2FC, y = minusLog10P, color = Significance).
# 3. GEOMS (Geometric Objects): The visual marks (geom_point, geom_line, geom_bar).

# Let's create a reusable, professional volcano plotting function:
build_paper_volcano <- function(df, condition_title, up_color = "#e41a1c", down_color = "#377eb8") {
  
  # Landmark genes highlighted in Hunt et al. Figure 2:
  landmark_genes <- c("Trim63", "Fbxo32", "Foxo1", "Foxo3", "Cdkn1a", "Cyr61", "Gadd45a", "S100a9", "Il1b")
  highlight_df <- df[df$Gene %in% landmark_genes & df$Significance != "NS", ]
  
  # Layer 1: Base ggplot with background non-significant points (grey, small, semi-transparent)
  p <- ggplot(df, aes(x = log2FC, y = minusLog10P)) +
    
    # Non-significant genes: grey85, alpha=0.3 to avoid visual crowding
    geom_point(data = subset(df, Significance == "NS"),
               color = "grey85", size = 0.8, alpha = 0.35) +
    
    # Significant genes: colored by direction
    geom_point(data = subset(df, Significance != "NS"),
               aes(color = Significance), size = 1.2, alpha = 0.85) +
    
    # Manual color palette: Up = red/green, Down = blue, NS = grey
    scale_color_manual(
      values = c("Up" = up_color, "Down" = down_color),
      labels = c("Up" = "Upregulated (p<0.05, FC>2)", "Down" = "Downregulated (p<0.05, FC<0.5)")
    ) +
    
    # Layer 2: Threshold guidelines
    # Vertical dashed lines at log2FC = -1 and +1 (2-fold change)
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    # Horizontal dashed line at p = 0.05 (-log10(0.05) ≈ 1.301)
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    
    # Layer 3: Highlight and label canonical atrogenes
    geom_point(data = highlight_df,
               shape = 21, size = 2.5, fill = "yellow", color = "black", stroke = 0.8) +
    geom_text(data = highlight_df,
              aes(label = Gene),
              vjust = -0.7, hjust = 0.5, fontface = "bold.italic", size = 3.2, color = "black") +
    
    # Layer 4: Mathematical axis labels using expression()
    labs(
      title = condition_title,
      subtitle = sprintf("Total: %s genes | Up: %s | Down: %s",
                         format(nrow(df), big.mark = ","),
                         sum(df$Significance == "Up"),
                         sum(df$Significance == "Down")),
      x = expression(bold(log[2]~"(Fold Change vs Control)")),
      y = expression(bold(-log[10]~"(P-value)"))
    ) +
    
    # Layer 5: Publication styling with theme_classic()
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey30"),
      axis.title = element_text(face = "bold", size = 10.5),
      axis.text = element_text(color = "black", size = 9),
      legend.position = "none" # Keep individual panels clean
    ) +
    coord_cartesian(xlim = c(-5, 6), ylim = c(0, 16))
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 3: GENERATING INDIVIDUAL PANELS
# ------------------------------------------------------------------------------
cat("Rendering Volcano Panel A: Dexamethasone...\n")
p_dexa <- build_paper_volcano(clean_dexa, "Dexamethasone (Dex)", up_color = "#e41a1c", down_color = "#377eb8")

cat("Rendering Volcano Panel B: LLC Cancer Cachexia...\n")
p_llc <- build_paper_volcano(clean_llc, "LLC Cancer Cachexia", up_color = "#e41a1c", down_color = "#377eb8")

cat("Rendering Volcano Panel C: Sarcopenic Aging...\n")
p_aging <- build_paper_volcano(clean_aging, "Aging (24m vs 6m)", up_color = "#e41a1c", down_color = "#377eb8")

# ------------------------------------------------------------------------------
# STEP 4: COMBINING PANELS USING PATCHWORK (THE '|' OPERATOR)
# ------------------------------------------------------------------------------
# Patchwork allows combining ggplot objects using intuitive math operators:
# '|' means place side-by-side (horizontal)
# '/' means stack on top of each other (vertical)
# '+' automatically packs them into a grid

cat("Stitching 3-panel Figure 2A-C with patchwork...\n")
volcano_panel <- (p_dexa | p_llc | p_aging) +
  plot_annotation(
    title = "Figure 2A-C Reproduction: Transcriptomic Skewing Across 3 Muscle Atrophy Modes",
    subtitle = "Dexamethasone & Aging show strong upward skewing, whereas LLC Cachexia is balanced/downward skewed.",
    caption = "Source Data: Hunt et al. Cell Reports 37(6): 109971 (2021) | GSE159952",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey25"),
      plot.caption = element_text(size = 8.5, color = "grey40", face = "italic")
    )
  )
print(volcano_panel) #Send to plot viewer
# ------------------------------------------------------------------------------
# STEP 5: SAVING PUBLICATION-GRADE IMAGES WITH ggsave()
# ------------------------------------------------------------------------------
# ggsave() handles exact dimensions (width, height), units ("in", "cm", "mm"),
# and resolution (dpi = 300 for journal submission).
output_file <- file.path(fig_dir, "Tutorial_Figure2_Volcano_3Panel.png")
cat(sprintf("Saving 300 DPI publication image: %s\n", output_file))


ggsave(
  filename = output_file,
  plot = volcano_panel,
  width = 16,
  height = 5.5,
  dpi = 300
)

cat("Tutorial 02 completed successfully! Volcano plot generated in tutorial/figures/.\n")
