# ==============================================================================
# TUTORIAL 03: Barplots, Category Ordering, & Classic Atrogene Profiles
# Topic: geom_col(), position_dodge(), factor levels, and publication aesthetics
# Reproduces: Figure 2D (Stimulus Specificity) & Figure 2F (Classic Atrogenes)
# ==============================================================================

repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
custom_lib <- file.path(repo_root, "r", "r_libs")
if (dir.exists(custom_lib)) .libPaths(c(custom_lib, .libPaths()))

library(ggplot2)
library(patchwork)

tutorial_dir <- file.path(repo_root, "tutorial")
fig_dir <- file.path(tutorial_dir, "figures")
results_dir <- file.path(repo_root, "r", "results")

cat("=================================================================\n")
cat("Starting Barplots & Atrogene Profiling Tutorial\n")
cat("=================================================================\n")

# ------------------------------------------------------------------------------
# STEP 1: UNDERSTANDING geom_col() VS geom_bar()
# ------------------------------------------------------------------------------
# In ggplot2:
# - geom_bar() counts rows automatically (stat = "count"). Useful when you have raw observations.
# - geom_col() plots values that are ALREADY calculated in your table (stat = "identity").
#   For scientific data where you have pre-calculated fold changes or means, ALWAYS use geom_col()!

# ------------------------------------------------------------------------------
# STEP 2: FIGURE 2D REPRODUCTION (STIMULUS SPECIFICITY)
# ------------------------------------------------------------------------------
# Load the pre-computed overlap summary from r/results
overlap_csv <- file.path(results_dir, "Figure2D_Overlap_Summary.csv")
overlap_raw <- read.csv(overlap_csv, stringsAsFactors = FALSE)

# Let's reshape this table into a "tidy" format suitable for ggplot2:
# Tidy rule: Each variable is a column, each observation is a row.
overlap_tidy <- data.frame(
  Condition = rep(c("Dexamethasone", "LLC Cancer", "Aging"), 2),
  Direction = c(rep("Upregulated", 3), rep("Downregulated", 3)),
  Unique_Count = c(1191, 148, 367, 362, 351, 38),
  Shared_Count = c(1266 - 1191, 175 - 148, 430 - 367, 391 - 362, 381 - 351, 40 - 38)
)

# Reshape to long format for stacked bar charts using clean Base R rbind:
overlap_long <- rbind(
  data.frame(
    Condition = overlap_tidy$Condition,
    Direction = overlap_tidy$Direction,
    Class = "Stimulus-Specific (>85%)",
    Count = overlap_tidy$Unique_Count
  ),
  data.frame(
    Condition = overlap_tidy$Condition,
    Direction = overlap_tidy$Direction,
    Class = "Shared with other modes",
    Count = overlap_tidy$Shared_Count
  )
)

# Control display order using factor levels
overlap_long$Condition <- factor(overlap_long$Condition, levels = c("Dexamethasone", "LLC Cancer", "Aging"))
overlap_long$Class <- factor(overlap_long$Class, levels = c("Stimulus-Specific (>85%)", "Shared with other modes"))

p_overlap <- ggplot(overlap_long, aes(x = Condition, y = Count, fill = Class)) +
  geom_col(position = "stack", width = 0.65, color = "black", linewidth = 0.3) +
  facet_wrap(~ Direction, scales = "free_y") +
  scale_fill_manual(values = c("Stimulus-Specific (>85%)" = "#2b83ba", "Shared with other modes" = "#d7191c")) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(
    title = "Figure 2D Reproduction: Overwhelming Stimulus-Specificity",
    subtitle = ">85% to 94% of DEGs are unique to their inducing atrophy mode",
    x = "",
    y = "Number of Differentially Expressed Genes",
    fill = ""
  )

# ------------------------------------------------------------------------------
# STEP 3: FIGURE 2F REPRODUCTION (CLASSIC ATROGENES DISCONNECT)
# ------------------------------------------------------------------------------
# Load classic atrogenes table
atrogenes_csv <- file.path(results_dir, "Figure2F_Classic_Atrogenes.csv")
atrogenes_df <- read.csv(atrogenes_csv, stringsAsFactors = FALSE)

# Let's reshape into long format:
atrogenes_long <- rbind(
  data.frame(Gene = atrogenes_df$Gene, Condition = "Dexamethasone", log2FC = atrogenes_df$Dex_log2FC),
  data.frame(Gene = atrogenes_df$Gene, Condition = "LLC Cancer", log2FC = atrogenes_df$LLC_log2FC),
  data.frame(Gene = atrogenes_df$Gene, Condition = "Aging (Sarcopenia)", log2FC = atrogenes_df$Aging_log2FC)
)

# Important beginner tip: In ggplot2, categorical order is governed by FACTOR LEVELS.
# Without factor levels, R sorts alphabetically!
# Let's order genes by their Dexamethasone fold change:
gene_order <- atrogenes_df$Gene[order(-atrogenes_df$Dex_log2FC)]
atrogenes_long$Gene <- factor(atrogenes_long$Gene, levels = gene_order)
atrogenes_long$Condition <- factor(atrogenes_long$Condition, levels = c("Dexamethasone", "LLC Cancer", "Aging (Sarcopenia)"))

p_atrogenes <- ggplot(atrogenes_long, aes(x = Gene, y = log2FC, fill = Condition)) +
  # position_dodge(width = 0.8) places bars for different conditions side-by-side
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "black", linewidth = 0.25) +
  # Zero baseline line
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Significance threshold guideline (+1.0 log2FC = 2-fold change)
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = c("Dexamethasone" = "#e41a1c", "LLC Cancer" = "#377eb8", "Aging (Sarcopenia)" = "#4daf4a")) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 10),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey25"),
    legend.position = "top"
  ) +
  labs(
    title = "Figure 2F Reproduction: Classic Atrogenes Fail in Sarcopenic Aging",
    subtitle = "MuRF1 (Trim63) & Atrogin-1 (Fbxo32) are induced in Dex & LLC, but UNCHANGED in Aging",
    x = "Canonical Atrogene",
    y = expression(bold(log[2]~"(Fold Change vs Control)")),
    fill = "Condition:"
  )

# ------------------------------------------------------------------------------
# STEP 4: COMBINE AND EXPORT
# ------------------------------------------------------------------------------
cat("Combining Figure 2D and 2F with patchwork...\n")
combined_barplots <- p_overlap / p_atrogenes +
  plot_layout(heights = c(1, 1.2))
print(combined_barplots)

output_barplots <- file.path(fig_dir, "Tutorial_Figure2D_2F_Barplots.png")
cat(sprintf("Saving: %s\n", output_barplots))

ggsave(
  filename = output_barplots,
  plot = combined_barplots,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Tutorial 03 completed successfully!\n")
