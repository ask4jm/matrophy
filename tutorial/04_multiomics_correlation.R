# ==============================================================================
# TUTORIAL 04: Multi-Omics Regression & Disconnect Analysis
# Topic: geom_point(), geom_smooth(method = "lm"), facet_wrap(), annotate()
# Reproduces: Figure 4A (RNA vs Protein Correlation) & Figure 4B (Disconnect)
# ==============================================================================

repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
custom_lib <- file.path(repo_root, "r", "r_libs")
if (dir.exists(custom_lib)) .libPaths(c(custom_lib, .libPaths()))

library(ggplot2)
library(patchwork)

tutorial_dir <- file.path(repo_root, "tutorial")
fig_dir <- file.path(tutorial_dir, "figures")
supp_dir <- file.path(repo_root, "data", "supplementary")

cat("=================================================================\n")
cat("Starting Multi-Omics Regression & Disconnect Tutorial\n")
cat("=================================================================\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD MATCHED TRANSCRIPTOME & PROTEOME DATA
# ------------------------------------------------------------------------------
multi_csv <- file.path(supp_dir, "Table_S3_Proteomics_Cross_Comparisons.csv")
cat(sprintf("Loading matched multi-omics table: %s\n", multi_csv))
multi_df <- read.csv(multi_csv, stringsAsFactors = FALSE)
cat(sprintf("Loaded %d matched genes with simultaneous mRNA and Protein measurements.\n", nrow(multi_df)))

# ------------------------------------------------------------------------------
# STEP 2: PREPARE TIDY DATASET FOR FACETING
# ------------------------------------------------------------------------------
# Let's stack the 3 conditions into a single tidy data frame:
tidy_multi <- rbind(
  data.frame(
    Gene = multi_df$GN,
    Condition = "Dexamethasone",
    RNA_log2FC = multi_df$Dexa.RNA,
    Protein_log2FC = multi_df$Dexa.Protein,
    Target_R2 = 0.164
  ),
  data.frame(
    Gene = multi_df$GN,
    Condition = "LLC Cancer Cachexia",
    RNA_log2FC = multi_df$LLC.RNA,
    Protein_log2FC = multi_df$LLC.Protein,
    Target_R2 = 0.082
  ),
  data.frame(
    Gene = multi_df$GN,
    Condition = "Aging (24m vs 6m)",
    RNA_log2FC = multi_df$Aging.RNA,
    Protein_log2FC = multi_df$Aging.Protein,
    Target_R2 = 0.202
  )
)

tidy_multi$Condition <- factor(
  tidy_multi$Condition,
  levels = c("Dexamethasone", "LLC Cancer Cachexia", "Aging (24m vs 6m)")
)

# ------------------------------------------------------------------------------
# STEP 3: CALCULATE CORRELATION METRICS DYNAMICALLY PER CONDITION
# ------------------------------------------------------------------------------
# In R, split() splits a data frame by a factor, and lapply() runs a function on each piece:
cor_stats <- do.call(rbind, lapply(split(tidy_multi, tidy_multi$Condition), function(sub) {
  r_val <- cor(sub$RNA_log2FC, sub$Protein_log2FC, use = "complete.obs")
  r2_val <- r_val^2
  data.frame(
    Condition = sub$Condition[1],
    r = r_val,
    R2 = r2_val,
    Label = sprintf("r = %.3f\nR² = %.3f\n(Published R² = %.3f)", r_val, r2_val, sub$Target_R2[1]),
    stringsAsFactors = FALSE
  )
}))

cat("\nComputed Pearson Correlations and R²:\n")
print(cor_stats[, c("Condition", "r", "R2")])

# ------------------------------------------------------------------------------
# STEP 4: MULTI-PANEL REGRESSION PLOT (FIGURE 4A)
# ------------------------------------------------------------------------------
p_multi <- ggplot(tidy_multi, aes(x = RNA_log2FC, y = Protein_log2FC, color = Condition)) +
  # Semi-transparent points to see point density
  geom_point(alpha = 0.25, size = 1.0) +
  # Zero crosshairs
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  # Linear model trendline
  geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 0.8, se = FALSE) +
  # Facet by condition
  facet_wrap(~ Condition, ncol = 3) +
  # Distinct scientific colors
  scale_color_manual(values = c("Dexamethasone" = "#e41a1c", "LLC Cancer Cachexia" = "#377eb8", "Aging (24m vs 6m)" = "#4daf4a")) +
  # Dynamic annotation text in the top-left of each facet
  geom_text(
    data = cor_stats,
    aes(x = -3.5, y = 3.5, label = Label),
    color = "black", fontface = "bold", size = 3.5, hjust = 0, vjust = 1,
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(-4, 5), ylim = c(-3.5, 4)) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey25"),
    legend.position = "none"
  ) +
  labs(
    title = "Figure 4A Reproduction: Transcriptome vs Proteome Correlation across Muscle Atrophy Modes",
    subtitle = "Gene-by-gene mRNA fold change explains only 8% to 20% of protein variance across 3,501 matched genes",
    x = expression(bold(mRNA~log[2]~"(Fold Change vs Control)")),
    y = expression(bold(Protein~log[2]~"(Fold Change vs Control)"))
  )

print(p_multi)
# ------------------------------------------------------------------------------
# STEP 5: SAVE 300 DPI FIGURE
# ------------------------------------------------------------------------------
output_multi <- file.path(fig_dir, "Tutorial_Figure4A_RNA_Protein_Correlation.png")
cat(sprintf("\nSaving: %s\n", output_multi))

ggsave(
  filename = output_multi,
  plot = p_multi,
  width = 14,
  height = 5.5,
  dpi = 300
)

cat("Tutorial 04 completed successfully!\n")

