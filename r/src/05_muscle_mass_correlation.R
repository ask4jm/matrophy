# ==============================================================================
# Script: 05_muscle_mass_correlation.R
# Purpose: Reproduce Figure 7B — Correlation of Protein Abundance with Muscle Mass
# Paper: Hunt et al. (Cell Reports 2021)
# ==============================================================================

cat("=================================================================\n")
cat("Step 6: Proteome Correlation with Muscle Mass (Figure 7B)\n")
cat("=================================================================\n")

# Determine repository root
repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
supp_dir <- file.path(repo_root, "data", "supplementary")
results_dir <- file.path(repo_root, "r", "results")
figures_dir <- file.path(repo_root, "r", "figures")

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

corr_csv <- file.path(supp_dir, "Table_S3_Proteomics_Aging_Muscle_Correlation.csv")

if (!file.exists(corr_csv)) {
  stop(sprintf("Aging-Muscle Correlation data not found at: %s", corr_csv))
}

cat(sprintf("Loading Aging-Muscle correlation dataset: %s\n", corr_csv))
corr_df <- read.csv(corr_csv, stringsAsFactors = FALSE)
cat(sprintf("Loaded %d proteins with muscle mass correlation metrics.\n", nrow(corr_df)))

# ------------------------------------------------------------------------------
# 1. IDENTIFY TOP CORRELATED PROTEINS
# ------------------------------------------------------------------------------
# Top negative correlations (proteins increasing as muscle mass decreases)
top_neg <- head(corr_df[order(corr_df$ave.correlation), ], 10)
# Top positive correlations (proteins decreasing as muscle mass decreases)
top_pos <- head(corr_df[order(-corr_df$ave.correlation), ], 10)

cat("\nTop 5 Negative Correlating Proteins (Inversely related to muscle mass):\n")
print(top_neg[1:5, c("GN", "raw.TA.mass.correlation", "TA.tibia.correlation", "ave.correlation")])

write.csv(rbind(top_neg, top_pos), file.path(results_dir, "Figure7B_Top_Muscle_Mass_Correlating_Proteins.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. PLOT FIGURE 7B: TOP CORRELATING PROTEINS BARPLOT
# ------------------------------------------------------------------------------
fig7b_path <- file.path(figures_dir, "Figure7B_Muscle_Mass_Correlation.png")
cat(sprintf("\nGenerating Figure 7B: %s...\n", fig7b_path))

png(fig7b_path, width = 3200, height = 2000, res = 300)
par(mar = c(6, 5, 4, 2), mfrow = c(1, 2))

# Panel 1: Top 10 Inversely Correlated with Muscle Mass (Markers of Muscle Loss)
barplot(top_neg$ave.correlation, names.arg = top_neg$GN,
        col = "#e41a1c", ylim = c(-1.0, 0),
        las = 2, font.axis = 2, cex.names = 0.95,
        ylab = "Pearson Correlation with Muscle Mass (r)",
        main = "Top Proteins Inversely Correlated\nwith Muscle Mass (Hunt et al. Fig 7B)",
        font.main = 2)
abline(h = 0, col = "black")
abline(h = -0.9, col = "gray60", lty = 2)

# Panel 2: Distribution of All Muscle Mass Correlations
hist(corr_df$ave.correlation, breaks = 30, col = "#377eb8", border = "white",
     xlab = "Average Correlation with Muscle Mass",
     ylab = "Number of Proteins",
     main = "Proteome-Wide Muscle Mass\nCorrelation Distribution (N = 1,326)",
     font.main = 2)
abline(v = 0, col = "black", lwd = 1.5)
abline(v = mean(corr_df$ave.correlation, na.rm = TRUE), col = "gold", lwd = 2, lty = 2)
legend("topleft", legend = c("Proteins", sprintf("Mean r = %.2f", mean(corr_df$ave.correlation, na.rm = TRUE))),
       fill = c("#377eb8", NA), lty = c(NA, 2), lwd = c(NA, 2), col = c(NA, "gold"),
       bty = "o", bg = "white")

dev.off()
cat(sprintf("Saved: %s\n", fig7b_path))

cat("\n05_muscle_mass_correlation.R completed successfully!\n")

