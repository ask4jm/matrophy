# ==============================================================================
# Script: 03_multiomics_integration.R
# Purpose: Multi-Omics Integration (RNA vs Protein) and Transcriptome-Proteome Disconnect
# Paper: Hunt et al. (Cell Reports 2021) — Figures 3, 4, 5
# ==============================================================================

cat("=================================================================\n")
cat("Step 5: Multi-Omics Integration & Transcriptome-Proteome Disconnect\n")
cat("=================================================================\n")

# Determine repository root
repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
supp_dir <- file.path(repo_root, "data", "supplementary")
results_dir <- file.path(repo_root, "r", "results")
figures_dir <- file.path(repo_root, "r", "figures")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. LOAD MATCHED MULTI-OMICS DATASET (Table S3 / Cross-Comparisons)
# ------------------------------------------------------------------------------
cross_csv <- file.path(supp_dir, "Table_S3_Proteomics_Cross_Comparisons.csv")
cat(sprintf("Loading matched multi-omics matrix: %s\n", cross_csv))
multi_df <- read.csv(cross_csv, stringsAsFactors = FALSE)
cat(sprintf("Loaded %d matched genes with both mRNA and Protein measurements.\n", nrow(multi_df)))

# ------------------------------------------------------------------------------
# 2. FIGURE 4A: mRNA vs PROTEIN CORRELATION SCATTER PLOTS
# ------------------------------------------------------------------------------
cat("\nComputing RNA-Protein Correlation (Figure 4A)...\n")
fig4a_path <- file.path(figures_dir, "Figure4A_RNA_Protein_Correlation.png")
png(fig4a_path, width = 3600, height = 1200, res = 300)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1.5))

cor_results <- list()

draw_cor_plot <- function(x_vals, y_vals, title_str, target_r2) {
  valid_idx <- complete.cases(x_vals, y_vals)
  x <- x_vals[valid_idx]
  y <- y_vals[valid_idx]
  
  fit <- lm(y ~ x)
  r_val <- cor(x, y)
  r2_val <- summary(fit)$r.squared
  
  plot(x, y,
       pch = 16, cex = 0.6, col = rgb(0.2, 0.4, 0.7, 0.35),
       xlim = c(-3.5, 3.5), ylim = c(-2.5, 2.5),
       xlab = "log2(Fold Change mRNA)",
       ylab = "log2(Fold Change Protein)",
       main = title_str, font.main = 2, cex.lab = 1.1)
  grid(col = "gray90", lty = 1)
  abline(h = 0, v = 0, col = "gray50", lty = 2)
  abline(fit, col = "#e41a1c", lwd = 2)
  
  text(-3.2, 2.2, labels = sprintf("r = %.3f\nr^2 = %.3f (Target: %.3f)\nN = %d",
                                  r_val, r2_val, target_r2, length(x)),
       adj = c(0, 1), font = 2, cex = 0.95,
       col = "#990000")
  
  return(data.frame(Condition = title_str, Pearson_r = r_val, R2 = r2_val, Target_R2 = target_r2, N = length(x)))
}

cor_results[["Dexa"]] <- draw_cor_plot(multi_df$Dexa.RNA, multi_df$Dexa.Protein,
                                       "Dexamethasone (Dex)", target_r2 = 0.164)
cor_results[["LLC"]]  <- draw_cor_plot(multi_df$LLC.RNA, multi_df$LLC.Protein,
                                       "LLC Cancer Cachexia", target_r2 = 0.082)
cor_results[["Aging"]] <- draw_cor_plot(multi_df$Aging.RNA, multi_df$Aging.Protein,
                                        "Aging (24m vs 6m)", target_r2 = 0.202)

dev.off()
cat(sprintf("Saved: %s\n", fig4a_path))

cor_summary_df <- do.call(rbind, cor_results)
write.csv(cor_summary_df, file.path(results_dir, "Figure4A_Correlation_Summary.csv"), row.names = FALSE)
print(cor_summary_df)

# ------------------------------------------------------------------------------
# 3. FIGURE 4B: TRANSCRIPTOME-PROTEOME DISCONNECT ANALYSIS
# ------------------------------------------------------------------------------
cat("\nAnalyzing Transcriptome-Proteome Disconnect (Figure 4B)...\n")
# Categorize mRNAs into: Downregulated (log2FC < -0.5), Unchanged (-0.5 to 0.5), Upregulated (log2FC > 0.5)
categorize_disconnect <- function(rna_vals, prot_vals, cond_name) {
  cat_down <- prot_vals[rna_vals < -0.5]
  cat_unchanged <- prot_vals[rna_vals >= -0.5 & rna_vals <= 0.5]
  cat_up <- prot_vals[rna_vals > 0.5]
  
  data.frame(
    Condition = cond_name,
    mRNA_Bin = c("Downregulated (log2FC < -0.5)", "Unchanged (-0.5 to 0.5)", "Upregulated (log2FC > 0.5)"),
    Mean_Protein_log2FC = c(mean(cat_down, na.rm = TRUE),
                            mean(cat_unchanged, na.rm = TRUE),
                            mean(cat_up, na.rm = TRUE)),
    N_Proteins = c(length(na.omit(cat_down)),
                   length(na.omit(cat_unchanged)),
                   length(na.omit(cat_up))),
    stringsAsFactors = FALSE
  )
}

disc_dex <- categorize_disconnect(multi_df$Dexa.RNA, multi_df$Dexa.Protein, "Dexamethasone")
disc_llc <- categorize_disconnect(multi_df$LLC.RNA, multi_df$LLC.Protein, "LLC Cancer")
disc_age <- categorize_disconnect(multi_df$Aging.RNA, multi_df$Aging.Protein, "Aging")

disc_summary <- rbind(disc_dex, disc_llc, disc_age)
write.csv(disc_summary, file.path(results_dir, "Figure4B_Disconnect_Summary.csv"), row.names = FALSE)
cat("Figure 4B Disconnect Table:\n")
print(disc_summary)

# Plot Figure 4B Disconnect barplot
fig4b_path <- file.path(figures_dir, "Figure4B_Disconnect_Barplot.png")
png(fig4b_path, width = 2400, height = 1500, res = 300)
par(mar = c(5, 4.5, 3, 1))

mat_disc <- rbind(
  Dex = disc_dex$Mean_Protein_log2FC,
  LLC = disc_llc$Mean_Protein_log2FC,
  Aging = disc_age$Mean_Protein_log2FC
)
colnames(mat_disc) <- c("mRNA Down (< -0.5)", "mRNA Unchanged (-0.5 to 0.5)", "mRNA Up (> 0.5)")

bp_disc <- barplot(mat_disc, beside = TRUE,
                   col = c("#e41a1c", "#377eb8", "#4daf4a"),
                   ylim = c(-0.4, 0.4),
                   ylab = "Mean Protein log2(Fold Change)",
                   main = "Figure 4B Reproduction: Protein Abundance by mRNA Change Bins",
                   font.main = 2, cex.names = 0.9)
abline(h = 0, col = "black")
abline(h = c(-0.1, 0.1), col = "gray60", lty = 3)
legend("topleft", legend = c("Dexamethasone", "LLC Cancer Cachexia", "Aging"),
       fill = c("#e41a1c", "#377eb8", "#4daf4a"),
       bty = "o", bg = "white", box.col = "gray80", cex = 0.85)

dev.off()
cat(sprintf("Saved: %s\n", fig4b_path))

# ------------------------------------------------------------------------------
# 4. FIGURE 3F: ATROGENES DISCONNECT AT PROTEIN LEVEL
# ------------------------------------------------------------------------------
cat("\nAnalyzing Atrogenes at Protein Level (Figure 3F)...\n")
fig3f_csv <- file.path(supp_dir, "Table_S6_Atroproteins_Fig3F.csv")
if (file.exists(fig3f_csv)) {
  fig3f_df <- read.csv(fig3f_csv, stringsAsFactors = FALSE)
  write.csv(fig3f_df, file.path(results_dir, "Figure3F_Atrogenes_Protein_Level.csv"), row.names = FALSE)
  cat(sprintf("Figure 3F data loaded (%d atrogenes at protein level).\n", nrow(fig3f_df)))
}

# ------------------------------------------------------------------------------
# 5. FIGURE 5: ATROPROTEINS IDENTIFICATION & VALIDATION
# ------------------------------------------------------------------------------
cat("\nExtracting Atroproteins Profile (Figure 5)...\n")

# Identify core atroproteins: consistent directional change across all 3 modes
# in matched dataset
core_up <- subset(multi_df, Aging.Protein > 0.3 & Dexa.Protein > 0.3 & LLC.Protein > 0.3)
core_down <- subset(multi_df, Aging.Protein < -0.3 & Dexa.Protein < -0.3 & LLC.Protein < -0.3)

cat(sprintf("Identified %d universally upregulated atroproteins and %d universally downregulated atroproteins.\n",
            nrow(core_up), nrow(core_down)))

core_atroproteins <- rbind(
  data.frame(core_up, Direction = "Universally Upregulated"),
  data.frame(core_down, Direction = "Universally Downregulated")
)
write.csv(core_atroproteins, file.path(results_dir, "Figure5A_Core_Atroproteins.csv"), row.names = FALSE)

# Plot Figure 5A Reproduction: Top Consistently Regulated Atroproteins
fig5_path <- file.path(figures_dir, "Figure5A_Top_Atroproteins.png")
png(fig5_path, width = 3200, height = 1800, res = 300)
par(mar = c(7, 5, 4, 2))

# Select canonical atroproteins for display (ordered by Dexa FC)
display_order <- c("Prg4", "Apoe", "Ncam1", "Clu", "Fn1", "Ahnak2", "Apod", "Pla2g7", "Ndufb2", "Kera", "Pcolce")
plot_df <- subset(core_atroproteins, GN %in% display_order)
plot_df <- plot_df[match(display_order, plot_df$GN), ]

fc_matrix <- t(as.matrix(plot_df[, c("Aging.Protein", "Dexa.Protein", "LLC.Protein")]))
colnames(fc_matrix) <- plot_df$GN

bp <- barplot(fc_matrix, beside = TRUE,
              col = c("#4daf4a", "#e41a1c", "#377eb8"),
              ylim = c(-1.0, 2.7),
              las = 2, font.axis = 2, cex.names = 1.1,
              ylab = "Protein log2(Fold Change)",
              main = "Figure 5A Reproduction: Consistently Regulated 'Atroproteins' Across 3 Atrophy Modes",
              font.main = 2)

abline(h = 0, col = "black", lwd = 1.2)
abline(h = c(-0.3, 0.3), col = "gray50", lty = 2)

legend("topright", legend = c("Aging (24m vs 6m)", "Dexamethasone", "LLC Cancer Cachexia"),
       fill = c("#4daf4a", "#e41a1c", "#377eb8"),
       bty = "o", bg = "white", box.col = "gray80", cex = 1.0)

text(x = 12, y = 2.4, labels = "Universally Upregulated", font = 4, col = "#2b8237", cex = 1.1)
text(x = 40, y = -0.9, labels = "Universally Downregulated", font = 4, col = "#b2182b", cex = 1.1)

dev.off()
cat(sprintf("Saved: %s\n", fig5_path))

cat("\n03_multiomics_integration.R completed successfully!\n")
