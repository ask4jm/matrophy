# ==============================================================================
# Script: 02_plot_figures.R
# Purpose: Generate publication-quality figures reproducing Hunt et al. (Cell Reports 2021) in Base R
# Dependencies: None (Uses standard base R graphics and stats)
# ==============================================================================

cat("=========================================================\n")
cat("Generating Reproduction Figures with Base R Graphics...\n")
cat("=========================================================\n")

# Path resolution
if (file.exists("r/results/Dexamethasone_DEGs.csv")) {
  repo_root <- getwd()
} else if (file.exists("../results/Dexamethasone_DEGs.csv")) {
  repo_root <- normalizePath("..")
} else if (file.exists("../../r/results/Dexamethasone_DEGs.csv")) {
  repo_root <- normalizePath("../..")
} else {
  repo_root <- "C:/ME/matrophy"
}

results_dir <- file.path(repo_root, "r", "results")
figures_dir <- file.path(repo_root, "r", "figures")
data_path <- file.path(repo_root, "data", "GSE159952_all.tpm.txt.gz")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Load differential expression results
dex_df <- read.csv(file.path(results_dir, "Dexamethasone_DEGs.csv"), stringsAsFactors = FALSE)
llc_df <- read.csv(file.path(results_dir, "LLC_Cancer_DEGs.csv"), stringsAsFactors = FALSE)
age_df <- read.csv(file.path(results_dir, "Aging_DEGs.csv"), stringsAsFactors = FALSE)
cyr_df <- read.csv(file.path(results_dir, "siCYR61_DEGs.csv"), stringsAsFactors = FALSE)

# -------------------------------------------------------------
# 1. FIGURE 2A-C: VOLCANO PLOTS
# -------------------------------------------------------------
cat("Plotting Figure 2 Volcano Plots...\n")
fig2_path <- file.path(figures_dir, "Figure2_VolcanoPlots.png")
png(fig2_path, width = 3600, height = 1200, res = 300)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1.5), oma = c(0, 0, 1, 0))

draw_volcano <- function(df, title_str) {
  capped_p <- pmin(df$minus_log10_p, 15)
  
  plot(df$log2FC, capped_p,
       pch = 16, cex = 0.5, col = "gray80",
       xlim = c(-8, 8), ylim = c(0, 16),
       xlab = "log2(Fold Change)", ylab = "-log10(p-value)",
       main = title_str, font.main = 2, cex.lab = 1.1, cex.axis = 1.0)
  
  grid(col = "gray92", lty = 1)
  
  # Non-significant
  points(df$log2FC[!df$is_sig], capped_p[!df$is_sig], pch = 16, cex = 0.5, col = "gray80")
  # Upregulated (orange/red)
  points(df$log2FC[df$is_up], capped_p[df$is_up], pch = 16, cex = 0.7, col = "#d95f02")
  # Downregulated (blue)
  points(df$log2FC[df$is_down], capped_p[df$is_down], pch = 16, cex = 0.7, col = "#1f78b4")
  
  # Threshold lines
  abline(h = -log10(0.05), lty = 2, col = "gray40", lwd = 1)
  abline(v = c(-1.0, 1.0), lty = 2, col = "gray40", lwd = 1)
  
  # Annotate classic atrogenes and markers
  key_genes <- c("Trim63", "Fbxo32", "Cyr61", "Psmd11", "Foxo1", "Foxo3", "S100a9", "Il1b", "Cdkn1a", "Ctsl")
  sub_keys <- df[df$gene %in% key_genes & df$is_sig, ]
  if (nrow(sub_keys) > 0) {
    points(sub_keys$log2FC, pmin(sub_keys$minus_log10_p, 15), pch = 16, cex = 1.3, col = "black")
    text(sub_keys$log2FC, pmin(sub_keys$minus_log10_p, 15), labels = sub_keys$gene,
         pos = 4, offset = 0.3, cex = 0.75, font = 2)
  }
  
  legend("topleft",
         legend = c(sprintf("Up (%d)", sum(df$is_up)),
                    sprintf("Down (%d)", sum(df$is_down)),
                    "Not Sig"),
         col = c("#d95f02", "#1f78b4", "gray70"),
         pch = 16, bty = "o", bg = "white", box.col = "gray80", cex = 0.85)
}

draw_volcano(dex_df, "Dexamethasone vs Control")
draw_volcano(llc_df, "LLC Cancer Cachexia vs Control")
draw_volcano(age_df, "Aging (24m vs 6m)")

dev.off()
cat(sprintf("Saved: %s\n", fig2_path))

# -------------------------------------------------------------
# 2. FIGURE 2D: STIMULUS-SPECIFIC VS SHARED DEGs
# -------------------------------------------------------------
cat("Plotting Figure 2D Overlap Bar Plot...\n")
fig2d_path <- file.path(figures_dir, "Figure2D_Stimulus_Specific_vs_Shared.png")
png(fig2d_path, width = 2400, height = 1500, res = 300)

par(mar = c(6, 4.5, 3, 1))

up_dex <- dex_df$gene[dex_df$is_up]
up_llc <- llc_df$gene[llc_df$is_up]
up_age <- age_df$gene[age_df$is_up]

down_dex <- dex_df$gene[dex_df$is_down]
down_llc <- llc_df$gene[llc_df$is_down]
down_age <- age_df$gene[age_df$is_down]

all_3_up <- intersect(intersect(up_dex, up_llc), up_age)
all_3_down <- intersect(intersect(down_dex, down_llc), down_age)

calc_shares <- function(s_main, s_other1, s_other2, all3) {
  unique_main <- setdiff(s_main, union(s_other1, s_other2))
  shared2 <- setdiff(intersect(s_main, union(s_other1, s_other2)), all3)
  return(c(Unique = length(unique_main), Shared_2 = length(shared2), Shared_3 = length(all3)))
}

mat_shares <- cbind(
  calc_shares(up_dex, up_llc, up_age, all_3_up),
  calc_shares(up_llc, up_dex, up_age, all_3_up),
  calc_shares(up_age, up_dex, up_llc, all_3_up),
  calc_shares(down_dex, down_llc, down_age, all_3_down),
  calc_shares(down_llc, down_dex, down_age, all_3_down),
  calc_shares(down_age, down_dex, down_llc, all_3_down)
)
colnames(mat_shares) <- c("Dex Up", "LLC Up", "Aging Up", "Dex Down", "LLC Down", "Aging Down")

bp <- barplot(mat_shares, beside = FALSE,
              col = c("#4daf4a", "#377eb8", "#e41a1c"),
              ylim = c(0, max(colSums(mat_shares)) * 1.15),
              las = 2, font.axis = 2, cex.names = 0.9,
              ylab = "Number of DEGs (p < 0.05, |log2FC| > 1)",
              main = "Figure 2D Reproduction: Stimulus-Specific vs Shared Transcriptional Programs",
              font.main = 2)

legend("topright",
       legend = c("Stimulus-Specific (Unique)", "Shared by 2 Stimuli", "Shared by All 3 Stimuli"),
       fill = c("#4daf4a", "#377eb8", "#e41a1c"),
       bty = "o", bg = "white", box.col = "gray80", cex = 0.85)

# Add % unique label above each bar
for (i in 1:ncol(mat_shares)) {
  tot <- sum(mat_shares[, i])
  uniq <- mat_shares[1, i]
  pct <- round((uniq / tot) * 100, 1)
  text(bp[i], tot + 30, labels = paste0(pct, "%\nUnique"), cex = 0.75, font = 2, col = "#2b6329")
}

dev.off()
cat(sprintf("Saved: %s\n", fig2d_path))

# -------------------------------------------------------------
# 3. FIGURE 2F: CLASSIC ATROGENES
# -------------------------------------------------------------
cat("Plotting Figure 2F Classic Atrogenes Bar Plot...\n")
fig2f_path <- file.path(figures_dir, "Figure2F_Classic_Atrogenes.png")
png(fig2f_path, width = 3000, height = 1500, res = 300)

par(mar = c(5, 4.5, 3, 1))

atrogene_df <- read.csv(file.path(results_dir, "Figure2F_Classic_Atrogenes.csv"), stringsAsFactors = FALSE)
fc_mat <- rbind(
  Dex = atrogene_df$Dex_log2FC,
  LLC = atrogene_df$LLC_log2FC,
  Aging = atrogene_df$Aging_log2FC
)
colnames(fc_mat) <- atrogene_df$Gene

bp_genes <- barplot(fc_mat, beside = TRUE,
                    col = c("#e41a1c", "#377eb8", "#4daf4a"),
                    ylim = c(-1.5, 3.0),
                    las = 2, font.axis = 2, cex.names = 0.9,
                    ylab = "log2(Fold Change)",
                    main = "Figure 2F Reproduction: Expression of Classic Atrogenes Across Distinct Atrophy Stimuli",
                    font.main = 2)

abline(h = 0, col = "black", lwd = 1)
abline(h = c(-1, 1), col = "gray60", lty = 3)

legend("topright",
       legend = c("Dexamethasone", "LLC Cancer Cachexia", "Aging (24m vs 6m)"),
       fill = c("#e41a1c", "#377eb8", "#4daf4a"),
       bty = "o", bg = "white", box.col = "gray80", cex = 0.85)

# Add significance asterisks
for (i in 1:nrow(atrogene_df)) {
  # Dex
  if (atrogene_df$Dex_pval[i] < 0.05) {
    y_val <- atrogene_df$Dex_log2FC[i] + ifelse(atrogene_df$Dex_log2FC[i] >= 0, 0.1, -0.2)
    text(bp_genes[1, i], y_val, "*", font = 2, cex = 1.1)
  }
  # LLC
  if (atrogene_df$LLC_pval[i] < 0.05) {
    y_val <- atrogene_df$LLC_log2FC[i] + ifelse(atrogene_df$LLC_log2FC[i] >= 0, 0.1, -0.2)
    text(bp_genes[2, i], y_val, "*", font = 2, cex = 1.1)
  }
  # Aging
  if (atrogene_df$Aging_pval[i] < 0.05) {
    y_val <- atrogene_df$Aging_log2FC[i] + ifelse(atrogene_df$Aging_log2FC[i] >= 0, 0.1, -0.2)
    text(bp_genes[3, i], y_val, "*", font = 2, cex = 1.1)
  }
}

dev.off()
cat(sprintf("Saved: %s\n", fig2f_path))

# -------------------------------------------------------------
# 4. FIGURE S1C: PRINCIPAL COMPONENT ANALYSIS (PCA)
# -------------------------------------------------------------
cat("Plotting Figure S1C Principal Component Analysis...\n")
figs1c_path <- file.path(figures_dir, "FigureS1C_PCA.png")
png(figs1c_path, width = 2400, height = 2000, res = 300)

gz_con <- gzfile(data_path, "rt")
tpm_full <- read.delim(gz_con, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
close(gz_con)

sample_cols <- setdiff(colnames(tpm_full), "ID")
log_expr <- log2(as.matrix(tpm_full[, sample_cols]) + 1.0)

# Filter top 5000 variable genes
gene_vars <- apply(log_expr, 1, var)
top5k <- order(gene_vars, decreasing = TRUE)[1:5000]
expr_top <- log_expr[top5k, ]

# Run PCA
pca_res <- prcomp(t(expr_top), center = TRUE, scale. = FALSE)
var_pct <- summary(pca_res)$importance[2, 1:2] * 100

pca_df <- data.frame(
  Sample = sample_cols,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
)

pca_df$Group <- "Other"
pca_df$Group[grepl("24month", pca_df$Sample)] <- "Aging (24-Month)"
pca_df$Group[grepl("6month", pca_df$Sample)] <- "Young (6-Month)"
pca_df$Group[grepl("^dexamethasone", pca_df$Sample)] <- "Dexamethasone"
pca_df$Group[grepl("Control_dexamethasone", pca_df$Sample)] <- "Control (Dex)"
pca_df$Group[grepl("^LLC", pca_df$Sample)] <- "LLC Cachexia"
pca_df$Group[grepl("Control_LLC", pca_df$Sample)] <- "Control (LLC)"
pca_df$Group[grepl("SICYR61", pca_df$Sample)] <- "siCYR61"
pca_df$Group[grepl("SINT", pca_df$Sample)] <- "siNT (Control)"

write.csv(pca_df, file.path(results_dir, "FigureS1C_PCA_Coordinates.csv"), row.names = FALSE)

group_colors <- c(
  "Aging (24-Month)" = "#e41a1c", "Young (6-Month)" = "#377eb8",
  "Dexamethasone" = "#ff7f00", "Control (Dex)" = "#fdbf6f",
  "LLC Cachexia" = "#984ea3", "Control (LLC)" = "#cab2d6",
  "siCYR61" = "#4daf4a", "siNT (Control)" = "#b2df8a"
)
group_pchs <- c(
  "Aging (24-Month)" = 16, "Young (6-Month)" = 16,
  "Dexamethasone" = 15, "Control (Dex)" = 15,
  "LLC Cachexia" = 17, "Control (LLC)" = 17,
  "siCYR61" = 18, "siNT (Control)" = 18
)

par(mar = c(4.5, 4.5, 3, 8), xpd = TRUE)
plot(pca_df$PC1, pca_df$PC2,
     col = group_colors[pca_df$Group], pch = group_pchs[pca_df$Group],
     cex = 1.8, lwd = 1.5,
     xlab = sprintf("PC1 (%.1f%% Variance)", var_pct[1]),
     ylab = sprintf("PC2 (%.1f%% Variance)", var_pct[2]),
     main = "Figure S1C Reproduction: PCA of Skeletal Muscle RNA-seq Profiles",
     font.main = 2)
grid(col = "gray90", lty = 1)
points(pca_df$PC1, pca_df$PC2,
       col = group_colors[pca_df$Group], pch = group_pchs[pca_df$Group],
       cex = 1.8, lwd = 1.5)

legend("topright", inset = c(-0.35, 0),
       legend = names(group_colors),
       col = group_colors, pch = group_pchs,
       pt.cex = 1.3, bty = "o", bg = "white", box.col = "gray80", cex = 0.75)

dev.off()
cat(sprintf("Saved: %s\n", figs1c_path))

# -------------------------------------------------------------
# 5. FIGURE 6F, 6H: CYR61 / WNT PATHWAYS
# -------------------------------------------------------------
cat("Plotting Figure 6 Cyr61 Knockdown and Wnt Pathway...\n")
fig6_path <- file.path(figures_dir, "Figure6_Cyr61_Knockdown_Wnt.png")
png(fig6_path, width = 3000, height = 1300, res = 300)

par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 2))

# 6F: Volcano
draw_volcano(cyr_df, "Figure 6F: siCYR61 vs siNT in Aged Muscle")

# 6H: Wnt barplot
wnt_df <- read.csv(file.path(results_dir, "Figure6_Wnt_Pathways_siCYR61.csv"), stringsAsFactors = FALSE)
bp_wnt <- barplot(wnt_df$log2FC_siCYR61, names.arg = wnt_df$Gene,
                  col = "#7570b3", ylim = c(-0.5, 1.5), las = 2,
                  font.axis = 2, cex.names = 0.85,
                  ylab = "log2(Fold Change: siCYR61 / siNT)",
                  main = "Figure 6H: Modulation of Wnt Pathway Components",
                  font.main = 2)
abline(h = 0, col = "black")
abline(h = c(-1, 1), col = "gray60", lty = 3)

for (i in 1:nrow(wnt_df)) {
  if (wnt_df$pval_siCYR61[i] < 0.05) {
    y_val <- wnt_df$log2FC_siCYR61[i] + ifelse(wnt_df$log2FC_siCYR61[i] >= 0, 0.08, -0.15)
    text(bp_wnt[i], y_val, "*", font = 2, cex = 1.2)
  }
}

dev.off()
cat(sprintf("Saved: %s\n", fig6_path))

cat("\nAll R reproduction figures generated successfully in r/figures/!\n")
