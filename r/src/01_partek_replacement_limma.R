# ==============================================================================
# Script: 01_partek_replacement_limma.R
# Purpose: Replace Partek Genomics Suite with Bioconductor limma for Differential Expression
# Paper: Hunt et al. (Cell Reports 2021)
# ==============================================================================

cat("=================================================================\n")
cat("Step 1 & 2: Replacing Partek with Bioconductor limma\n")
cat("=================================================================\n")

# Load local libraries
if (dir.exists("r/r_libs")) {
  .libPaths(c(normalizePath("r/r_libs"), .libPaths()))
} else if (dir.exists("r_libs")) {
  .libPaths(c(normalizePath("r_libs"), .libPaths()))
} else if (dir.exists("../r_libs")) {
  .libPaths(c(normalizePath("../r_libs"), .libPaths()))
}

suppressPackageStartupMessages({
  library(limma)
})

# Path resolution
repo_root <- if (file.exists("data/GSE159952_all.tpm.txt.gz")) getwd() else normalizePath("..")
data_path <- file.path(repo_root, "data", "GSE159952_all.tpm.txt.gz")
results_dir <- file.path(repo_root, "r", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cat(sprintf("Loading TPM matrix: %s\n", data_path))
gz_con <- gzfile(data_path, "rt")
tpm_data <- read.delim(gz_con, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
close(gz_con)

cat(sprintf("Total genes: %d across %d samples\n", nrow(tpm_data), ncol(tpm_data) - 1))

# Filter unexpressed genes
expr_cols <- setdiff(colnames(tpm_data), "ID")
active_idx <- which(rowSums(tpm_data[, expr_cols]) > 0)
filtered_data <- tpm_data[active_idx, ]
cat(sprintf("Actively expressed genes (TPM sum > 0): %d\n", nrow(filtered_data)))

# Log2 transformation: log2(TPM + 1)
log_expr <- as.matrix(log2(filtered_data[, expr_cols] + 1.0))
rownames(log_expr) <- filtered_data$ID

# Build Experimental Design Matrix (Replacing Partek ANOVA model)
sample_names <- colnames(log_expr)
groups <- factor(sapply(sample_names, function(s) {
  if (grepl("24month", s)) return("Aging_24m")
  if (grepl("6month", s)) return("Young_6m")
  if (grepl("^dexamethasone", s)) return("Dexamethasone")
  if (grepl("Control_dexamethasone", s)) return("Control_Dex")
  if (grepl("^LLC", s)) return("LLC_Cancer")
  if (grepl("Control_LLC", s)) return("Control_LLC")
  if (grepl("SICYR61", s)) return("siCYR61")
  if (grepl("SINT", s)) return("siNT")
  return("Other")
}))

cat("\nExperimental Design Groups:\n")
print(table(groups))

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

# Fit linear model across all 23,193 genes
cat("\nFitting linear models with limma::lmFit()...\n")
fit <- lmFit(log_expr, design)

# Define contrast matrix corresponding to the 4 comparisons in the paper
contrast_matrix <- makeContrasts(
  Dexamethasone_vs_Control = Dexamethasone - Control_Dex,
  LLC_vs_Control           = LLC_Cancer - Control_LLC,
  Aging_vs_Young           = Aging_24m - Young_6m,
  siCYR61_vs_siNT          = siCYR61 - siNT,
  levels = design
)

cat("Computing empirical Bayes statistics with limma::eBayes()...\n")
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract results for each condition
comparisons <- c(
  "Dexamethasone" = "Dexamethasone_vs_Control",
  "LLC_Cancer"    = "LLC_vs_Control",
  "Aging"         = "Aging_vs_Young",
  "siCYR61"       = "siCYR61_vs_siNT"
)

limma_results <- list()

for (comp_name in names(comparisons)) {
  coef_name <- comparisons[[comp_name]]
  tab <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "none")
  
  res_df <- data.frame(
    gene = rownames(tab),
    log2FC = tab$logFC,
    AveExpr = tab$AveExpr,
    t_stat = tab$t,
    p_value = tab$P.Value,
    adj_p_val = tab$adj.P.Val,
    minus_log10_p = -log10(pmax(tab$P.Value, 1e-300)),
    is_up = (tab$P.Value < 0.05) & (tab$logFC > 1.0),
    is_down = (tab$P.Value < 0.05) & (tab$logFC < -1.0),
    stringsAsFactors = FALSE
  )
  res_df$is_sig <- res_df$is_up | res_df$is_down
  
  n_up <- sum(res_df$is_up)
  n_down <- sum(res_df$is_down)
  cat(sprintf("\n%s (limma empirical Bayes):\n", comp_name))
  cat(sprintf("  Upregulated (p < 0.05, log2FC > 1): %d\n", n_up))
  cat(sprintf("  Downregulated (p < 0.05, log2FC < -1): %d\n", n_down))
  cat(sprintf("  Total DEGs: %d\n", n_up + n_down))
  
  out_file <- file.path(results_dir, paste0(comp_name, "_limma_DEGs.csv"))
  write.csv(res_df, out_file, row.names = FALSE)
  limma_results[[comp_name]] <- res_df
}

cat("\nlimma execution complete. Results saved in r/results/*_limma_DEGs.csv\n")

