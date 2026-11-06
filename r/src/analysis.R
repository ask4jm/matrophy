# ==============================================================================
# R Script for Reproducing Hunt et al. (Cell Reports 2021) RNA-seq Analysis
# Dataset: GSE159952 (Mus musculus Tibialis Anterior muscle)
# ==============================================================================

cat("Loading GSE159952 TPM expression dataset...\n")
tpm_file <- "data/GSE159952_all.tpm.txt.gz"

if (!file.exists(tpm_file)) {
  stop("Dataset not found. Please ensure GSE159952_all.tpm.txt.gz is in the data/ directory.")
}

# Read compressed tab-delimited file
data <- read.delim(gzfile(tpm_file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d genes across %d samples.\n", nrow(data), ncol(data) - 1))

# Helper function to compute differential expression
compute_de <- function(data, test_cols, ctrl_cols, condition_name) {
  cat(sprintf("\nComputing DE for %s...\n", condition_name))
  test_mat <- as.matrix(data[, test_cols])
  ctrl_mat <- as.matrix(data[, ctrl_cols])
  
  mean_test <- rowMeans(test_mat)
  mean_ctrl <- rowMeans(ctrl_mat)
  log2fc <- log2(mean_test + 1e-4) - log2(mean_ctrl + 1e-4)
  
  # Two-tailed Student's t-test per gene
  pvals <- sapply(1:nrow(data), function(i) {
    t <- tryCatch(t.test(test_mat[i, ], ctrl_mat[i, ], var.equal = TRUE)$p.value,
                  error = function(e) 1.0)
    if (is.na(t)) 1.0 else t
  })
  
  res <- data.frame(
    gene = data$ID,
    mean_test = mean_test,
    mean_ctrl = mean_ctrl,
    log2FC = log2fc,
    p_value = pvals,
    is_up = (pvals < 0.05) & (log2fc > 1.0),
    is_down = (pvals < 0.05) & (log2fc < -1.0)
  )
  
  n_up <- sum(res$is_up)
  n_down <- sum(res$is_down)
  cat(sprintf("  %s: %d Upregulated, %d Downregulated (p < 0.05, |log2FC| > 1)\n", condition_name, n_up, n_down))
  return(res)
}

# Run for Dexamethasone
dex_res <- compute_de(
  data,
  test_cols = paste0("dexamethasone_", 1:4),
  ctrl_cols = paste0("Control_dexamethasone_", 1:4),
  condition_name = "Dexamethasone"
)

# Run for LLC Cancer Cachexia
llc_res <- compute_de(
  data,
  test_cols = paste0("LLC_", 1:4),
  ctrl_cols = paste0("Control_LLC_", 1:4),
  condition_name = "LLC_Cancer"
)

# Run for Aging
age_res <- compute_de(
  data,
  test_cols = paste0("24month_", 1:3),
  ctrl_cols = paste0("6month_", 1:3),
  condition_name = "Aging"
)

cat("\nR differential expression execution completed successfully!\n")

