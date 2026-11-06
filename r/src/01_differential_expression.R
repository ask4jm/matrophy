# ==============================================================================
# Script: 01_differential_expression.R
# Purpose: Reproduce differential expression analysis of skeletal muscle atrophy
# Paper: Hunt et al. (Cell Reports 2021)
# ==============================================================================

cat("=========================================================\n")
cat("Starting R Differential Expression Analysis (Hunt et al. 2021)\n")
cat("=========================================================\n")

# Determine repository root
initial_dir <- getwd()
# Find data directory
if (file.exists("data/GSE159952_all.tpm.txt.gz")) {
  repo_root <- getwd()
} else if (file.exists("../../data/GSE159952_all.tpm.txt.gz")) {
  repo_root <- normalizePath("../..")
} else if (file.exists("../data/GSE159952_all.tpm.txt.gz")) {
  repo_root <- normalizePath("..")
} else {
  repo_root <- "C:/ME/matrophy"
}

data_path <- file.path(repo_root, "data", "GSE159952_all.tpm.txt.gz")
results_dir <- file.path(repo_root, "r", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cat(sprintf("Loading expression matrix: %s\n", data_path))
gz_con <- gzfile(data_path, "rt")
tpm_data <- read.delim(gz_con, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
close(gz_con)

cat(sprintf("Total genes loaded: %d across %d samples\n", nrow(tpm_data), ncol(tpm_data) - 1))

# Filter unexpressed genes (sum of TPM > 0)
expr_cols <- setdiff(colnames(tpm_data), "ID")
row_sums <- rowSums(tpm_data[, expr_cols])
filtered_data <- tpm_data[row_sums > 0, ]
cat(sprintf("Genes with active expression: %d\n", nrow(filtered_data)))

# Vectorized two-sample t-test function
fast_ttest_de <- function(mat, test_cols, ctrl_cols, condition_name) {
  cat(sprintf("\nAnalyzing %s (%d test vs %d ctrl samples)...\n", 
              condition_name, length(test_cols), length(ctrl_cols)))
  
  X1 <- as.matrix(mat[, test_cols])
  X2 <- as.matrix(mat[, ctrl_cols])
  n1 <- ncol(X1)
  n2 <- ncol(X2)
  
  m1 <- rowMeans(X1)
  m2 <- rowMeans(X2)
  
  # Sample variances
  v1 <- rowSums((X1 - m1)^2) / (n1 - 1)
  v2 <- rowSums((X2 - m2)^2) / (n2 - 1)
  
  # Pooled variance
  df_val <- n1 + n2 - 2
  sp2 <- ((n1 - 1) * v1 + (n2 - 1) * v2) / df_val
  se <- sqrt(sp2 * (1 / n1 + 1 / n2))
  
  # t-statistic and two-tailed p-value
  diff_mean <- m1 - m2
  t_stat <- ifelse(se > 0, diff_mean / se, 0)
  p_val <- 2 * pt(-abs(t_stat), df = df_val)
  p_val[is.na(p_val)] <- 1.0
  p_val[se == 0] <- 1.0
  
  # Log2 fold change: log2(mean_test + 1e-4) - log2(mean_ctrl + 1e-4)
  log2fc <- log2(m1 + 1e-4) - log2(m2 + 1e-4)
  
  res_df <- data.frame(
    gene = mat$ID,
    mean_test = m1,
    mean_ctrl = m2,
    log2FC = log2fc,
    p_value = p_val,
    minus_log10_p = -log10(pmax(p_val, 1e-300)),
    is_up = (p_val < 0.05) & (log2fc > 1.0),
    is_down = (p_val < 0.05) & (log2fc < -1.0),
    stringsAsFactors = FALSE
  )
  res_df$is_sig <- res_df$is_up | res_df$is_down
  
  n_up <- sum(res_df$is_up)
  n_down <- sum(res_df$is_down)
  cat(sprintf("  %s Upregulated (p < 0.05, log2FC > 1): %d\n", condition_name, n_up))
  cat(sprintf("  %s Downregulated (p < 0.05, log2FC < -1): %d\n", condition_name, n_down))
  cat(sprintf("  Total DEGs: %d\n", n_up + n_down))
  
  out_csv <- file.path(results_dir, paste0(condition_name, "_DEGs.csv"))
  write.csv(res_df, out_csv, row.names = FALSE)
  return(res_df)
}

# 1. Dexamethasone vs Control
dex_res <- fast_ttest_de(
  filtered_data,
  test_cols = paste0("dexamethasone_", 1:4),
  ctrl_cols = paste0("Control_dexamethasone_", 1:4),
  condition_name = "Dexamethasone"
)

# 2. LLC Cancer Cachexia vs Control
llc_res <- fast_ttest_de(
  filtered_data,
  test_cols = paste0("LLC_", 1:4),
  ctrl_cols = paste0("Control_LLC_", 1:4),
  condition_name = "LLC_Cancer"
)

# 3. Aging (24m vs 6m)
age_res <- fast_ttest_de(
  filtered_data,
  test_cols = paste0("24month_", 1:3),
  ctrl_cols = paste0("6month_", 1:3),
  condition_name = "Aging"
)

# 4. siCYR61 Knockdown vs siNT
cyr_res <- fast_ttest_de(
  filtered_data,
  test_cols = paste0("SICYR61_", 1:4),
  ctrl_cols = paste0("SINT_", 1:4),
  condition_name = "siCYR61"
)

# -------------------------------------------------------------
# Overlap Analysis (Figure 2D)
# -------------------------------------------------------------
cat("\nComputing Overlap and Stimulus-Specificity...\n")
up_dex <- dex_res$gene[dex_res$is_up]
up_llc <- llc_res$gene[llc_res$is_up]
up_age <- age_res$gene[age_res$is_up]

down_dex <- dex_res$gene[dex_res$is_down]
down_llc <- llc_res$gene[llc_res$is_down]
down_age <- age_res$gene[age_res$is_down]

all_3_up <- intersect(intersect(up_dex, up_llc), up_age)
all_3_down <- intersect(intersect(down_dex, down_llc), down_age)

overlap_df <- data.frame(
  Category = c("Upregulated (p<0.05, log2FC>1)", "Downregulated (p<0.05, log2FC<-1)"),
  Total_Dex = c(length(up_dex), length(down_dex)),
  Total_LLC = c(length(up_llc), length(down_llc)),
  Total_Aging = c(length(up_age), length(down_age)),
  Unique_Dex = c(length(setdiff(up_dex, union(up_llc, up_age))),
                 length(setdiff(down_dex, union(down_llc, down_age)))),
  Unique_LLC = c(length(setdiff(up_llc, union(up_dex, up_age))),
                 length(setdiff(down_llc, union(down_dex, down_age)))),
  Unique_Aging = c(length(setdiff(up_age, union(up_dex, up_llc))),
                   length(setdiff(down_age, union(down_dex, down_llc)))),
  Shared_All_3 = c(length(all_3_up), length(all_3_down))
)

write.csv(overlap_df, file.path(results_dir, "Figure2D_Overlap_Summary.csv"), row.names = FALSE)
cat("Figure 2D Overlap Table:\n")
print(overlap_df)

writeLines(all_3_up, file.path(results_dir, "shared_all_3_upregulated_genes.txt"))
writeLines(all_3_down, file.path(results_dir, "shared_all_3_downregulated_genes.txt"))

# -------------------------------------------------------------
# Classic Atrogenes Analysis (Figure 2F)
# -------------------------------------------------------------
cat("\nExtracting Classic Atrogenes Profile (Figure 2F)...\n")
classic_atrogenes <- c(
  "Trim63", "Fbxo32", "Psmd11", "Fbxo30", "Fbxo21", "Foxo1", "Foxo3",
  "Ctsl", "Cdkn1a", "Gadd45a", "Mstn", "Casp3", "Ubr5", "S100a9", "Il1b", "Cyr61"
)

atrogene_rows <- lapply(classic_atrogenes, function(g) {
  d_row <- dex_res[dex_res$gene == g, ]
  l_row <- llc_res[llc_res$gene == g, ]
  a_row <- age_res[age_res$gene == g, ]
  
  data.frame(
    Gene = g,
    Dex_log2FC = if (nrow(d_row) > 0) d_row$log2FC[1] else 0,
    Dex_pval = if (nrow(d_row) > 0) d_row$p_value[1] else 1,
    LLC_log2FC = if (nrow(l_row) > 0) l_row$log2FC[1] else 0,
    LLC_pval = if (nrow(l_row) > 0) l_row$p_value[1] else 1,
    Aging_log2FC = if (nrow(a_row) > 0) a_row$log2FC[1] else 0,
    Aging_pval = if (nrow(a_row) > 0) a_row$p_value[1] else 1,
    stringsAsFactors = FALSE
  )
})
atrogene_df <- do.call(rbind, atrogene_rows)
write.csv(atrogene_df, file.path(results_dir, "Figure2F_Classic_Atrogenes.csv"), row.names = FALSE)
print(atrogene_df[, c("Gene", "Dex_log2FC", "LLC_log2FC", "Aging_log2FC")])

# -------------------------------------------------------------
# Cyr61 Wnt Pathway Targets (Figure 6H)
# -------------------------------------------------------------
cat("\nExtracting Cyr61 Wnt Pathway Targets (Figure 6H)...\n")
wnt_targets <- c("Wnt6", "Wnt7b", "Sfrp5", "Rspo1", "Kremen2", "Wnt4", "Fzd1", "Fzd4", "Dkk3", "Axin2", "Ctnnb1")

wnt_rows <- lapply(wnt_targets, function(g) {
  r <- cyr_res[cyr_res$gene == g, ]
  if (nrow(r) > 0) {
    data.frame(
      Gene = g,
      log2FC_siCYR61 = r$log2FC[1],
      pval_siCYR61 = r$p_value[1],
      minus_log10_p = r$minus_log10_p[1],
      is_sig = r$is_sig[1],
      stringsAsFactors = FALSE
    )
  } else NULL
})
wnt_df <- do.call(rbind, wnt_rows)
write.csv(wnt_df, file.path(results_dir, "Figure6_Wnt_Pathways_siCYR61.csv"), row.names = FALSE)
print(wnt_df)

cat("\n01_differential_expression.R completed successfully!\n")
