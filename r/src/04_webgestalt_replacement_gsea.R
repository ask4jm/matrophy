# ==============================================================================
# Script: 04_webgestalt_replacement_gsea.R
# Purpose: Replace WebGestalt with R-based Functional Over-Representation Analysis
# Paper: Hunt et al. (Cell Reports 2021)
# ==============================================================================

cat("=================================================================\n")
cat("Step 6: Replacing WebGestalt with R Functional Enrichment\n")
cat("=================================================================\n")

# Path resolution
repo_root <- if (file.exists("r/results/Dexamethasone_DEGs.csv")) getwd() else normalizePath("..")
results_dir <- file.path(repo_root, "r", "results")
figures_dir <- file.path(repo_root, "r", "figures")

# Define key Biological Process Gene Sets highlighted in Hunt et al. Table S2 & text
pathway_gene_sets <- list(
  "Collagen & Extracellular Matrix" = c(
    "Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2", "Col6a1", "Col6a2", "Col6a3",
    "Lama2", "Lamb1", "Fn1", "Postn", "Mmp2", "Mmp9", "Mmp14", "Timp1", "Timp2", "Lox"
  ),
  "Inflammatory & Immune Response (IL1/IL6)" = c(
    "Il1b", "Il6", "Il1r1", "Nfkb1", "Nfkb2", "Rela", "Stat3", "S100a8", "S100a9",
    "Ccl2", "Cxcl1", "Cxcl2", "Tlr2", "Tlr4", "Ptgs2", "Icam1", "Vcam1"
  ),
  "Lipid & Lipoprotein Metabolic Process" = c(
    "Lpl", "Fabp3", "Fabp4", "Cd36", "Cpt1b", "Acadl", "Acadvl", "Hadha", "Hadhb",
    "Apoe", "Apod", "Pparcc1a", "Pparg", "Plin2", "Plin5", "Lipe"
  ),
  "Proteasome & Ubiquitin-Dependent Degradation" = c(
    "Trim63", "Fbxo32", "Fbxo30", "Fbxo21", "Psmd11", "Psma1", "Psmb1", "Psmc1",
    "Ube2c", "Ube2d1", "Ube4b", "Ubr5", "Ctsl", "Casp3", "Foxo1", "Foxo3"
  ),
  "Mitochondrial Respiration & Gene Expression" = c(
    "Ndufs1", "Ndufs2", "Ndufv1", "Sdha", "Sdhb", "Uqcrc1", "Cox4i1", "Cox5a",
    "Atp5a1", "Atp5b", "Tfam", "Ppargc1a", "Polg", "Mrps12", "Mrpl12"
  ),
  "Wnt Signaling & Muscle Regeneration" = c(
    "Wnt4", "Wnt6", "Wnt7b", "Sfrp5", "Rspo1", "Kremen2", "Fzd1", "Fzd4", "Dkk3",
    "Axin2", "Ctnnb1", "Myod1", "Myog", "Myf5", "Myf6", "Cyr61"
  )
)

# Load DEG lists
dex_df <- read.csv(file.path(results_dir, "Dexamethasone_DEGs.csv"), stringsAsFactors = FALSE)
llc_df <- read.csv(file.path(results_dir, "LLC_Cancer_DEGs.csv"), stringsAsFactors = FALSE)
age_df <- read.csv(file.path(results_dir, "Aging_DEGs.csv"), stringsAsFactors = FALSE)

all_genes <- dex_df$gene
n_background <- length(all_genes)

# Over-representation test (Hypergeometric test / Fisher's exact test)
test_enrichment <- function(sig_genes, gene_set, background_n) {
  k <- length(intersect(sig_genes, gene_set)) # overlap
  m <- length(gene_set)                       # size of pathway
  n <- background_n - m                       # non-pathway genes
  N <- length(sig_genes)                      # size of query
  
  p_val <- phyper(k - 1, m, n, N, lower.tail = FALSE)
  return(c(Overlap = k, Pathway_Size = m, P_Value = p_val))
}

conditions <- list(
  "Dexamethasone (Dex)" = dex_df$gene[dex_df$is_up],
  "LLC Cancer Cachexia"  = llc_df$gene[llc_df$is_up],
  "Aging (Sarcopenia)"  = age_df$gene[age_df$is_up]
)

enrich_results <- list()

for (cond_name in names(conditions)) {
  sig_genes <- conditions[[cond_name]]
  rows <- lapply(names(pathway_gene_sets), function(pname) {
    pgenes <- pathway_gene_sets[[pname]]
    res <- test_enrichment(sig_genes, pgenes, n_background)
    data.frame(
      Condition = cond_name,
      Pathway = pname,
      Overlap = res[["Overlap"]],
      Pathway_Size = res[["Pathway_Size"]],
      P_Value = res[["P_Value"]],
      Enrichment_Score = (res[["Overlap"]] / length(sig_genes)) / (res[["Pathway_Size"]] / n_background),
      stringsAsFactors = FALSE
    )
  })
  enrich_results[[cond_name]] <- do.call(rbind, rows)
}

final_enrich_df <- do.call(rbind, enrich_results)
final_enrich_df$minus_log10_p <- -log10(pmax(final_enrich_df$P_Value, 1e-20))
write.csv(final_enrich_df, file.path(results_dir, "WebGestalt_Replacement_Enrichment.csv"), row.names = FALSE)

cat("\nFunctional Enrichment Results (WebGestalt replacement):\n")
print(final_enrich_df[, c("Condition", "Pathway", "Overlap", "P_Value", "Enrichment_Score")])

# Plotting Pathway Enrichment Dot Plot (Replacing WebGestalt Online Plot)
png(file.path(figures_dir, "WebGestalt_Replacement_Enrichment_Plot.png"), width = 2800, height = 1400, res = 300)
par(mar = c(5, 16, 3, 2))

pathways <- unique(final_enrich_df$Pathway)
conds <- unique(final_enrich_df$Condition)

# Reshape into matrix of -log10(p)
p_mat <- matrix(0, nrow = length(pathways), ncol = length(conds),
                dimnames = list(pathways, conds))

for (i in 1:nrow(final_enrich_df)) {
  p_mat[final_enrich_df$Pathway[i], final_enrich_df$Condition[i]] <- final_enrich_df$minus_log10_p[i]
}

bp <- barplot(t(p_mat), beside = TRUE, horiz = TRUE,
              col = c("#e41a1c", "#377eb8", "#4daf4a"),
              xlab = "-log10(Hypergeometric p-value)",
              main = "WebGestalt Replacement: Pathway Enrichment Across Atrophy Modes",
              las = 1, cex.names = 0.8, font.main = 2)

abline(v = -log10(0.05), col = "gray40", lty = 2)
legend("bottomright", legend = conds, fill = c("#e41a1c", "#377eb8", "#4daf4a"),
       bty = "o", bg = "white", box.col = "gray80", cex = 0.8)

dev.off()
cat(sprintf("\nSaved plot: %s\n", file.path(figures_dir, "WebGestalt_Replacement_Enrichment_Plot.png")))
cat("WebGestalt replacement analysis completed successfully!\n")

