# ==============================================================================
# TUTORIAL 01: R Basics & Biological Data Loading
# Topic: Variables, Vectors, DataFrames, and Loading Transcriptomic Tables
# ==============================================================================
cat("\014") # cat("\f") Clear the console
ls() # rm(list=ls()) Delete all current variables
gc() # garbage collector - free memory
# ------------------------------------------------------------------------------
# SECTION 1: R FUNDAMENTALS
# ------------------------------------------------------------------------------

# In R, the "<-" symbol is the assignment operator (stores a value in a variable)
my_sample <- "Dexa_Treated_Muscle"
num_replicates <- 4
is_control <- FALSE

cat("Sample:", my_sample, "\nNumber of Replicates:", num_replicates, "\n") # concats text + values

# A VECTOR is a 1-dimensional sequence of values of the same type.
# Created with the c() function (c stands for 'combine').
gene_names <- c("Trim63", "Fbxo32", "Foxo1", "Cyr61", "Acta1")
fold_changes <- c(1.05, 0.92, 0.55, 1.17, -0.05)
p_values <- c(0.0012, 0.0045, 0.0320, 0.0008, 0.8500)

# Vectors are indexed with square brackets []. NOTE: R is 1-indexed (starts at 1, not 0!)
cat("First gene:", gene_names[1], "\n")
cat("Top 3 genes:", gene_names[1:3], "\n")
cat(gene_names[-c(2,3)]) # Omit values

# Vectorized operations: R applies math to every element simultaneously!
# Let's calculate -log10 of the p-values in one single calculation:
neg_log10_p <- -log10(p_values)
cat("-log10(P-values):", round(neg_log10_p, 2), "\n") # Round to value

# ------------------------------------------------------------------------------
# SECTION 2: WHAT IS A DATA FRAME?
# ------------------------------------------------------------------------------
# A data.frame is a 2D table where each column is a vector (variables)
# and each row represents an observation (e.g. a gene).

toy_de_table <- data.frame(
  Gene = gene_names,
  log2FC = fold_changes,
  P_Value = p_values,
  minusLog10P = -log10(p_values),
  stringsAsFactors = FALSE
)

cat("\n--- Preview of Toy Data Frame ---\n")
print(toy_de_table)

view(toy_de_table)
View(toy_de_table)

find("print")

# Check where View comes from:
find("View")
# [1] "package:utils"

# Check where view comes from (after loading tidyverse):
find("view")
# [1] "package:tibble"

# Accessing columns with the $ operator:
cat("\nJust the fold changes: ", toy_de_table$log2FC, "\n")

# Subsetting (filtering) rows:
# Syntax: table[ row_condition , column_selection ]
# e.g., Find genes where P_Value < 0.05 AND log2FC > 1:
sig_up <- toy_de_table[toy_de_table$P_Value < 0.05 & toy_de_table$log2FC > 1, ]
cat("\nSignificantly Upregulated Genes (P < 0.05 & log2FC > 1):\n")
print(sig_up)

# ------------------------------------------------------------------------------
# SECTION 3: LOADING REAL DATA FROM HUNT ET AL. (CELL REPORTS 2021)
# ------------------------------------------------------------------------------
# Let's locate the repo's data directory
repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
supp_dir  <- file.path(repo_root, "data", "supplementary")

cat("\nLoading real published datasets from:", supp_dir, "\n")

# 1. Loading CSV tables with Base R's read.csv()
# check.names = FALSE ensures column names like 'P value' don't get changed to 'P.value'
dexa_path <- file.path(supp_dir, "Table_S1_RNAseq_Dexamethasone.csv")
llc_path  <- file.path(supp_dir, "Table_S1_RNAseq_LLC.csv")
aging_path <- file.path(supp_dir, "Table_S1_RNAseq_Aging.csv")

if (!file.exists(dexa_path)) {
  stop("CSV files not found! Ensure data/supplementary contains the exported sheets.")
}

dexa_df  <- read.csv(dexa_path, stringsAsFactors = FALSE, check.names = FALSE)
llc_df   <- read.csv(llc_path, stringsAsFactors = FALSE, check.names = FALSE)
aging_df <- read.csv(aging_path, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Loaded Dexamethasone: %d genes across %d columns\n", nrow(dexa_df), ncol(dexa_df)))
cat(sprintf("Loaded LLC Cancer:    %d genes across %d columns\n", nrow(llc_df), ncol(llc_df)))
cat(sprintf("Loaded Aging:         %d genes across %d columns\n", nrow(aging_df), ncol(aging_df)))

# ------------------------------------------------------------------------------
# SECTION 4: INSPECTING AND PREPARING COLUMNS
# ------------------------------------------------------------------------------
# Check column names:
cat("\nDexamethasone Columns:\n")
print(colnames(dexa_df))

# Notice the critical columns:
# - Gene identifier: "Gene symbol" or "Gene name" or "ID"
# - Fold Change: "log2fC" or "log2FC"
# - Significance: "P value"
# - Minus log10 P: "-log10P"

# Let's write a small helper function to standardize column names across datasets:
standardize_de_table <- function(df, condition_name) {
  # Find fold change column (handles both 'log2fC' and 'log2FC')
  fc_col <- grep("log2f?c", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  pval_col <- grep("p[ -]?val", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  gene_col <- grep("gene|symbol|id", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  clean_df <- data.frame(
    Gene = df[[gene_col]],
    log2FC = as.numeric(df[[fc_col]]),
    P_Value = as.numeric(df[[pval_col]]),
    Condition = condition_name,
    stringsAsFactors = FALSE
  )
  
  # Remove NAs and compute -log10(P)
  clean_df <- clean_df[!is.na(clean_df$log2FC) & !is.na(clean_df$P_Value), ]
  # Avoid log10(0) by capping minimum p-value at 1e-50
  clean_df$P_Value_capped <- pmax(clean_df$P_Value, 1e-50)
  clean_df$minusLog10P <- -log10(clean_df$P_Value_capped)
  
  # Add publication significance classification:
  # Up: p < 0.05 and log2FC > 1.0 (2-fold up)
  # Down: p < 0.05 and log2FC < -1.0 (2-fold down)
  clean_df$Significance <- "NS"
  clean_df$Significance[clean_df$P_Value < 0.05 & clean_df$log2FC > 1.0] <- "Up"
  clean_df$Significance[clean_df$P_Value < 0.05 & clean_df$log2FC < -1.0] <- "Down"
  
  # Convert to an ordered factor for consistent plotting legends
  clean_df$Significance <- factor(clean_df$Significance, levels = c("Up", "Down", "NS"))
  
  return(clean_df)
}

clean_dexa  <- standardize_de_table(dexa_df, "Dexamethasone")
clean_llc   <- standardize_de_table(llc_df, "LLC Cancer")
clean_aging <- standardize_de_table(aging_df, "Aging (24m vs 6m)")

cat("\nSummary of Significance in Dexamethasone:\n")
print(table(clean_dexa$Significance))

cat("\nSummary of Significance in LLC Cancer:\n")
print(table(clean_llc$Significance))

cat("\nSummary of Significance in Aging:\n")
print(table(clean_aging$Significance))

cat("\nTutorial 01 completed successfully! Ready for Tutorial 02 (Volcano plots).\n")
