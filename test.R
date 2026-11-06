# ==============================================================================
# Quick Start / Interactive Test Script
# Run the complete R reproduction pipeline for Hunt et al. (Cell Reports 2021)
# ==============================================================================

# Ensure working directory is repository root
if (!dir.exists("r/src") && dir.exists("../r/src")) {
  setwd("..")
}

cat("Running R reproduction pipeline from root...\n")
source("r/src/run_all.R")
