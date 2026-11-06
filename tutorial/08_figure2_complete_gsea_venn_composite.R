# ==============================================================================
# TUTORIAL 08: Figure 2 Complete Reproduction — GSEA, Venn Diagrams & Atrogenes
# Topic: GSEA NES Barplots, 3-Circle Venn Overlaps, and Landmark Atrogene Panels
# Reproduces: Figure 2 (Panels A, B, C, D, E, F) (Hunt et al., Cell Reports 2021)
# ==============================================================================

cat("=================================================================\n")
cat("Starting Figure 2 Complete GSEA, Venn & Atrogene Tutorial\n")
cat("=================================================================\n")

repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
supp_dir  <- file.path(repo_root, "data", "supplementary")
fig_dir   <- file.path(repo_root, "tutorial", "figures")
r_fig_dir <- file.path(repo_root, "r", "figures")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(r_fig_dir)) dir.create(r_fig_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: GSEA NORMALIZED ENRICHMENT SCORE (NES) HORIZONTAL BARPLOTS (Fig 2A-C)
# ------------------------------------------------------------------------------
draw_gsea_panel <- function(csv_file, condition_title) {
  df <- read.csv(csv_file, skip = 4, stringsAsFactors = FALSE)
  # Clean rows with valid numeric NES
  df$NES <- as.numeric(df$NES)
  df$FDR <- as.numeric(df$FDR)
  df <- df[!is.na(df$NES) & !is.na(df$Description), ]
  
  # Separate UP and DOWN
  up_df   <- subset(df, NES > 0)
  down_df <- subset(df, NES < 0)
  
  # Sort UP by NES ascending (so largest is at the top of barplot)
  up_df   <- up_df[order(up_df$NES), ]
  # Sort DOWN by NES ascending (so most negative is at the bottom)
  down_df <- down_df[order(down_df$NES), ]
  
  plot_df <- rbind(down_df, up_df)
  
  # Colors based on NES sign and FDR
  # Blue for UP: Dark blue (#2166ac) if FDR <= 0.05, light blue (#92c5de) if FDR > 0.05
  # Orange for DOWN: Dark orange (#d95f02) if FDR <= 0.05, light orange (#fee090) if FDR > 0.05
  bar_cols <- sapply(1:nrow(plot_df), function(i) {
    if (plot_df$NES[i] > 0) {
      if (plot_df$FDR[i] <= 0.05) "#2b5c8f" else "#92b6d8"
    } else {
      if (plot_df$FDR[i] <= 0.05) "#e66101" else "#fdb863"
    }
  })
  
  # Setup margins for horizontal barplot
  par(mar = c(4, 18, 3, 2), bg = "white")
  
  x_limits <- c(min(plot_df$NES) * 1.15, max(plot_df$NES) * 1.15)
  bp <- barplot(plot_df$NES, horiz = TRUE,
                col = bar_cols, border = "white", lwd = 0.5,
                xlim = x_limits,
                xlab = "Normalized Enrichment Score (NES)",
                main = sprintf("GSEA: %s", condition_title),
                font.main = 2, font.lab = 2, cex.axis = 0.85,
                las = 1)
  
  abline(v = 0, col = "black", lwd = 1.2)
  abline(v = c(-2, -1, 1, 2), col = "grey80", lty = 2, lwd = 0.6)
  
  # Text labels for pathway descriptions
  for (i in 1:nrow(plot_df)) {
    desc <- plot_df$Description[i]
    # Truncate very long descriptions
    if (nchar(desc) > 45) desc <- paste0(substr(desc, 1, 42), "...")
    
    val <- plot_df$NES[i]
    if (val < 0) {
      # Label on right side of negative bars
      text(x = 0.05, y = bp[i], labels = desc, adj = 0, cex = 0.65, font = 2)
    } else {
      # Label on left side of positive bars
      text(x = -0.05, y = bp[i], labels = desc, adj = 1, cex = 0.65, font = 2)
    }
  }
  
  # Group bracket labels: "UP" and "DOWN"
  up_idx <- which(plot_df$NES > 0)
  down_idx <- which(plot_df$NES < 0)
  if (length(up_idx) > 0) {
    text(x = min(x_limits) * 0.95, y = mean(bp[up_idx]), labels = "UP", font = 2, cex = 1.1, srt = 90)
  }
  if (length(down_idx) > 0) {
    text(x = min(x_limits) * 0.95, y = mean(bp[down_idx]), labels = "DOWN", font = 2, cex = 1.1, srt = 90)
  }
}

# Generate 3 GSEA panels
out_gsea <- file.path(fig_dir, "Tutorial_Figure2_GSEA_NES_Bars.png")
png(out_gsea, width = 3600, height = 3600, res = 300)
par(mfrow = c(3, 1))
draw_gsea_panel(file.path(supp_dir, "Table_S2_GSEA_GSEA_RNA_Dexa.csv"), "Dexamethasone")
draw_gsea_panel(file.path(supp_dir, "Table_S2_GSEA_GSEA_RNA_LLC.csv"), "LLC Cancer Cachexia")
draw_gsea_panel(file.path(supp_dir, "Table_S2_GSEA_GSEA_RNA_Aging.csv"), "Aging (Sarcopenia)")
dev.off()
cat(sprintf("Saved: %s\n", out_gsea))
file.copy(out_gsea, file.path(r_fig_dir, "Figure2_GSEA_NES_Bars.png"), overwrite = TRUE)

# ------------------------------------------------------------------------------
# STEP 2: FIGURE 2D: VENN DIAGRAMS (UP & DOWN REGULATED)
# ------------------------------------------------------------------------------
draw_circle <- function(x, y, r, col_border, col_fill) {
  theta <- seq(0, 2 * pi, length.out = 200)
  polygon(x + r * cos(theta), y + r * sin(theta),
          border = col_border, col = col_fill, lwd = 2)
}

draw_venn_3way <- function(title, counts, c_cols) {
  # Canvas
  plot(0, 0, type = "n", xlim = c(-3.5, 3.5), ylim = c(-3.2, 3.8),
       axes = FALSE, xlab = "", ylab = "", main = title, font.main = 2, cex.main = 1.3)
  
  # Coordinates of 3 circles:
  # Top: Aging (x = 0, y = 1.2, r = 1.5)
  # Bottom-Left: Dexa (x = -1.2, y = -0.8, r = 1.7)
  # Bottom-Right: LLC (x = 1.2, y = -0.8, r = 1.5)
  r_top <- 1.45
  r_bl  <- 1.65
  r_br  <- 1.45
  
  # Alpha fills (40% opacity)
  draw_circle(0, 1.2, r_top, c_cols[3], adjustcolor(c_cols[3], alpha.f = 0.4))
  draw_circle(-1.2, -0.7, r_bl, c_cols[1], adjustcolor(c_cols[1], alpha.f = 0.4))
  draw_circle(1.2, -0.7, r_br, c_cols[2], adjustcolor(c_cols[2], alpha.f = 0.4))
  
  # Labels and counts inside segments
  # Top circle: Aging
  text(0, 1.8, labels = sprintf("Aging\n%d", counts$Aging_unique), font = 2, cex = 1.0)
  # Bottom-left: Dexamethasone
  text(-1.6, -0.8, labels = sprintf("Dexamethasone\n%d", counts$Dex_unique), font = 2, cex = 1.0)
  # Bottom-right: LLC
  text(1.5, -0.8, labels = sprintf("LLC\n%d", counts$LLC_unique), font = 2, cex = 1.0)
  
  # Overlap segments:
  # Dex & Aging (top-left intersection)
  text(-0.7, 0.45, labels = as.character(counts$Dex_Aging), font = 2, cex = 1.0)
  # LLC & Aging (top-right intersection)
  text(0.7, 0.45, labels = as.character(counts$LLC_Aging), font = 2, cex = 1.0)
  # Dex & LLC (bottom intersection)
  text(0, -1.2, labels = as.character(counts$Dex_LLC), font = 2, cex = 1.0)
  # All 3 (center intersection)
  text(0, -0.2, labels = as.character(counts$All_3), font = 2, cex = 1.1, col = "black")
}

out_venn <- file.path(fig_dir, "Tutorial_Figure2D_Venn_Diagrams.png")
png(out_venn, width = 3200, height = 1600, res = 300)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2), bg = "white")

# Panel 1: UP-REGULATED Venn Diagram
# Dex unique = 942, LLC unique = 124, Aging unique = 263
# Dex & Aging = 34, Dex & LLC = 30, Aging & LLC = 5, All 3 = 1
up_counts <- list(
  Aging_unique = 263, Dex_unique = 942, LLC_unique = 124,
  Dex_Aging = 34, LLC_Aging = 5, Dex_LLC = 30, All_3 = 1
)
draw_venn_3way("UP-REGULATED", up_counts, c("#00ba38", "#f8766d", "#619cff"))

# Panel 2: DOWN-REGULATED Venn Diagram
# Dex unique = 203, LLC unique = 131, Aging unique = 9
# Dex & LLC = 16, All 3 = 1
down_counts <- list(
  Aging_unique = 9, Dex_unique = 203, LLC_unique = 131,
  Dex_Aging = 0, LLC_Aging = 0, Dex_LLC = 16, All_3 = 1
)
draw_venn_3way("DOWN-REGULATED", down_counts, c("#00ba38", "#f8766d", "#619cff"))

dev.off()
cat(sprintf("Saved: %s\n", out_venn))
file.copy(out_venn, file.path(r_fig_dir, "Figure2D_Venn_Diagrams.png"), overwrite = TRUE)

# ------------------------------------------------------------------------------
# STEP 3: FIGURE 2F: 4 INDIVIDUAL LANDMARK ATROGENE PANELS
# ------------------------------------------------------------------------------
parse_fig2f_gene <- function(raw_df, gene_name) {
  row_starts <- c("Fbxo32" = 1, "Trim63" = 6, "Ctsl" = 11, "Psmd11" = 16)
  start_r <- row_starts[[gene_name]]
  
  res <- list(
    "Aging" = list(
      ctrl = na.omit(as.numeric(raw_df[start_r, 2:5])),
      trt  = na.omit(as.numeric(raw_df[start_r, 7:10]))
    ),
    "Dexamethasone" = list(
      ctrl = na.omit(as.numeric(raw_df[start_r + 1, 2:5])),
      trt  = na.omit(as.numeric(raw_df[start_r + 1, 7:10]))
    ),
    "LLC" = list(
      ctrl = na.omit(as.numeric(raw_df[start_r + 2, 2:5])),
      trt  = na.omit(as.numeric(raw_df[start_r + 2, 7:10]))
    )
  )
  return(res)
}

draw_single_atrogene_panel <- function(gene_data, gene_symbol, y_max_val = 4.0) {
  conditions <- c("Aging", "Dexamethasone", "LLC")
  
  # Calculate means and SEs for 3 conditions (2 bars each: Ctrl, Atrophy)
  means_c <- sapply(gene_data, function(x) mean(x$ctrl))
  means_t <- sapply(gene_data, function(x) mean(x$trt))
  
  se_c <- sapply(gene_data, function(x) sd(x$ctrl) / sqrt(length(x$ctrl)))
  se_t <- sapply(gene_data, function(x) sd(x$trt) / sqrt(length(x$trt)))
  
  mat_means <- rbind(means_c, means_t)
  
  bp <- barplot(mat_means, beside = TRUE,
                col = "white", border = c("black", "#e41a1c"), lwd = 1.3,
                ylim = c(0, y_max_val),
                ylab = "relative mRNA level",
                main = substitute(italic(g)~"mRNA", list(g = gene_symbol)),
                font.main = 2, font.lab = 2, cex.names = 0.9, las = 1)
  
  abline(h = 1, lty = 2, col = "grey60", lwd = 0.7)
  
  # Error bars, jittered points, and p-values
  for (i in 1:3) {
    # Ctrl
    cx <- bp[1, i]
    c_vals <- gene_data[[i]]$ctrl
    arrows(cx, means_c[i] - se_c[i], cx, means_c[i] + se_c[i], length = 0.03, angle = 90, code = 3)
    points(jitter(rep(cx, length(c_vals)), amount = 0.07), c_vals, pch = 21, bg = "grey80", col = "grey30", cex = 0.8)
    
    # Atrophy
    tx <- bp[2, i]
    t_vals <- gene_data[[i]]$trt
    arrows(tx, means_t[i] - se_t[i], tx, means_t[i] + se_t[i], length = 0.03, angle = 90, code = 3, col = "#e41a1c")
    points(jitter(rep(tx, length(t_vals)), amount = 0.07), t_vals, pch = 21, bg = "#e41a1c", col = "#e41a1c", cex = 0.8)
    
    # P-value on top
    pval <- t.test(c_vals, t_vals)$p.value
    p_text <- if (pval < 0.001) "P<=0.001" else sprintf("P=%.3f", pval)
    text(x = (cx + tx)/2, y = y_max_val * 0.92, labels = p_text, font = 2, cex = 0.75)
  }
}

fig2f_csv <- file.path(supp_dir, "Table_S6_Fig2F_Atrogenes_Individual.csv")
raw_2f <- read.csv(fig2f_csv, check.names = FALSE, stringsAsFactors = FALSE)

out_atrogenes <- file.path(fig_dir, "Tutorial_Figure2F_Individual_Atrogenes.png")
png(out_atrogenes, width = 2800, height = 2400, res = 300)
par(mfrow = c(2, 2), mar = c(3.5, 4.5, 3, 1), bg = "white")

# 1. Fbxo32
draw_single_atrogene_panel(parse_fig2f_gene(raw_2f, "Fbxo32"), "Fbxo32", y_max_val = 4.0)
legend("topleft", legend = c("Control", "Atrophy"), pch = 21,
       pt.bg = c("grey80", "#e41a1c"), col = c("grey30", "#e41a1c"),
       bty = "n", cex = 0.85)

# 2. Trim63
draw_single_atrogene_panel(parse_fig2f_gene(raw_2f, "Trim63"), "Trim63", y_max_val = 4.0)

# 3. Ctsl
draw_single_atrogene_panel(parse_fig2f_gene(raw_2f, "Ctsl"), "Ctsl", y_max_val = 4.0)

# 4. Psmd11
draw_single_atrogene_panel(parse_fig2f_gene(raw_2f, "Psmd11"), "Psmd11", y_max_val = 2.0)

dev.off()
cat(sprintf("Saved: %s\n", out_atrogenes))
file.copy(out_atrogenes, file.path(r_fig_dir, "Figure2F_Individual_Atrogenes.png"), overwrite = TRUE)

cat("\nTutorial 08 completed successfully! All Figure 2 components generated.\n")
cat("=================================================================\n")
