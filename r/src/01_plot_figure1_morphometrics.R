# ==============================================================================
# TUTORIAL 07: Figure 1 Reproduction — Morphometric & Myofiber Phenotyping
# Topic: Physiological & Histological Phenotyping across Dexamethasone, LLC, and Aging
# Reproduces: Figure 1A, Figure 1B, and Figure 1C (Hunt et al., Cell Reports 2021)
# ==============================================================================

cat("=================================================================\n")
cat("Starting Figure 1 Morphometrics & Myofiber Phenotyping Tutorial\n")
cat("=================================================================\n")

repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
supp_dir  <- file.path(repo_root, "data", "supplementary")
fig_dir   <- file.path(repo_root, "tutorial", "figures")
r_fig_dir <- file.path(repo_root, "r", "figures")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(r_fig_dir)) dir.create(r_fig_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# HELPER 1: DRAW INDIVIDUAL MORPHOMETRIC BARPLOT WITH JITTER & STATS
# ------------------------------------------------------------------------------
draw_morph_bar <- function(ctrl_vals, trt_vals, title, ylab, trt_name, trt_col, ylim_user = NULL) {
  ctrl_vals <- na.omit(as.numeric(ctrl_vals))
  trt_vals  <- na.omit(as.numeric(trt_vals))
  
  m_ctrl <- mean(ctrl_vals)
  m_trt  <- mean(trt_vals)
  
  se_ctrl <- sd(ctrl_vals) / sqrt(length(ctrl_vals))
  se_trt  <- sd(trt_vals) / sqrt(length(trt_vals))
  
  # Two-tailed Student's t-test p-value
  pval <- t.test(ctrl_vals, trt_vals)$p.value
  p_label <- if (pval < 0.001) "P<=0.001" else if (pval < 0.01) sprintf("P=%.3f", pval) else sprintf("P=%.3f", pval)
  
  # Percent change
  pct_change <- ((m_trt - m_ctrl) / m_ctrl) * 100
  pct_label <- sprintf("%+.0f%%", pct_change)
  
  # Determine y-limits
  all_max <- max(c(ctrl_vals, trt_vals, m_ctrl + se_ctrl, m_trt + se_trt), na.rm = TRUE)
  y_max <- if (!is.null(ylim_user)) ylim_user[2] else if (is.finite(all_max) && all_max > 0) all_max * 1.25 else 100
  y_min <- 0
  
  # Draw empty barplot
  bp <- barplot(c(m_ctrl, m_trt), names.arg = c("Ctrl", trt_name),
                col = "white", border = c("black", trt_col), lwd = 1.5,
                ylim = c(y_min, y_max),
                ylab = ylab, main = title,
                font.main = 2, font.lab = 2, cex.names = 0.9, cex.axis = 0.85,
                las = 1)
  
  # Error bars (SEM)
  arrows(bp[1], m_ctrl - se_ctrl, bp[1], m_ctrl + se_ctrl, length = 0.04, angle = 90, code = 3, lwd = 1.2)
  arrows(bp[2], m_trt - se_trt, bp[2], m_trt + se_trt, length = 0.04, angle = 90, code = 3, lwd = 1.2, col = trt_col)
  
  # Individual animal data points (jittered dots)
  set.seed(42)
  points(jitter(rep(bp[1], length(ctrl_vals)), amount = 0.08), ctrl_vals,
         pch = 21, bg = "grey80", col = "grey30", cex = 0.8)
  points(jitter(rep(bp[2], length(trt_vals)), amount = 0.08), trt_vals,
         pch = 21, bg = trt_col, col = trt_col, cex = 0.8)
  
  # P-value at the top
  text(x = mean(bp), y = y_max * 0.93, labels = p_label, font = 2, cex = 0.8)
  
  # Percentage change at bottom-right
  text(x = bp[2] + 0.35, y = y_max * 0.08, labels = pct_label, font = 2, cex = 0.85, adj = 1)
}

# ------------------------------------------------------------------------------
# HELPER 2: PARSE FIBER TYPE DATA (2A, 2X, 2B)
# ------------------------------------------------------------------------------
parse_fiber_type_csv <- function(csv_path) {
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  # Row 1 = 2A, Row 2 = 2X, Row 3 = 2B
  # Determine where Treatment begins
  n_cols <- ncol(df)
  # Find column where second header or non-empty string occurs
  header_names <- colnames(df)
  trt_col_idx <- which(grepl("Dexa|LLC|Old", header_names, ignore.case = TRUE))[1]
  
  if (is.na(trt_col_idx)) {
    # If not in header, find in row 1
    trt_col_idx <- which(grepl("Dexa|LLC|Old", as.character(df[1, ])))[1]
    if (is.na(trt_col_idx)) trt_col_idx <- ceiling(n_cols / 2) + 1
  }
  
  res <- list()
  for (r in 1:3) {
    fiber_name <- as.character(df[r, 1])
    if (is.na(fiber_name) || fiber_name == "") fiber_name <- c("2A", "2X", "2B")[r]
    
    ctrl_v <- na.omit(as.numeric(df[r, 2:(trt_col_idx - 1)]))
    trt_v  <- na.omit(as.numeric(df[r, trt_col_idx:n_cols]))
    
    res[[fiber_name]] <- list(ctrl = ctrl_v, trt = trt_v)
  }
  return(res)
}

# ------------------------------------------------------------------------------
# HELPER 3: DRAW GROUPED FIBER TYPE BARPLOT (2A, 2X, 2B)
# ------------------------------------------------------------------------------
draw_fiber_type_grouped <- function(fiber_data, title, ylab, trt_name, trt_col, ylim_user = NULL) {
  means_ctrl <- sapply(fiber_data, function(x) mean(x$ctrl))
  means_trt  <- sapply(fiber_data, function(x) mean(x$trt))
  
  se_ctrl <- sapply(fiber_data, function(x) sd(x$ctrl) / sqrt(length(x$ctrl)))
  se_trt  <- sapply(fiber_data, function(x) sd(x$trt) / sqrt(length(x$trt)))
  
  # Matrix for grouped barplot (2 rows: Ctrl, Trt; 3 cols: 2A, 2X, 2B)
  mat_means <- rbind(means_ctrl, means_trt)
  
  all_max <- max(mat_means + rbind(se_ctrl, se_trt), na.rm = TRUE)
  y_max <- if (!is.null(ylim_user)) ylim_user[2] else if (is.finite(all_max) && all_max > 0) all_max * 1.25 else 100
  
  bp <- barplot(mat_means, beside = TRUE,
                col = "white", border = c("black", trt_col), lwd = 1.5,
                ylim = c(0, y_max),
                ylab = ylab, main = title,
                font.main = 2, font.lab = 2, cex.names = 0.9, cex.axis = 0.85,
                las = 1)
  
  # Error bars and jittered points
  p_vals <- numeric(3)
  for (i in 1:3) {
    # Ctrl
    c_x <- bp[1, i]
    c_vals <- fiber_data[[i]]$ctrl
    arrows(c_x, means_ctrl[i] - se_ctrl[i], c_x, means_ctrl[i] + se_ctrl[i], length = 0.03, angle = 90, code = 3)
    points(jitter(rep(c_x, length(c_vals)), amount = 0.08), c_vals, pch = 21, bg = "grey80", col = "grey30", cex = 0.7)
    
    # Trt
    t_x <- bp[2, i]
    t_vals <- fiber_data[[i]]$trt
    arrows(t_x, means_trt[i] - se_trt[i], t_x, means_trt[i] + se_trt[i], length = 0.03, angle = 90, code = 3, col = trt_col)
    points(jitter(rep(t_x, length(t_vals)), amount = 0.08), t_vals, pch = 21, bg = trt_col, col = trt_col, cex = 0.7)
    
    # P-value
    pval <- t.test(c_vals, t_vals)$p.value
    p_vals[i] <- pval
  }
  
  # Draw bottom p-values matching Hunt et al. layout
  p_text <- sprintf("P=  %.3f  %.3f  %.3f", p_vals[1], p_vals[2], p_vals[3])
  mtext(p_text, side = 1, line = 2.2, cex = 0.65, font = 2)
}

# ------------------------------------------------------------------------------
# STEP 4: PLOT FIGURE 1A, 1B, 1C
# ------------------------------------------------------------------------------
plot_figure_1_condition <- function(cond_id, cond_name, trt_label, trt_color, out_filename) {
  cat(sprintf("Generating %s...\n", cond_name))
  
  # Explicit file mapping
  files_map <- list(
    "Dexa" = list(
      body = "Fig1A_dexabody.csv",
      tibia = "Fig1A_dexatibia.csv",
      ta = "Fig1A_dexaTA.csv",
      ta_body = "Fig1A_dexaTA_body.csv",
      ta_tibia = "Fig1A_dexaTA_tibia.csv",
      all_fnum = "Fig1A_dexa_all_fiber_number.csv",
      all_fsize = "Fig1A_dexa_all_fiber_size.csv",
      fnum = "Fig1A_dexa_fiber_number.csv",
      ftype = "Fig1A_dexa_fiber_type.csv",
      fsize = "Fig1A_dexa_fiber_size.csv"
    ),
    "LLC" = list(
      body = "Fig1B_LLC_body.csv",
      tibia = "Fig1B_LLC_tibia.csv",
      ta = "Fig1B_LLC_TA.csv",
      ta_body = "Fig1B_LLC_TA_body.csv",
      ta_tibia = "Fig1B_LLC_TA_tibia.csv",
      all_fnum = "Fig1B_LLC_total_fiber.csv",
      all_fsize = "Fig1B_LLC_total_size.csv",
      fnum = "Fig1B_LLC_fiber_number.csv",
      ftype = "Fig1B_LLC_fiber_type.csv",
      fsize = "Fig1B_LLC_fiber_size.csv"
    ),
    "Aging" = list(
      body = "Fig1C_old_body.csv",
      tibia = "Fig1C_old_tibia.csv",
      ta = "Fig1C_old_TA.csv",
      ta_body = "Fig1C_old_TA_body.csv",
      ta_tibia = "Fig1C_old_TA_tibia.csv",
      all_fnum = "Fig1C_old_all_fiber_number.csv",
      all_fsize = "Fig1C_old_all_fiber_size.csv",
      fnum = "Fig1C_old_fiber_number.csv",
      ftype = "Fig1C_old_fiber_type.csv",
      fsize = "Fig1C_old_fiber_size.csv"
    )
  )[[cond_id]]
  
  png(out_filename, width = 3600, height = 1800, res = 300)
  par(mfrow = c(2, 5), mar = c(3.5, 4.2, 3, 1), bg = "white")
  
  # Row 1: 5 Gross Morphometric Parameters
  # 1. Body mass
  df_body <- read.csv(file.path(supp_dir, files_map$body), check.names = FALSE)
  draw_morph_bar(df_body[[1]], df_body[[2]], "Body mass", "body mass (g)", trt_label, trt_color)
  
  # 2. Tibia length
  df_tibia <- read.csv(file.path(supp_dir, files_map$tibia), check.names = FALSE)
  draw_morph_bar(df_tibia[[1]], df_tibia[[2]], "Tibia length", "Tibia length (mm)", trt_label, trt_color)
  
  # 3. TA mass
  df_ta <- read.csv(file.path(supp_dir, files_map$ta), check.names = FALSE)
  draw_morph_bar(df_ta[[1]], df_ta[[2]], "TA mass", "TA mass (mg)", trt_label, trt_color)
  
  # 4. TA / body mass
  df_tabody <- read.csv(file.path(supp_dir, files_map$ta_body), check.names = FALSE)
  draw_morph_bar(df_tabody[[1]], df_tabody[[2]], "TA/body mass", "TA mass/body mass (mg/g)", trt_label, trt_color)
  
  # 5. TA / tibia length
  df_tatibia <- read.csv(file.path(supp_dir, files_map$ta_tibia), check.names = FALSE)
  draw_morph_bar(df_tatibia[[1]], df_tatibia[[2]], "TA/tibia length", "TA mass/tibia length\n(mg/mm)", trt_label, trt_color)
  
  # Row 2: 5 Histological / Myofiber Parameters
  # 6. All myofiber number
  df_fnum <- read.csv(file.path(supp_dir, files_map$all_fnum), check.names = FALSE)
  draw_morph_bar(df_fnum[[1]], df_fnum[[2]], "All myofiber number", "total myofiber number", trt_label, trt_color)
  
  # 7. All myofiber size
  df_fsize <- read.csv(file.path(supp_dir, files_map$all_fsize), check.names = FALSE)
  draw_morph_bar(df_fsize[[1]], df_fsize[[2]], "All myofiber size", "Feret's minimal diameter (µm)", trt_label, trt_color)
  
  # 8. Myofiber number by type
  fiber_num_data <- parse_fiber_type_csv(file.path(supp_dir, files_map$fnum))
  draw_fiber_type_grouped(fiber_num_data, "Myofiber number", "number of myofibers", trt_label, trt_color)
  
  # 9. Myofiber proportion by type
  fiber_prop_data <- parse_fiber_type_csv(file.path(supp_dir, files_map$ftype))
  draw_fiber_type_grouped(fiber_prop_data, "Myofiber proportion", "% of myofibers", trt_label, trt_color, ylim_user = c(0, 100))
  
  # 10. Myofiber size by type
  fiber_size_data <- parse_fiber_type_csv(file.path(supp_dir, files_map$fsize))
  draw_fiber_type_grouped(fiber_size_data, "Myofiber size", "Feret's minimal diameter (µm)", trt_label, trt_color)
  
  dev.off()
  cat(sprintf("Saved: %s\n", out_filename))
}

# Generate Figure 1A: Dexamethasone
plot_figure_1_condition(
  cond_id = "Dexa",
  cond_name = "Figure 1A: Dexamethasone Atrophy Phenotype",
  trt_label = "Dexa",
  trt_color = "#e41a1c",
  out_filename = file.path(fig_dir, "Tutorial_Figure1A_Dexamethasone_Morphometrics.png")
)
file.copy(file.path(fig_dir, "Tutorial_Figure1A_Dexamethasone_Morphometrics.png"),
          file.path(r_fig_dir, "Figure1A_Dexamethasone_Morphometrics.png"), overwrite = TRUE)

# Generate Figure 1B: LLC Cancer Cachexia
plot_figure_1_condition(
  cond_id = "LLC",
  cond_name = "Figure 1B: LLC Cancer Cachexia Phenotype",
  trt_label = "LLC",
  trt_color = "#377eb8",
  out_filename = file.path(fig_dir, "Tutorial_Figure1B_LLC_Morphometrics.png")
)
file.copy(file.path(fig_dir, "Tutorial_Figure1B_LLC_Morphometrics.png"),
          file.path(r_fig_dir, "Figure1B_LLC_Morphometrics.png"), overwrite = TRUE)

# Generate Figure 1C: Aging (Sarcopenia)
plot_figure_1_condition(
  cond_id = "Aging",
  cond_name = "Figure 1C: Sarcopenic Aging Phenotype",
  trt_label = "Old",
  trt_color = "#4daf4a",
  out_filename = file.path(fig_dir, "Tutorial_Figure1C_Aging_Morphometrics.png")
)
file.copy(file.path(fig_dir, "Tutorial_Figure1C_Aging_Morphometrics.png"),
          file.path(r_fig_dir, "Figure1C_Aging_Morphometrics.png"), overwrite = TRUE)

cat("\nTutorial 07 completed successfully! Figure 1A, 1B, and 1C generated.\n")
cat("=================================================================\n")
