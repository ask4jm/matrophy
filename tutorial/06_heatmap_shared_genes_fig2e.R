# ==============================================================================
# TUTORIAL 06: Heatmap of Shared DEGs Across Atrophy Modes
# Topic: Genomic Heatmaps, Red-White-Blue Color Palettes, NA Shading & Annotations
# Reproduces: Figure 2E (Hunt et al., Cell Reports 2021)
# ==============================================================================

repo_root <- if (file.exists("data")) getwd() else if (file.exists("../data")) normalizePath("..") else normalizePath("../..")
supp_dir  <- file.path(repo_root, "data", "supplementary")
fig_dir   <- file.path(repo_root, "tutorial", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

cat("=================================================================\n")
cat("Starting Figure 2E Heatmap Tutorial\n")
cat("=================================================================\n")

# ------------------------------------------------------------------------------
# STEP 1: LOAD SOURCE DATA
# ------------------------------------------------------------------------------
fig2e_csv <- file.path(supp_dir, "Table_S6_Fig2E_Heatmap.csv")
if (!file.exists(fig2e_csv)) {
  stop("Missing Table_S6_Fig2E_Heatmap.csv! Ensure it is exported from Supplement 7.")
}

raw_df <- read.csv(fig2e_csv, check.names = FALSE, stringsAsFactors = FALSE)
cat(sprintf("Loaded Figure 2E source data: %d genes across %d columns\n", nrow(raw_df), ncol(raw_df)))

gene_symbols <- raw_df[[1]]
# Extract numeric matrix of 21 samples:
# - Aging: 3 Control, 3 Old
# - Dexa: 4 Control, 4 Dexa
# - LLC: 3 Control, 4 LLC
data_mat <- as.matrix(raw_df[, 2:22])

# Convert fold-change to log2 space:
# 0.25 -> -2 (Red)
# 1.0  ->  0 (White)
# 4.0  -> +2 (Blue)
log2_mat <- log2(data_mat)

# Cap fold change at -2 and +2 to match the paper's color scale (0.25 to 4.0 fold change)
log2_mat_capped <- pmax(pmin(log2_mat, 2.0), -2.0)

# Define the 6 functional categories and row boundaries (86 genes total)
# 1. ALL UP: rows 1:9 (9 genes)
# 2. OLD AND DEXA UP: rows 10:40 (31 genes)
# 3. OLD AND LLC UP: rows 41:45 (5 genes)
# 4. DEXA AND LLC UP: rows 46:69 (24 genes)
# 5. DEXA AND LLC DOWN: rows 70:85 (16 genes)
# 6. ALL DOWN: row 86 (1 gene: Ptx4)
category_splits <- list(
  "ALL UP" = 1:9,
  "OLD AND DEXA UP" = 10:40,
  "OLD AND LLC UP" = 41:45,
  "DEXA AND LLC UP" = 46:69,
  "DEXA AND LLC DOWN" = 70:85,
  "ALL DOWN" = 86:86
)

# ------------------------------------------------------------------------------
# STEP 2: DEFINE PUBLICATION COLOR PALETTE
# ------------------------------------------------------------------------------
# In Hunt et al. Fig 2E:
# Red (#d73027) = Downregulated (FC = 0.25)
# White (#ffffff) = Unchanged (FC = 1.0)
# Blue (#2166ac) = Upregulated (FC = 4.0)
# Dark Grey (#666666) = Not detected / missing (NA)
n_colors <- 101
palette_fun <- colorRampPalette(c("#b2182b", "#d6604d", "#f4a582", "#ffffff", "#92c5de", "#4393c3", "#2166ac"))
hm_colors <- palette_fun(n_colors)

get_color <- function(val) {
  if (is.na(val)) return("#707070") # Dark grey for missing/unexpressed in cohort
  # Map [-2, +2] to 1..101
  idx <- round((val - (-2)) / (2 - (-2)) * (n_colors - 1)) + 1
  idx <- max(1, min(n_colors, idx))
  return(hm_colors[idx])
}

# ------------------------------------------------------------------------------
# STEP 3: RENDER HIGH-RESOLUTION FIGURE 2E HEATMAP
# ------------------------------------------------------------------------------
out_png <- file.path(fig_dir, "Tutorial_Figure2E_Heatmap_Shared_DEGs.png")
cat(sprintf("Rendering Figure 2E Heatmap: %s...\n", out_png))

png(out_png, width = 2400, height = 4800, res = 300)

# Setup canvas margins: left for gene names, right for category labels
par(mar = c(3, 8.5, 5, 10), bg = "white")

n_rows <- nrow(log2_mat_capped)
n_cols <- ncol(log2_mat_capped)

# Initialize blank plot
plot(0, 0, type = "n",
     xlim = c(0, n_cols), ylim = c(0, n_rows),
     xaxs = "i", yaxs = "i", axes = FALSE,
     xlab = "", ylab = "", main = "")

# Draw cells from top to bottom (row 1 at top)
for (r in 1:n_rows) {
  y_bottom <- n_rows - r
  y_top <- y_bottom + 1
  
  for (c in 1:n_cols) {
    x_left <- c - 1
    x_right <- c
    
    val <- log2_mat_capped[r, c]
    cell_col <- get_color(val)
    
    rect(x_left, y_bottom, x_right, y_top, col = cell_col, border = NA)
  }
  
  # Gene symbol label on the left margin
  text(x = -0.3, y = y_bottom + 0.5,
       labels = gene_symbols[r],
       adj = 1, cex = 0.55, font = 3, col = "black", xpd = TRUE)
}

# ------------------------------------------------------------------------------
# STEP 4: DRAW GRIDLINES, COLUMN LABELS, AND GROUP DIVIDERS
# ------------------------------------------------------------------------------
# Sample group boundaries:
# Aging: cols 1:3 (Ctrl), 4:6 (Old)
# Dexa:  cols 7:10 (Ctrl), 11:14 (Dexa)
# LLC:   cols 15:17 (Ctrl), 18:21 (LLC)

# Thin white column dividers
for (c in 0:n_cols) {
  abline(v = c, col = "white", lwd = 0.4)
}

# Group separator lines (solid black)
group_boundaries <- c(0, 3, 6, 10, 14, 17, 21)
for (b in group_boundaries) {
  abline(v = b, col = "black", lwd = 1.2)
}

# Column Headers at the top
axis_y <- n_rows + 0.5
text(x = 1.5, y = axis_y, labels = "Control", font = 2, cex = 0.8, xpd = TRUE)
text(x = 4.5, y = axis_y, labels = "Old", font = 2, cex = 0.8, xpd = TRUE)

text(x = 8.0, y = axis_y, labels = "Control", font = 2, cex = 0.8, xpd = TRUE)
text(x = 12.0, y = axis_y, labels = "Dexa", font = 2, cex = 0.8, xpd = TRUE)

text(x = 15.5, y = axis_y, labels = "Control", font = 2, cex = 0.8, xpd = TRUE)
text(x = 19.0, y = axis_y, labels = "LLC", font = 2, cex = 0.8, xpd = TRUE)

# Horizontal Category Dividers & Category Labels on the Right
for (cat_name in names(category_splits)) {
  rows_in_cat <- category_splits[[cat_name]]
  top_row <- min(rows_in_cat)
  bot_row <- max(rows_in_cat)
  
  # Top divider line
  y_line_top <- n_rows - (top_row - 1)
  abline(h = y_line_top, col = "black", lwd = 1.2)
  
  # Bottom divider line
  y_line_bot <- n_rows - bot_row
  abline(h = y_line_bot, col = "black", lwd = 1.2)
  
  # Category bracket/label text on the right
  y_mid <- (y_line_top + y_line_bot) / 2
  text(x = n_cols + 0.8, y = y_mid,
       labels = cat_name,
       adj = 0, font = 2, cex = 0.75, col = "black", xpd = TRUE)
}

# Outer border
box(lwd = 1.5)

# Panel label 'E' at top left
text(x = -6.5, y = n_rows + 2.2, labels = "E", font = 2, cex = 1.8, xpd = TRUE)

# ------------------------------------------------------------------------------
# STEP 5: COLOR BAR LEGEND (Top Right)
# ------------------------------------------------------------------------------
# Add colorbar at top right (x = 17 to 22, y = n_rows + 2 to n_rows + 3)
cb_x_start <- 16.5
cb_width <- 5.5
cb_y <- n_rows + 2.2
cb_h <- 0.8

text(x = cb_x_start + cb_width / 2, y = cb_y + 1.2,
     labels = "Fold Change", font = 2, cex = 0.65, xpd = TRUE)

cb_n <- length(hm_colors)
for (i in 1:cb_n) {
  x0 <- cb_x_start + (i - 1) * (cb_width / cb_n)
  x1 <- cb_x_start + i * (cb_width / cb_n)
  rect(x0, cb_y, x1, cb_y + cb_h, col = hm_colors[i], border = NA, xpd = TRUE)
}
rect(cb_x_start, cb_y, cb_x_start + cb_width, cb_y + cb_h, border = "black", lwd = 0.6, xpd = TRUE)

# Legend tick labels
text(x = cb_x_start, y = cb_y - 0.6, labels = "0.25", cex = 0.55, font = 2, xpd = TRUE)
text(x = cb_x_start + cb_width / 2, y = cb_y - 0.6, labels = "1", cex = 0.55, font = 2, xpd = TRUE)
text(x = cb_x_start + cb_width, y = cb_y - 0.6, labels = "4", cex = 0.55, font = 2, xpd = TRUE)

dev.off()
cat(sprintf("Figure 2E Heatmap saved successfully: %s\n", out_png))
cat("=================================================================\n")
