## =========================================================
## Figure 4 supplement: PC-covariate association heatmaps
## before and after batch correction
## =========================================================
## Requires: batch_correction.R to have been run first
## (it saves PC_covariate_association_results.RData).

rm(list = ls())

library(ComplexHeatmap)
library(circlize)
library(grid)

outDir <- paste0(
  "/Users",
  "/debda22",
  "/Projects",
  "/epi",
  "/results-new",
  "/batch-correction-new"
)

# Loads: assoc_results (list: before, slide, slide_plus_array), col_fun
load(file = paste0(outDir, "/PC_covariate_association_results.RData"))

# Plot heatmap with legend restricted to 0.05 -> 0
plot_assoc_heatmap <- function(res, title, file) {
  sig_mat   <- res$heatmap_sig_mat
  assoc_mat <- res$heatmap_assoc_mat

  pdf(file)
  par(mfrow = c(1,1))

  ht <- ComplexHeatmap::pheatmap(
    sig_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = assoc_mat,
    main = title,
    fontsize_number = 12,
    column_names_side = "top",
    row_names_side = "left",
    angle_col = "45",
    silent = FALSE,
    col = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (sig_mat[i, j] < 0.05) {
        grid.rect(x, y, width, height, gp = gpar(col = "red", lwd = 2, fill = NA))
      }
    },
    heatmap_legend_param = list(
      title = "p-value",
      at = c(0.001, 0.01, 0.05),      # show only meaningful p-value breaks
      labels = c("0.001", "0.01", "0.05")
    )
  )
  # explicit draw(): auto-printing doesn't happen inside a function
  ComplexHeatmap::draw(ht)

  dev.off()
}

# Before batch correction
plot_assoc_heatmap(
  assoc_results$before,
  title = "PC-Covariate Associations (No Batch Correction)",
  file  = paste0(outDir, "/heatmap_before_batch_correction.pdf")
)

# After batch correction (Slide + Array)
plot_assoc_heatmap(
  assoc_results$slide_plus_array,
  title = "PC-Covariate Associations (Batch Correction: Slide + Array)",
  file  = paste0(outDir, "/heatmap_slide_plus_array_batch_correction.pdf")
)

# Optional intermediate step (Slide only), uncomment if wanted:
# plot_assoc_heatmap(
#   assoc_results$slide,
#   title = "PC-Covariate Associations (Batch Correction: Slide)",
#   file  = paste0(outDir, "/heatmap_slide_batch_correction.pdf")
# )

sink(paste0(outDir, "/session.info.figure4.supplement.txt"))
print(sessionInfo())
sink()
