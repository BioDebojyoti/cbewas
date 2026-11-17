# ──────────────────────────────────────────────────────────────────────────────
# Supplementary Figure 3: ZFP57 core-CpG box plots for ANY exposure
# (same layout as Figure 3, but the exposure is a parameter)
#
# Usage
#   Rscript Supplementary_Figure3.R pets
#   Rscript Supplementary_Figure3.R pets infection.during.pregnancy   # several
# or, inside an R session, before sourcing:
#   exposure_to_plot <- "pets"; source("Supplementary_Figure3.R")
#
# Exposure names are the column names in `meta`.
# Output: Supplementary_Figure3_<exposure>.pdf  (one file per exposure)
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
load(resultsFile)
source("helpers_zfp57_boxplot.R")

# ── Choose exposure(s) ────────────────────────────────────────────────────────
if (!exists("exposure_to_plot")) {
  args <- commandArgs(trailingOnly = TRUE)
  exposure_to_plot <- if (length(args) >= 1) args else "infection.during.pregnancy"
}

# Gene and CpG set (change goi/cpgs here to plot something else)
goi  <- "ZFP57"
cpgs <- core_cpgs_t1

# Exposures that were analysed (binary, enough samples per cell)
available <- names(all_dmc_results)

for (exp_name in exposure_to_plot) {
  
  if (!exp_name %in% colnames(meta)) {
    stop(sprintf("'%s' is not a column of meta. Analysed exposures:\n  %s",
                 exp_name, paste(available, collapse = "\n  ")))
  }
  if (!exp_name %in% available) {
    warning(sprintf("'%s' was not part of the limma analysis (not binary or too few samples). Plotting anyway.",
                    exp_name))
  }
  
  cat("Plotting", exp_name, "\n")
  
  p <- boxplot_exposure_dmr_cpg_goi(
    cpgs_for_goi = cpgs,
    goi          = goi,
    x            = exp_name
  )
  
  if (is.null(p)) {
    message("  Skipped: '", exp_name, "' does not have exactly 2 levels after dropping NAs.")
    next
  }
  
  print(p)
  
  ggplot2::ggsave(
    filename = file.path(outDir, paste0("Supplementary_Figure3_", exp_name, ".pdf")),
    plot = p, dpi = 300, height = 12, width = 18, units = "in"
  )
}
