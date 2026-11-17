# ──────────────────────────────────────────────────────────────────────────────
# Figure 4: Box plots of ZFP57 core CpGs (infection during pregnancy)
# Output: Figure4.pdf
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
load(resultsFile)
source("helpers_zfp57_boxplot.R")

fig3 <- boxplot_exposure_dmr_cpg_goi(
  cpgs_for_goi = core_cpgs_t1,
  goi          = "ZFP57",
  x            = "infection.during.pregnancy"
)

fig3

ggplot2::ggsave(
  filename = file.path(outDir, "Figure4.pdf"),
  plot = fig3, dpi = 300, height = 12, width = 18, units = "in"
)
