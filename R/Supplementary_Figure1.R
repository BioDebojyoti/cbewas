# ──────────────────────────────────────────────────────────────────────────────
# Supplementary Figure (nominal only): Interaction EWAS signal across all
# prenatal exposures — CpGs with nominal p < 0.05 only (no |Δβ| filter).
# Same as Figure1.R without the |Δβ| > 0.10 stage.
# Output: Supplementary_Figure1.pdf
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
library(showtext)
load(resultsFile)

sfig_df <- ewas_log %>%
  dplyr::mutate(
    exposure_label = gsub("\\.", " ", gsub("\\.during\\.pregnancy", "", exposure)),
    min_cell = pmin(n_T1D_exp, n_T1D_unexp, n_ctrl_exp, n_ctrl_unexp),
    sample_distribution = case_when(
      min_cell < SAMPLE_SPARSE_THR   ~ "Sparse (bias risk)",
      min_cell < SAMPLE_MODERATE_THR ~ "Constrained (variance risk)",
      TRUE                           ~ "Adequate"
    ),
    sample_distribution = factor(
      sample_distribution,
      levels = c("Adequate", "Constrained (variance risk)", "Sparse (bias risk)")
    ),
    n_cpgs = n_nominal_p05
  )

sfig <- ggplot(
  sfig_df,
  aes(x = reorder(exposure_label, n_cpgs), y = n_cpgs, alpha = sample_distribution)
) +
  geom_col(fill = "#2C7BB6", width = 0.65) +
  scale_alpha_manual(
    values = c(
      "Adequate"                    = 1.0,
      "Constrained (variance risk)" = 0.75,
      "Sparse (bias risk)"          = 0.5
    ),
    name = "Sample distribution adequacy: "
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title         = element_text(size = 18, face = "bold"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title         = element_text(size = 12, face = "bold"),
    strip.text         = element_text(size = 12, face = "bold"),
    axis.text          = element_text(size = 12, face = "bold"),
    plot.caption       = element_text(size = 12)
  ) +
  labs(
    x        = NULL,
    y        = "Number of CpGs (p < 0.05)",
    title    = "Interaction EWAS signal across 22 prenatal exposures",
    subtitle = "CpGs with nominal p < 0.05  |  Transparency encodes sample distribution adequacy"
  )

sfig

# showtext handles fonts containing Greek glyphs
showtext_auto()
font_add("Arial", regular = "/Library/Fonts/Arial.ttf")

ggplot2::ggsave(
  filename = file.path(outDir, "Supplementary_Figure1.pdf"),
  plot = sfig, dpi = 300, height = 12, width = 18, units = "in"
)

showtext_auto(FALSE)
