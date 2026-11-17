# ──────────────────────────────────────────────────────────────────────────────
# Figure 1: Interaction EWAS signal across all prenatal exposures
# Output: Figure1.pdf
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
library(showtext)
load(resultsFile)

fig1_df <- ewas_log %>%
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
    n_signal = n_nominal_p05 - 10000
  ) %>%
  tidyr::pivot_longer(
    cols      = c(n_nominal_p05, n_nom_interaction_deltabeta),
    names_to  = "filter_stage",
    values_to = "n_cpgs"
  ) %>%
  dplyr::mutate(
    filter_stage = factor(
      filter_stage,
      levels = c("n_nominal_p05", "n_nom_interaction_deltabeta"),
      labels = c("p < 0.05", "p < 0.05 & |\u0394\u03B2| > 0.10")
    )
  )

fig1 <- ggplot(
  fig1_df,
  aes(x = reorder(exposure_label, n_cpgs), y = n_cpgs,
      fill = filter_stage, alpha = sample_distribution)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(
    values = c(
      "p < 0.05"                         = "#2C7BB6",
      "p < 0.05 & |\u0394\u03B2| > 0.10" = "#D7191C"
    ),
    name = "CpG filter"
  ) +
  scale_alpha_manual(
    values = c(
      "Adequate"                    = 1.0,
      "Constrained (variance risk)" = 0.75,
      "Sparse (bias risk)"          = 0.5
    ),
    name = "Sample distribution adequacy: "
  ) +
  facet_grid(~ filter_stage, scale = "free") +
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
    y        = "Number of CpGs",
    title    = "Interaction EWAS signal across 22 prenatal exposures",
    subtitle = "Transparency encodes sample distribution adequacy"
  )

fig1

# showtext handles fonts containing Greek glyphs
showtext_auto()
font_add("Arial", regular = "/Library/Fonts/Arial.ttf")

ggplot2::ggsave(
  filename = file.path(outDir, "Figure1.pdf"),
  plot = fig1, dpi = 300, height = 12, width = 18, units = "in"
)

showtext_auto(FALSE)
