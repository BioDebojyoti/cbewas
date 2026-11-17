# ──────────────────────────────────────────────────────────────────────────────
# Figure 3: ZFP57 heatmap — interaction Δβ per CpG per exposure
# Output: Figure3.pdf
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
load(resultsFile)

heatmap_df <- zfp57_results %>%
  dplyr::mutate(
    exposure_label = gsub("\\.", " ", gsub("\\.during\\.pregnancy", "", exposure)),
    nominal        = P.Value < 0.05,
    sig_label      = case_when(
      P.Value < 0.01 ~ "**",
      P.Value < 0.05 ~ "*",
      TRUE           ~ "")
  ) %>%
  dplyr::left_join(
    zfp57_cpg_consistency %>% dplyr::select(cpg, n_exposures_nominal),
    by = "cpg"
  ) %>%
  # Order CpGs by consistency then effect size
  dplyr::mutate(
    cpg = factor(cpg,
                 levels = zfp57_cpg_consistency %>%
                   dplyr::arrange(desc(n_exposures_nominal), desc(mean_itx_deltabeta)) %>%
                   dplyr::pull(cpg))
  )

fig2 <- heatmap_df %>%
  ggplot(aes(x = exposure_label, y = cpg, fill = interaction_deltabeta)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 4, colour = "white", fontface = "bold") +
  # Highlight tier 1 CpGs with a border
  geom_tile(
    data = heatmap_df %>% dplyr::filter(cpg %in% core_cpgs_t1 & exposure %in% common_exposures),
    colour = "black", linewidth = 1.0, fill = NA
  ) +
  scale_fill_gradientn(
    colours  = c("#2C7BB6", "white", "#D7191C"),
    values   = scales::rescale(c(-0.224, 0, 0.309)),
    limits   = c(-0.224, 0.309),
    na.value = "grey80",
    name     = expression("Interaction " * Delta * beta)
  ) +
  scale_x_discrete(position = "top") +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(size = 18, face = "bold"),
    axis.text.x     = element_text(size = 12, face = "bold", angle = 30, hjust = 0),
    axis.text.y     = element_text(size = 12, face = "bold", family = "mono"),
    axis.title.y    = element_text(size = 12, face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "right",
    plot.caption    = element_text(size = 9, face = "bold")
  ) +
  labs(
    x        = NULL,
    y        = "CpG site",
    title    = expression("ZFP57 — Interaction " * Delta * beta * " per CpG per exposure"),
    subtitle = "* p<0.05  ** p<0.01",
    caption  = "Black border = CpGs significant (nominal p-value) in ≥ 4/6 exposures"
  )

fig2

ggplot2::ggsave(
  filename = file.path(outDir, "Figure3.pdf"),
  plot = fig2, dpi = 300, height = 12, width = 18, units = "in"
)
