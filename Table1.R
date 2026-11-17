# ──────────────────────────────────────────────────────────────────────────────
# Table 1 — sample distribution and epigenetic signal across exposures
# Output: Table1.docx
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
load(resultsFile)

table1 <- ewas_log %>%
  dplyr::select(all_of(c(
    "exposure", "n_T1D_exp", "n_T1D_unexp", "n_ctrl_exp", "n_ctrl_unexp",
    "n_nominal_p05", "n_dmr_filtered"
  ))) %>%
  dplyr::mutate(
    exposure = stringr::str_to_sentence(
      gsub("\\.", " ", gsub("\\.during\\.pregnancy", "", exposure))
    ),
    min_cell = pmin(n_T1D_exp, n_T1D_unexp, n_ctrl_exp, n_ctrl_unexp),
    sample_distribution = case_when(
      min_cell < SAMPLE_SPARSE_THR   ~ "Sparse (bias risk)",
      min_cell < SAMPLE_MODERATE_THR ~ "Constrained (variance risk)",
      TRUE                           ~ "Adequate"
    ),
    sample_distribution = factor(
      sample_distribution,
      levels = c("Adequate", "Constrained (variance risk)", "Sparse (bias risk)")
    )
  ) %>%
  dplyr::arrange(sample_distribution) %>%
  dplyr::left_join(
    interaction_limma_pvalue_uniform_distribution_metrics %>%
      dplyr::select(all_of(c("exposure", "lambda"))),
    by = "exposure"
  ) %>%
  dplyr::rename(
    Exposure       = exposure,
    `T1D (exp)`    = n_T1D_exp,
    `T1D (unexp)`  = n_T1D_unexp,
    `Ctrl (exp)`   = n_ctrl_exp,
    `Ctrl (unexp)` = n_ctrl_unexp,
    `CpGs†`        = n_nominal_p05,
    `DMRs‡`        = n_dmr_filtered,
    `Min. n`       = min_cell,
    `Cell quality` = sample_distribution,
    "\u03BB"       = lambda
  )

# Names must match the `status` strings below exactly (including the symbol prefix)
color_map <- c(
  "⚠ Notable inflation in ≥1 exposure" = "#AA1111",  # Dark Red
  "~ Mild inflation in ≥1 exposure"    = "#D97706",  # Amber/Orange
  "↓ Deflated across exposures"        = "#2563EB",  # Blue
  "✓ Well-calibrated"                  = "#16A34A"   # Green
)

table1_flex <- table1 %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    lambda_values = list(as.numeric(c_across(contains("\u03BB")))),
    status = case_when(
      any(unlist(lambda_values) > 1.20, na.rm = TRUE) ~ "⚠ Notable inflation in ≥1 exposure",
      any(unlist(lambda_values) > 1.10, na.rm = TRUE) ~ "~ Mild inflation in ≥1 exposure",
      all(unlist(lambda_values) < 0.90, na.rm = TRUE) ~ "↓ Deflated across exposures",
      TRUE                                            ~ "✓ Well-calibrated"
    ),
    row_color = unname(color_map[status])
  ) %>%
  dplyr::ungroup() %>%
  flextable::flextable() %>%
  flextable::set_caption(
    caption = paste0(
      "Sample distribution and epigenetic signal across prenatal exposures. ",
      "E = exposed, U = unexposed."
    )
  ) %>%
  flextable::align(j = 1, align = "left", part = "body") %>%
  flextable::align(j = 2:9, align = "center", part = "body") %>%
  flextable::align(align = "center", part = "header") %>%
  flextable::padding(padding.top = 4, padding.bottom = 4, part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 11, part = "all") %>%
  flextable::bold(j = 9, bold = FALSE, part = "body") %>%
  flextable::color(
    i = ~ `Cell quality` == "Sparse (bias risk)",
    color = "grey60",
    part = "body"
  ) %>%
  flextable::color(
    j = "λ",
    color = .$body$dataset$row_color,
    part = "body"
  ) %>%
  flextable::footnote(
    i = 1, j = 8,
    value = as_paragraph("Nominally significant CpGs (p < 0.05) for the interaction term StatusT1D X Exposure."),
    ref_symbols = "*",
    part = "header"
  ) %>%
  flextable::footnote(
    i = 1, j = 9,
    value = as_paragraph(paste0("DMRs meeting HMFDR < 0.05, | interaction Δβ| > ", DELTA_ITX_THRESH, ", and ≥3 constituent CpGs.")),
    ref_symbols = "†",
    part = "header"
  ) %>%
  flextable::autofit()

save_as_docx(
  table1_flex %>% flextable::delete_columns(j = c(11, 13)),
  path = file.path(outDir, "Table1.docx"),
  pr_section = sect_properties
)
