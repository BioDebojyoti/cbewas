# ──────────────────────────────────────────────────────────────────────────────
# Table 2 — filtered DMRs across exposures
# Output: Table2.docx
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
load(resultsFile)

table2 <- all_dmr_results %>%
  dplyr::bind_rows(.id = "Exposure") %>%
  dplyr::arrange(Exposure) %>%
  dplyr::mutate(
    Exposure = stringr::str_to_sentence(gsub(
      "\\.", " ", gsub(
        "\\.during.pregnancy", "",
        gsub("severe.life", "severe life", gsub(
          "potatoe", "potato",
          gsub("analgesics.med", "analgesics",
               gsub("x.rayed", "X-rayed", Exposure))
        ))
      )
    ))
  ) %>%
  dplyr::mutate(
    across(c("Stouffer", "HMFDR", "Fisher"),
           function(x) {
             dplyr::if_else(
               x > 0 & x < 0.001,
               formatC(x, digits = 3, format = "e"),
               formatC(x, digits = 3, format = "f")
             )
           })
  ) %>%
  dplyr::mutate(across(c("maxdiff", "interaction_deltabeta"), \(x) round(x, 3))) %>%
  dplyr::select(-c("min_smoothed_fdr", "delta_ctrl", "delta_t1d", "pattern", "direction")) %>%
  dplyr::mutate(
    Exposure = ifelse(row_number() != 1 & lag(Exposure) == Exposure, NA_character_, Exposure)
  )

table2_flex <- table2 %>%
  flextable::flextable() %>%
  flextable::set_header_labels(interaction_deltabeta = "interaction \u0394\u03B2") %>%
  flextable::set_table_properties(opts_word = list(repeat_headers = TRUE)) %>%
  flextable::autofit()

save_as_docx(
  table2_flex,
  path = file.path(outDir, "Table2.docx"),
  pr_section = sect_properties
)
