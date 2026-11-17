# ──────────────────────────────────────────────────────────────────────────────
# analysis.R — interaction EWAS / DMR analysis (T1D vs Control)
# Produces: analysis_results.RData (+ dmr_results_T1D_vs_Control.html)
# Run this first; Figure*/Table* scripts depend on its output.
# ──────────────────────────────────────────────────────────────────────────────
rm(list = ls())
set.seed(123)

source("config.R")

library(limma)
library(DMRcate)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Gviz)
library(tictoc)
library(kableExtra)

# -----------------------------
# LOAD DATA
# -----------------------------
load(file = paste0(dataDir, "/EPIC_limma_ready_objects.RData"))

Mmat <- combat_data_slide_plus_array
rm(pd, combat_data_slide_plus_array)
gc()

meta <- readRDS(metaFile)

# ── Restrict to T1D vs Control only ──────────────────────────────────────────
meta <- meta %>%
  dplyr::filter(Status %in% c("T1D", "Con")) %>%
  droplevels()

rownames(meta) <- meta$Sample_Name
Mmat <- Mmat[, rownames(meta)]

stopifnot(identical(colnames(Mmat), rownames(meta)))

# ----------- create corresponding beta-value matrix for downstream analysis ---
Bmat <- 2^Mmat / (2^Mmat + 1)

exposures <- colnames(meta)[29:83]

# ---- exposures to test ------
limma_compatible <- c()
maybe_grouped    <- c()

ncol_required <- length(levels(meta$Status))
samples_required_per_level <- 3

for (grp in exposures) {
  curr_tab <- table(meta[[grp]], meta$Status, useNA = "no")
  if (nrow(curr_tab) == 2 & ncol(curr_tab) == ncol_required & all(curr_tab >= samples_required_per_level)) {
    limma_compatible <- c(limma_compatible, grp)
  } else {
    if (nrow(curr_tab) != 2 || ncol(curr_tab) != ncol_required) {
      message(paste0(grp, " is not usable and shoud be dropped!"))
    } else {
      maybe_grouped <- c(maybe_grouped, grp)
    }
  }
}

# -----------------------------
# SET FACTORS
# -----------------------------
meta$Status <- factor(meta$Status)
meta$Status <- relevel(meta$Status, ref = "Con")

# ── Results collectors ────────────────────────────────────────────────────────
all_dmc_results <- list()
all_dmr_results <- list()
all_granges     <- list()

ewas_log <- tibble(
  exposure           = character(),
  n_T1D_exp          = integer(),
  n_T1D_unexp        = integer(),
  n_ctrl_exp         = integer(),
  n_ctrl_unexp       = integer(),
  n_nominal_p05      = integer(),
  n_fdr_p05          = integer(),
  ran_dmrcate        = logical(),
  n_dmr_raw          = integer(),
  n_dmr_filtered     = integer()
)

# ── One-time setup: EPIC annotation probe ID → chr:pos mapping ───────────────
epic_anno <- as.data.frame(
  getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
) %>%
  tibble::rownames_to_column("cpg") %>%
  dplyr::mutate(chrpos = paste0(chr, ":", pos)) %>%
  dplyr::select(cpg, chrpos)

cat(sprintf("EPIC annotation loaded: %d probes\n", nrow(epic_anno)))

# ══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP over exposures
# ══════════════════════════════════════════════════════════════════════════════
tic("limma step: ")

for (x in limma_compatible) {
  
  cat(sprintf("\n\n---\n\n## Exposure: %s\n\n", x))
  
  # ── Drop samples with NA in this exposure ──────────────────────────────────
  keep      <- !is.na(meta[[x]])
  meta_loop <- meta[keep, ]
  Mmat_loop <- Mmat[, keep]
  Bmat_loop <- Bmat[, keep]
  
  # 1. DYNAMIC EXPOSURE REFERENCE
  e_levels <- unique(as.character(meta_loop[[x]]))
  if (length(e_levels) != 2) next
  
  if ("no" %in% e_levels) {
    e_ref <- "no"
  } else if ("low" %in% e_levels) {
    e_ref <- "low"
  } else if ("no milk" %in% e_levels) {
    e_ref <- "no milk"
  } else {
    e_ref <- e_levels[1]
  }
  meta_loop[[x]] <- relevel(factor(meta_loop[[x]]), ref = e_ref)
  
  # ── Sample counts per 2×2 cell ─────────────────────────────────────────────
  ref_status <- levels(meta_loop$Status)[1]
  alt_status <- levels(meta_loop$Status)[2]
  
  ref_exposure <- levels(meta_loop[[x]])[1]
  alt_exposure <- levels(meta_loop[[x]])[2]
  
  cell_counts <- table(meta_loop$Status, meta_loop[[x]])
  
  n_t1d_exp    <- cell_counts[alt_status, alt_exposure]
  n_t1d_unexp  <- cell_counts[alt_status, ref_exposure]
  n_ctrl_exp   <- cell_counts[ref_status, alt_exposure]
  n_ctrl_unexp <- cell_counts[ref_status, ref_exposure]
  
  cat(sprintf("  T1D:     exposed=%d  unexposed=%d\n", n_t1d_exp,  n_t1d_unexp))
  cat(sprintf("  Control: exposed=%d  unexposed=%d\n", n_ctrl_exp, n_ctrl_unexp))
  
  if (any(c(n_t1d_exp, n_t1d_unexp, n_ctrl_exp, n_ctrl_unexp) < 1L)) {
    cat("  Empty cell in 2x2 grid. Skipping.\n")
    next
  }
  
  # ── Build design ───────────────────────────────────────────────────────────
  design_formula <- as.formula(paste0(
    "~ Status * ", x,
    " + Sex + HLA_group + CD8T + CD4T + NK + Bcell + Mono + Gran"
  ))
  
  design <- tryCatch(
    model.matrix(design_formula, data = meta_loop),
    error = function(e) { message("  model.matrix failed: ", e$message); NULL }
  )
  if (is.null(design)) next
  
  all_cols <- colnames(design)
  coef_interaction <- all_cols[grep(":", all_cols)]
  exposure_col     <- stringr::str_split(coef_interaction, ":")[[1]][2]
  
  if (length(coef_interaction) != 1L) {
    cat("  Could not identify interaction column. Skipping.\n")
    next
  }
  cat(sprintf("  Interaction coef: %s\n", coef_interaction))
  
  # FIT 1 — M-values: robust inference (Du et al. 2010)
  fit_m  <- lmFit(Mmat_loop, design)
  fit_m  <- eBayes(fit_m, robust = TRUE)
  tt_m   <- topTable(fit_m, coef = coef_interaction,
                     number = Inf, sort.by = "none")
  
  # FIT 2 — β values: effect sizes only (do NOT use its p-values)
  fit_b  <- lmFit(Bmat_loop, design)
  fit_b  <- eBayes(fit_b, robust = TRUE)
  coefs_b <- coef(fit_b)
  
  # Cell-mean algebra (Status ref = Control, exposure ref = No):
  #   delta_ctrl = exposure_col ; delta_T1D = exposure_col + interaction
  delta_ctrl            <- coefs_b[, exposure_col]
  delta_t1d             <- coefs_b[, exposure_col] + coefs_b[, coef_interaction]
  interaction_deltabeta <- coefs_b[, coef_interaction]
  
  stopifnot(identical(rownames(tt_m), rownames(coefs_b)))
  
  tt_final <- tt_m %>%
    dplyr::mutate(
      delta_ctrl            = delta_ctrl,
      delta_t1d             = delta_t1d,
      interaction_deltabeta = interaction_deltabeta
    )
  
  all_dmc_results[[x]] <- tt_final
  
  n_nom <- sum(tt_final$P.Value < 0.05, na.rm = TRUE)
  cat(sprintf("  Nominal CpGs (p<0.05): %d\n", n_nom))
  
  n_nom_interaction_deltabeta <- sum(tt_final$P.Value < 0.05 & abs(tt_final$interaction_deltabeta) > DELTA_ITX_THRESH, na.rm = TRUE)
  cat(sprintf("  Nominal CpGs (p<0.05) & abs(interaction Δβ) > %3.2f : %d\n", DELTA_ITX_THRESH, n_nom_interaction_deltabeta))
  
  n_fdr <- sum(tt_final$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("  FDR-corrected CpGs (FDR<0.05): %d\n", n_fdr))
  
  log_row <- tibble(
    exposure                     = x,
    n_T1D_exp                    = n_t1d_exp,
    n_T1D_unexp                  = n_t1d_unexp,
    n_ctrl_exp                   = n_ctrl_exp,
    n_ctrl_unexp                 = n_ctrl_unexp,
    n_nominal_p05                = n_nom,
    n_nom_interaction_deltabeta  = n_nom_interaction_deltabeta,
    n_fdr_p05                    = n_fdr,
    ran_dmrcate                  = FALSE,
    n_dmr_raw                    = 0L,
    n_dmr_filtered               = 0L
  )
  
  if (n_nom < MIN_NOMINAL_CPGS) {
    cat("  Signal too weak for regional analysis. DMRcate skipped.\n")
    ewas_log <- bind_rows(ewas_log, log_row)
    next
  }
  
  # ── DMRcate — FULL M-value fit (no CpG subsetting) ─────────────────────────
  anno <- suppressWarnings(suppressMessages(
    cpg.annotate(
      object        = Mmat_loop,
      datatype      = "array",
      what          = "M",
      analysis.type = "differential",
      design        = design,
      coef          = coef_interaction,
      contrasts     = FALSE,
      arraytype     = "EPICv1",
      fdr           = 0.05
    )
  ))
  
  dmr <- tryCatch({
    dmrcate(
      anno,
      lambda     = 1000,
      C          = 2,
      pcutoff    = PCUT_DMRCATE,
      betacutoff = DELTA_ITX_THRESH,
      min.cpgs   = 3
    )
  }, error = function(e) {
    message("Skipping DMR: ", x)
    return(NULL)
  })
  
  if (is.null(dmr)) {
    dmr_result <- NULL
  } else {
    dmr_result <- tryCatch(
      suppressWarnings(suppressMessages(
        extractRanges(dmr, genome = "hg19")
      )),
      error = function(e) {
        message("  extractRanges error (likely empty DMR object): ", e$message)
        NULL
      }
    )
  }
  
  log_row$ran_dmrcate <- TRUE
  
  if (is.null(dmr_result) || length(dmr_result) == 0L) {
    cat("  DMRcate returned no regions.\n")
    ewas_log <- bind_rows(ewas_log, log_row)
    next
  } else {
    all_granges[[x]] <- dmr_result
  }
  
  dmr_df <- as.data.frame(dmr_result) %>%
    dplyr::rename(interaction_deltabeta = meandiff)
  
  log_row$n_dmr_raw <- nrow(dmr_df)
  
  # ── Attach region-level Δβ values via overlap with CpG-level results ──────
  anno_df <- as.data.frame(anno@ranges) %>%
    tibble::rownames_to_column("chrpos")
  
  deltabeta_df <- tt_final %>%
    tibble::rownames_to_column("cpg") %>%
    dplyr::select(cpg, delta_ctrl, delta_t1d, interaction_deltabeta) %>%
    dplyr::left_join(epic_anno, by = "cpg")
  
  anno_df_with_db <- anno_df %>%
    dplyr::left_join(
      deltabeta_df %>% dplyr::select(chrpos, delta_ctrl, delta_t1d, interaction_deltabeta),
      by = "chrpos"
    )
  
  stopifnot(nrow(anno_df_with_db) == nrow(anno_df))
  cat(sprintf("CpGs with Δβ attached: %d / %d\n",
              sum(!is.na(anno_df_with_db$interaction_deltabeta)),
              nrow(anno_df_with_db)))
  
  cpg_gr <- makeGRangesFromDataFrame(
    anno_df_with_db,
    seqnames.field     = "seqnames",
    start.field        = "start",
    end.field          = "end",
    keep.extra.columns = TRUE
  )
  
  ol <- findOverlaps(dmr_result, cpg_gr)
  
  region_stats <- tibble(
    dmr_idx    = queryHits(ol),
    delta_ctrl = cpg_gr$delta_ctrl[subjectHits(ol)],
    delta_t1d  = cpg_gr$delta_t1d[subjectHits(ol)]
  ) %>%
    dplyr::group_by(dmr_idx) %>%
    dplyr::summarise(
      across(everything(), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  dmr_df$delta_ctrl <- NA_real_
  dmr_df$delta_t1d  <- NA_real_
  
  dmr_df$delta_ctrl[region_stats$dmr_idx] <- region_stats$delta_ctrl
  dmr_df$delta_t1d[region_stats$dmr_idx]  <- region_stats$delta_t1d
  
  # ── Filter DMRs ────────────────────────────────────────────────────────────
  dmr_filtered <- dmr_df %>%
    dplyr::filter(
      HMFDR < HMFDR_THRESH,
      abs(interaction_deltabeta) > DELTA_ITX_THRESH
    ) %>%
    dplyr::mutate(
      pattern = case_when(
        sign(delta_t1d) == sign(delta_ctrl) & abs(delta_t1d) > abs(delta_ctrl) ~ "amplification",
        sign(delta_t1d) == sign(delta_ctrl) & abs(delta_t1d) < abs(delta_ctrl) ~ "attenuation",
        sign(delta_t1d) != sign(delta_ctrl)                                     ~ "reversal",
        TRUE                                                                     ~ "indeterminate"
      ),
      direction = case_when(
        pattern == "amplification" & delta_t1d > 0 ~ "both hyper; T1D > Control",
        pattern == "amplification" & delta_t1d < 0 ~ "both hypo; T1D > Control",
        pattern == "attenuation"   & delta_t1d > 0 ~ "both hyper; Control > T1D",
        pattern == "attenuation"   & delta_t1d < 0 ~ "both hypo; Control > T1D",
        pattern == "reversal" & delta_t1d > 0 ~ "T1D hyper / Control hypo",
        pattern == "reversal" & delta_t1d < 0 ~ "T1D hypo / Control hyper",
        TRUE                                   ~ "indeterminate"
      )
    )
  
  log_row$n_dmr_filtered <- nrow(dmr_filtered)
  ewas_log               <- bind_rows(ewas_log, log_row)
  
  if (nrow(dmr_filtered) > 0L) {
    all_dmr_results[[x]] <- dmr_filtered
  }
  
  cat(sprintf("  DMRs: raw=%d  filtered=%d\n",
              log_row$n_dmr_raw, log_row$n_dmr_filtered))
}

toc()

# ── Combined DMR table ────────────────────────────────────────────────────────
dmr_combined <- all_dmr_results %>%
  lapply(function(df) {
    df %>%
      mutate(
        overlapping.genes = sapply(overlapping.genes, function(g) {
          if (is.null(g) || all(is.na(g))) return(NA_character_)
          paste(unlist(g), collapse = "; ")
        })
      )
  }) %>% dplyr::bind_rows(.id = "exposure") %>%
  dplyr::select(
    exposure, seqnames, start, end, no.cpgs,
    interaction_deltabeta, delta_ctrl, delta_t1d,
    HMFDR, overlapping.genes,
    pattern, direction
  ) %>%
  dplyr::arrange(exposure, HMFDR) %>%
  dplyr::left_join(
    ewas_log %>%
      dplyr::mutate(
        min_cell = pmin(n_T1D_exp, n_T1D_unexp, n_ctrl_exp, n_ctrl_unexp)
      ) %>%
      dplyr::select(exposure, min_cell),
    by = "exposure"
  ) %>%
  dplyr::mutate(
    sample_distribution = case_when(
      min_cell < SAMPLE_SPARSE_THR   ~ "Sparse (Bias Risk)",
      min_cell < SAMPLE_MODERATE_THR ~ "Moderately Constrained (Variance Risk)",
      TRUE                           ~ "Adequate (Stable Cell Size)"
    )
  )

dmr_combined %>%
  dplyr::mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
  kable(caption = "DMRs across all exposures") %>%
  kable_styling(font_size = 10) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  scroll_box(width = "100%", height = "600px") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  save_kable(file.path(outDir, "dmr_results_T1D_vs_Control.html"))

# ── Genomic inflation factor (λ) per exposure (used by Table1.R) ──────────────
clean_name <- function(nm) gsub("\\.", " ", gsub("\\.during\\.pregnancy", "", nm))

interaction_limma_pvalue_uniform_distribution_metrics <- purrr::imap_dfr(
  all_dmc_results,
  function(df, nm) {
    p      <- df$P.Value[!is.na(df$P.Value)]
    chisq  <- qchisq(p, df = 1, lower.tail = FALSE)
    lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
    tibble(
      exposure = stringr::str_to_sentence(clean_name(nm)),
      lambda   = as.character(round(lambda, 3))
    )
  }
)

# ── Gene / exposure frequency metrics ─────────────────────────────────────────
gene_freq <- dmr_combined %>%
  tidyr::separate_rows(overlapping.genes, sep = ",\\s*") %>%
  dplyr::filter(!is.na(overlapping.genes), overlapping.genes != "NA",
                nchar(trimws(overlapping.genes)) > 0) %>%
  dplyr::mutate(overlapping.genes = trimws(overlapping.genes)) %>%
  dplyr::group_by(overlapping.genes) %>%
  dplyr::summarise(
    n_exposures = n_distinct(exposure),
    exposures   = paste(unique(exposure), collapse = "; "),
    .groups     = "drop"
  ) %>%
  dplyr::arrange(desc(n_exposures))

print(gene_freq)

exposure_freq <- dmr_combined %>%
  tidyr::separate_rows(overlapping.genes, sep = ",\\s*") %>%
  dplyr::filter(!is.na(overlapping.genes), overlapping.genes != "NA",
                nchar(trimws(overlapping.genes)) > 0) %>%
  dplyr::mutate(overlapping.genes = trimws(overlapping.genes)) %>%
  dplyr::group_by(exposure) %>%
  dplyr::summarise(
    n_genes = n_distinct(overlapping.genes),
    genes   = paste(unique(overlapping.genes), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_genes))

print(exposure_freq)

# ══════════════════════════════════════════════════════════════════════════════
# ZFP57 locus follow-up (feeds Figure2.R, Figure3.R, Supplementary_Figure2.R)
# ══════════════════════════════════════════════════════════════════════════════
get_gene_coordinates <- function(goi, dmr_combined) {
  tryCatch({
    filtered_df   <- dmr_combined %>% dplyr::filter(grepl(goi, overlapping.genes))
    unique_chroms <- unique(filtered_df$seqnames)
    if (length(unique_chroms) != 1) {
      message("Multiple chromosome found! Check data.")
      return(NULL)
    }
    filtered_df %>%
      dplyr::summarise(seqname = unique_chroms, start = min(start), end = max(end))
  }, error = function(e) NULL)
}

get_gene_cpgs <- function(gene_coordinates) {
  chr         <- gene_coordinates$seqname
  start_coord <- gene_coordinates$start
  stop_coord  <- gene_coordinates$end
  epic_anno %>%
    dplyr::filter(grepl(chr, chrpos)) %>%
    tidyr::separate(chrpos, into = c("cpg_chr", "cpg_pos"), sep = ":", convert = TRUE) %>%
    dplyr::filter(cpg_chr == chr, cpg_pos >= start_coord, cpg_pos <= stop_coord) %>%
    dplyr::pull(cpg)
}

get_gene_dmc_results <- function(goi_cpgs) {
  goi_results_list <- list()
  for (x in names(all_dmc_results)) {
    goi_results_list[[x]] <- all_dmc_results[[x]] %>%
      tibble::rownames_to_column("cpg") %>%
      dplyr::filter(cpg %in% goi_cpgs) %>%
      dplyr::select(all_of(c("cpg", "P.Value", "adj.P.Val",
                             "delta_ctrl", "delta_t1d", "interaction_deltabeta")))
  }
  dplyr::bind_rows(goi_results_list, .id = "exposure")
}

get_dmr_gene_dmc_summary <- function(goi_results) {
  goi_results %>%
    dplyr::group_by(exposure) %>%
    dplyr::summarise(
      n_cpgs_tested        = n(),
      n_nominal_p05        = sum(P.Value < 0.05,  na.rm = TRUE),
      n_nominal_p01        = sum(P.Value < 0.01,  na.rm = TRUE),
      n_nominal_p001       = sum(P.Value < 0.001, na.rm = TRUE),
      median_itx_deltabeta = median(interaction_deltabeta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(n_nominal_p001))
}

zfp57_dmr_coords <- get_gene_coordinates("ZFP57", dmr_combined)
zfp57_cpgs       <- get_gene_cpgs(zfp57_dmr_coords)
zfp57_results    <- get_gene_dmc_results(zfp57_cpgs)
zfp57_summary    <- get_dmr_gene_dmc_summary(zfp57_results)
print(zfp57_summary)

# Exposures in which ZFP57 DMR was found
common_exposures <- c(
  "infection.during.pregnancy",
  "worked.during.pregnancy",
  "severe.life.event.during.pregnancy",
  "antibiotics.during.pregnancy",
  "analgesics.med.during.pregnancy",
  "pets"
)

zfp57_cpg_consistency <- zfp57_results %>%
  dplyr::filter(exposure %in% common_exposures) %>%
  dplyr::group_by(cpg) %>%
  dplyr::summarise(
    n_exposures_nominal = sum(P.Value < 0.05, na.rm = TRUE),
    n_exposures_tested  = n(),
    mean_itx_deltabeta  = mean(interaction_deltabeta, na.rm = TRUE),
    sd_itx_deltabeta    = sd(interaction_deltabeta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_exposures_nominal))

print(zfp57_cpg_consistency, n = 22)

# Tier 1: nominal in ≥4 of 6 exposures; Tier 2: ≥3 of 6
core_cpgs_t1 <- zfp57_cpg_consistency %>% dplyr::filter(n_exposures_nominal >= 4) %>% dplyr::pull(cpg)
core_cpgs_t2 <- zfp57_cpg_consistency %>% dplyr::filter(n_exposures_nominal >= 3) %>% dplyr::pull(cpg)

cat(sprintf("Tier 1 core CpGs (≥4/6 exposures): %d\n", length(core_cpgs_t1)))
cat(sprintf("Tier 2 core CpGs (≥3/6 exposures): %d\n", length(core_cpgs_t2)))

# Keep only the ZFP57 β-values — the full Bmat is not needed downstream
# (zfp57_cpgs comes from the full EPIC annotation; Bmat only has probes that passed QC)
Bmat_zfp57 <- Bmat[intersect(zfp57_cpgs, rownames(Bmat)), , drop = FALSE]

# ══════════════════════════════════════════════════════════════════════════════
# SAVE everything the figure / table scripts need
# ══════════════════════════════════════════════════════════════════════════════
save(
  list = c(
    "meta", "ewas_log", "all_dmc_results", "all_dmr_results", "all_granges",
    "dmr_combined", "gene_freq", "exposure_freq", "epic_anno",
    "interaction_limma_pvalue_uniform_distribution_metrics",
    "zfp57_dmr_coords", "zfp57_cpgs", "zfp57_results", "zfp57_summary",
    "zfp57_cpg_consistency", "common_exposures", "core_cpgs_t1", "core_cpgs_t2",
    "Bmat_zfp57"
  ),
  file = resultsFile
)
cat("Saved:", resultsFile, "\n")
