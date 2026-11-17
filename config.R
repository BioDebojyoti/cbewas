# ──────────────────────────────────────────────────────────────────────────────
# config.R — shared paths, parameters and table settings
# Sourced by analysis.R and every Figure / Table script.
# Run all scripts with the working directory set to the folder containing them.
# ──────────────────────────────────────────────────────────────────────────────
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(purrr)
library(officer)
library(flextable)

# ── Paths ─────────────────────────────────────────────────────────────────────
dataDir <- "/Users/debda22/Projects/abis-epigenetic-patterns-at-birth/results-new/batch-correction-new"
metaFile <- "/Users/debda22/Projects/abis-epigenetic-patterns-at-birth/metadata/epigenetics.meta.data.RDS"
outDir  <- "/Users/debda22/Projects/abis-epigenetic-patterns-at-birth/results-test/interaction_dmc_and_dmr/T1D_vs_Con"
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

# Intermediate results written by analysis.R and read by all figure/table scripts
resultsFile <- file.path(outDir, "analysis_results.RData")

# ── Parameters ────────────────────────────────────────────────────────────────
PCUT_DMRCATE        <- 0.001  # DMRcate smoother input gate (not final significance)
HMFDR_THRESH        <- 0.05   # regional significance threshold
DELTA_ITX_THRESH    <- 0.10   # |interaction Δβ| threshold for DMR retention
MIN_NOMINAL_CPGS    <- 3      # skip DMRcate if fewer than this many CpGs at p<0.05
SAMPLE_SPARSE_THR   <- 4
SAMPLE_MODERATE_THR <- 8

# ── Word-table page settings (Table1.R, Table2.R) ─────────────────────────────
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape", height = 11.7, width = 8.3),
  page_margins = page_mar(top = 0.5, bottom = 0.5, left = 0.5, right = 0.5,
                          header = 0.3, footer = 0.3, gutter = 0),
  type = "continuous"
)

set_flextable_defaults(padding.top = 0, padding.bottom = 0, line_spacing = 1)
