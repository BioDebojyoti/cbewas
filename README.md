# Interaction EWAS / DMR analysis: T1D vs Control

Investigate whether the effect of each prenatal exposure on DNA methylation differs between
T1D cases and controls (`Status * exposure` interaction), at single-CpG level (limma)
and region level (DMRcate), and produces the figures and tables derived from it.

## Files

| File | Purpose | Outputs (in `outDir`) |
|---|---|---|
| `config.R` | Shared paths, thresholds, Word page settings. **Edit paths here.** | – |
| `analysis.R` | Full pipeline: limma + DMRcate loop, DMR table, λ per exposure, gene/exposure frequencies, ZFP57 follow-up | `analysis_results.RData`, `dmr_results_T1D_vs_Control.html` |
| `Figure1.R` | Interaction EWAS signal across exposures | `Figure1.pdf` |
| `Figure2.R` | Gene–exposure overlap: alluvial, chord, UpSet | `Figure2.pdf`, `Supplementary_Figure2.pdf` |
| `Figure3.R` | ZFP57 heatmap (interaction Δβ per CpG per exposure) | `Figure3.pdf` |
| `Figure4.R` | ZFP57 core-CpG box plots (infection during pregnancy) | `Figure4.pdf` |
| `Supplementary_Figure1.R` | Same as Figure 1 but only CpGs with nominal p < 0.05 (no \|Δβ\| stage) | `Supplementary_Figure1.pdf` |
| `Supplementary_Figure3.R` | ZFP57 core-CpG box plots across all common exposures | `Supplementary.Figure3.zfp57.across.exposures.pdf` |
| `Supplementary_Figure4.R` | PC-covariate association heatmaps before and after batch correction | `Supplementary_Figure4.pdf` |
| `Table1.R` | Sample distribution, nominal CpGs, DMRs, λ | `Table1.docx` |
| `Table2.R` | Filtered DMRs across exposures | `Table2.docx` |
| `helpers_zfp57_boxplot.R` | Box-plot function used by `Figure4.R` and `Supplementary_Figure3.R` | – |

`config.R` and `helpers_zfp57_boxplot.R` are sourced by the other scripts and are not run directly.
<!--| `Supplementary_Figure1.R` | P-value histograms and volcano plots per exposure | `Supplementary_Figure1_EWAS_PVALUE_HIST.pdf`, `Supplementary_Figure1_EWAS_PVALUE_HIST_AND_VOLCANO.pdf` | -->

## Supplementary Figure 3 (any exposure)

`Supplementary_Figure3.R` draws the Figure 4 style box plots (ZFP57 core CpGs) for any exposure.
Exposure names are the column names in `meta` (the same names used in `all_dmc_results`).

```sh
Rscript Supplementary_Figure3.R pets
Rscript Supplementary_Figure3.R pets infection.during.pregnancy   # several, one PDF each
```

or in an R session:

```r
exposure_to_plot <- "pets"
source("Supplementary_Figure3.R")
```

With no argument it uses `infection.during.pregnancy`. If the exposure has no ZFP57 DMR, the CpGs are
still plotted and the subtitle says no DMR was called. If it has several, the one with the lowest HMFDR is shown.
An unknown name stops with the list of analysed exposures. To plot another gene or CpG set, change `goi` and `cpgs` at the top of the script.

## Run order

1. `analysis.R` first. Everything else depends on `analysis_results.RData`.
2. Then any of the figure and table scripts, in any order (`Supplementary_Figure3.R` takes the exposure as an argument, see above).

Set the working directory to this folder (the scripts use `source("config.R")`).

```r
for (f in c("analysis.R",
            "Figure1.R", "Figure2.R", "Figure3.R", "Figure4.R",
            "Supplementary_Figure1.R",
            "Supplementary_Figure3.R",
            "Supplementary_Figure4.R",
            "Table1.R", "Table2.R")) {
  source(f)
}
```

`analysis.R` starts with `rm(list = ls())`, so it clears the session.

## Inputs

Paths are set in `config.R`:

- `dataDir/EPIC_limma_ready_objects.RData`: must contain `combat_data_slide_plus_array` (M-value matrix) and `pd`
- `metaFile`: sample metadata RDS (`Status`, `Sex`, `HLA_group`, cell-type proportions `CD8T CD4T NK Bcell Mono Gran`, `Sample_Name`, exposure columns 29–83)

## Method summary

- Samples restricted to `Status` in {T1D, Con}; control is the reference level.
- Binary exposures with at least 3 samples per cell level are tested (the exposure's reference level is `no`, `low`, `no milk`, or the first level).
- Model per exposure: `~ Status * exposure + Sex + HLA_group + CD8T + CD4T + NK + Bcell + Mono + Gran`.
- **M-values** give p-values (limma, `eBayes(robust = TRUE)`).
- **β-values** give effect sizes only (Δβ for controls, Δβ for T1D, interaction Δβ).
- **DMRcate** runs on the interaction coefficient (`lambda = 1000`, `C = 2`, `pcutoff = 0.001`, `min.cpgs = 3`). Regions are kept if HMFDR < 0.05 and |interaction Δβ| > 0.10.
- Each retained DMR is classified as amplification, attenuation or reversal from the T1D vs control Δβ.
- Exposures are flagged by smallest 2×2 cell size: sparse (< 4), constrained (< 8), adequate.

## Thresholds (in `config.R`)

| Parameter | Value | Meaning |
|---|---|---|
| `PCUT_DMRCATE` | 0.001 | DMRcate smoother input gate |
| `HMFDR_THRESH` | 0.05 | Regional significance |
| `DELTA_ITX_THRESH` | 0.10 | Minimum \|interaction Δβ\| |
| `MIN_NOMINAL_CPGS` | 3 | Skip DMRcate if fewer nominal CpGs |
| `SAMPLE_SPARSE_THR` / `SAMPLE_MODERATE_THR` | 4 / 8 | Cell-size cut-offs |

## `analysis_results.RData` contents

`meta`, `ewas_log`, `all_dmc_results`, `all_dmr_results`, `all_granges`, `dmr_combined`,
`gene_freq`, `exposure_freq`, `epic_anno`, `interaction_limma_pvalue_uniform_distribution_metrics`,
`zfp57_dmr_coords`, `zfp57_cpgs`, `zfp57_results`, `zfp57_summary`, `zfp57_cpg_consistency`,
`common_exposures`, `core_cpgs_t1`, `core_cpgs_t2`, `Bmat_zfp57`.

`all_dmc_results` holds genome-wide results for every exposure and is large (can be 1 GB or more).

## Requirements

R packages: limma, DMRcate, minfi, IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
IlluminaHumanMethylation450kanno.ilmn12.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene, Gviz,
dplyr, tidyr, tibble, stringr, purrr, ggplot2, ggrastr, patchwork, ggalluvial, circlize,
ComplexHeatmap, showtext, flextable, officer, kableExtra, tictoc.

`Figure1.R` registers Arial from `/Library/Fonts/Arial.ttf` (macOS path). Change it on other systems.

## Notes

- `common_exposures` (exposures with a ZFP57 DMR) is hardcoded in `analysis.R`.
- `Table1.R` removes columns by position (`delete_columns(j = c(11, 13))`). Update it if the column set changes.

## Authors
- Debojyoti Das ([BioDebojyoti](https://github.com/BioDebojyoti))
