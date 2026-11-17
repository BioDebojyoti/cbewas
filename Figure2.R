# ──────────────────────────────────────────────────────────────────────────────
# Figure 2: gene–exposure overlap (chord diagram, UpSet)
# Output: Figure2.pdf                (a: chord, b: UpSet)
#         Supplementary_Figure2.pdf  (alluvial only)
#         Figure2_option_c.pdf       (a: alluvial, b: chord, c: UpSet)
# ──────────────────────────────────────────────────────────────────────────────
source("config.R")
library(ggalluvial)
library(circlize)
library(ComplexHeatmap)
library(grid)
load(resultsFile)

# ── Alluvial plot (panel_a) ───────────────────────────────────────────────────
alluvial_data <- gene_freq %>%
  tidyr::separate_rows(exposures, sep = ";\\s*") %>%
  dplyr::select(overlapping.genes, exposures) %>%
  dplyr::mutate(exposures = stringr::str_to_sentence(gsub(
    "\\.", " ", gsub(
      "\\.during.pregnancy", "",
      gsub("severe.life", "severe life", gsub(
        "potatoe", "potato",
        gsub("analgesics.med", "analgesics", exposures)
      ))
    )
  ))) %>%
  dplyr::count(overlapping.genes, exposures)

panel_a <- alluvial_data %>%
  dplyr::group_by(exposures) %>%
  dplyr::mutate(total_n = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(total_n), exposures) %>%
  ggplot(aes(y = n, axis1 = exposures, axis2 = overlapping.genes)) +
  geom_alluvium(aes(fill = exposures), width = 1/16, alpha = 0.7, decreasing = FALSE) +
  geom_stratum(width = 1/16, color = "grey60", decreasing = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            angle = 0, size = 2.0, decreasing = FALSE) +
  scale_x_discrete(limits = c("Exposures", "Genes"), expand = c(.05, .05)) +
  theme_minimal() +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = 12, face = "bold"),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot2::ggsave(
  filename = file.path(outDir, "Supplementary_Figure2.pdf"),
  plot = panel_a, dpi = 300, height = 12, width = 18, units = "in"
)

# ── Chord diagram data: top 7 genes shared by >1 exposure ─────────────────────
top_genes <- gene_freq %>%
  dplyr::filter(n_exposures > 1) %>%
  dplyr::slice_max(n_exposures, n = 7, with_ties = FALSE)

chord_data <- top_genes %>%
  tidyr::separate_rows(exposures, sep = ";\\s*") %>%
  dplyr::select(overlapping.genes, exposures) %>%
  dplyr::mutate(exposures = stringr::str_to_sentence(gsub(
    "\\.", "\n", gsub(
      "\\.during.pregnancy", "",
      gsub("severe.life", "severe life", gsub(
        "potatoe", "potato",
        gsub("analgesics.med", "analgesics", exposures)
      ))
    )
  )))

print(chord_data)

genes_in_chord     <- unique(chord_data$overlapping.genes)
exposures_in_chord <- unique(chord_data$exposures)
grid.col <- c(
  setNames(rainbow(length(genes_in_chord)), genes_in_chord),          # genes
  setNames(rep("grey80", length(exposures_in_chord)), exposures_in_chord)  # exposures
)

# ── UpSet data ────────────────────────────────────────────────────────────────
upset_matrix <- gene_freq %>%
  dplyr::arrange(n_exposures) %>%
  tidyr::separate_rows(exposures, sep = ";\\s*") %>%
  dplyr::mutate(exposures = gsub("\\.during.pregnancy", "", exposures)) %>%
  dplyr::distinct(overlapping.genes, exposures) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::pivot_wider(names_from = exposures, values_from = value, values_fill = 0) %>%
  tibble::column_to_rownames("overlapping.genes") %>%
  as.matrix()

comb_mat <- make_comb_mat(upset_matrix)

# ══════════════════════════════════════════════════════════════════════════════
# Figure2.pdf — a: alluvial (top), b: chord (bottom-left), c: UpSet (bottom-right)
# ══════════════════════════════════════════════════════════════════════════════
pdf(file = file.path(outDir, "Figure2_b.pdf"), width = 16, height = 16)

# Panel b: base graphics (chord diagram), bounded to bottom-left quadrant
par(fig = c(0, 0.5, 0, 0.5), mar = c(2, 2, 2, 2))
par(cex = 0.8)
chordDiagram(x = chord_data, grid.col = grid.col, transparency = 0.4,
             annotationTrack = c("name", "grid"))
circos.clear()

# Panel a: ggplot via grid, top half
pushViewport(viewport(x = 0.5, y = 0.75, width = 1.0, height = 0.5))
grid.draw(ggplotGrob(panel_a))
popViewport()

# Panel c: ComplexHeatmap UpSet, bottom-right quadrant
pushViewport(viewport(x = 0.75, y = 0.25, width = 0.48, height = 0.45))
draw(UpSet(comb_mat), newpage = FALSE)   # newpage = FALSE keeps earlier panels
popViewport()

# Panel labels
grid.text("a", x = 0.02, y = 0.95, gp = gpar(fontsize = 20, fontface = "bold"))
grid.text("b", x = 0.02, y = 0.47, gp = gpar(fontsize = 20, fontface = "bold"))
grid.text("c", x = 0.52, y = 0.47, gp = gpar(fontsize = 20, fontface = "bold"))

dev.off()

# ══════════════════════════════════════════════════════════════════════════════
# Figure2.pdf — a: chord (left), b: UpSet (right)
# ══════════════════════════════════════════════════════════════════════════════
pdf(file = file.path(outDir, "Figure2.pdf"), width = 16, height = 8)

par(fig = c(0, 0.5, 0, 1.0), mar = c(2, 2, 2, 2))
par(cex = 0.8)
chordDiagram(x = chord_data, grid.col = grid.col, transparency = 0.4,
             annotationTrack = c("name", "grid"))
circos.clear()

pushViewport(viewport(x = 0.75, y = 0.50, width = 0.48, height = 0.95))
draw(UpSet(comb_mat), newpage = FALSE)
popViewport()

grid.text("a", x = 0.02, y = 0.95, gp = gpar(fontsize = 20, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.95, gp = gpar(fontsize = 20, fontface = "bold"))

dev.off()
