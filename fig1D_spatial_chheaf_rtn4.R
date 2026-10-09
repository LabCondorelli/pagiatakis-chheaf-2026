#!/usr/bin/env Rscript
## Fig. 1D - H&E, Rtn4 (measured), Chheaf (Tangram) and their co-expression on four sections.
## Rtn4 and Chheaf scaled per section; co-expression on the Fig. S1F scale, cap 2.6.
## Usage: Rscript fig1D_spatial_chheaf_rtn4.R
## Output: results/fig1D_spatial_chheaf_rtn4.{pdf,png}, results/fig1D_scales.tsv

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork); library(grid) })
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "00_config.R"))
source(file.path(here, "utils", "spatial_helpers.R"))
options(future.globals.maxSize = 1e16)

PANEL_IN <- c(3.56, 5.38)
PT <- c(GSM5355663_WT_Sham = 0.215, GSM5943190_WT_MI_day1_2 = 0.211,
        GSM5943192_WT_MI_day7_2 = 0.241, GSM5943195_WT_MI_day14_3 = 0.209)
COLUMNS <- c("Sham", "Day1", "Day7", "Day14")
ROWS <- c(he = "H&E Image", rtn4 = "Rtn4", chheaf = "Chheaf", coloc = "Rtn4 + Chheaf")
LAB_FRAC <- 0.085
PAL <- viridis::viridis(11)

vis <- subset(readRDS(VISIUM_RDS), orig.ident %in% PUB_SLIDES)
vis <- NormalizeData(vis, assay = "Spatial", verbose = FALSE)
DefaultAssay(vis) <- "Spatial"
for (sl in PUB_SLIDES) vis <- rotate_section(vis, sl, SLIDE_ROTATION[[sl]])

proj <- read.csv(TANGRAM_GENES, row.names = 1, check.names = FALSE)
if (!all(colnames(vis) %in% rownames(proj))) stop("Tangram spot ids do not match the Visium object")
proj <- proj[colnames(vis), ]
expr <- GetAssayData(vis, assay = "Spatial", layer = "data")

## z-scores over the four sections pooled.
coloc <- lapply(seq_len(nrow(PAIRS)), function(i)
   coexpression_score(proj[[tolower(PAIRS$linc[i])]], expr[PAIRS$gene[i], ]))
pooled_min <- min(unlist(coloc))

vis$rtn4   <- expr["Rtn4", ]
vis$chheaf <- proj[["gm12092"]]
vis$coloc  <- coloc[[which(PAIRS$linc == "Gm12092" & PAIRS$gene == "Rtn4")]]

spot_row <- function(feature, lims = NULL) wrap_plots(lapply(PUB_SLIDES, function(sl)
   SpatialFeaturePlot(vis, features = feature, images = sl, pt.size.factor = PT[[sl]],
                      image.alpha = 0, stroke = 0) +
      scale_fill_gradientn(colours = PAL, limits = lims, oob = scales::squish) +
      theme_void() + theme(legend.position = "none", plot.margin = margin(1, 1, 1, 1))), nrow = 1)

he_panel <- function(sl) {
   bc <- sub("_\\d+$", "", colnames(vis)[vis$orig.ident == sl])
   img <- as.raster(he_crop(file.path(RAW_DIR, sl), bc, SLIDE_ROTATION[[sl]]))
   cell <- grid.layout(1, 1, widths = unit(ncol(img) / nrow(img), "null"), heights = unit(1, "null"), respect = TRUE)
   wrap_elements(full = gTree(children = gList(rasterGrob(img, interpolate = TRUE)),
                              vp = vpStack(viewport(layout = cell), viewport(layout.pos.row = 1, layout.pos.col = 1))))
}

row_label <- function(txt) wrap_elements(full = gTree(children = gList(
   linesGrob(x = c(0.82, 0.82), y = c(0.04, 0.96), gp = gpar(lwd = 1.1)),
   textGrob(txt, x = 0.45, rot = 90, gp = gpar(fontsize = 7)))))
with_label <- function(key, body) wrap_plots(list(row_label(ROWS[[key]]), body), ncol = 2,
                                             widths = c(LAB_FRAC, 1 - LAB_FRAC))

fig <- wrap_plots(list(
   wrap_elements(full = textGrob("D", x = 0.005, y = 0.6, hjust = 0, gp = gpar(fontsize = 13))),
   wrap_plots(list(plot_spacer(), wrap_elements(full = column_header_grob(COLUMNS))), ncol = 2,
              widths = c(LAB_FRAC, 1 - LAB_FRAC)),
   with_label("he", wrap_plots(lapply(PUB_SLIDES, he_panel), nrow = 1)),
   with_label("rtn4", spot_row("rtn4")),
   with_label("chheaf", spot_row("chheaf")),
   with_label("coloc", spot_row("coloc", c(pooled_min, COEXPR_CAP))),
   wrap_elements(full = colourbar_grob(PAL, "Expression", bar_width_in = 0.62, title_side = "bottom"))),
   ncol = 1, heights = c(0.30, 0.32, 1, 1, 1, 1, 0.55)) +
   plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

out <- file.path(OUT_DIR, "fig1D_spatial_chheaf_rtn4")
ggsave(paste0(out, ".pdf"), fig, width = PANEL_IN[1], height = PANEL_IN[2])
ggsave(paste0(out, ".png"), fig, width = PANEL_IN[1], height = PANEL_IN[2], dpi = 400)

scales_tab <- do.call(rbind, lapply(PUB_SLIDES, function(sl) {
   i <- vis$orig.ident == sl
   data.frame(section = sl, spots = sum(i),
              rtn4_min = min(vis$rtn4[i]), rtn4_max = max(vis$rtn4[i]),
              chheaf_min = min(vis$chheaf[i]), chheaf_max = max(vis$chheaf[i]),
              coloc_max = max(vis$coloc[i]), coloc_pct_above_cap = 100 * mean(vis$coloc[i] > COEXPR_CAP))
}))
scales_tab$coloc_scale <- sprintf("%.3f-%.1f", pooled_min, COEXPR_CAP)
scales_tab[] <- lapply(scales_tab, function(x) if (is.double(x)) signif(x, 4) else x)
write.table(scales_tab, file.path(OUT_DIR, "fig1D_scales.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(scales_tab, row.names = FALSE)
cat(sprintf("co-expression scale %.3f - %.1f; written %s.{pdf,png}\n", pooled_min, COEXPR_CAP, out))
