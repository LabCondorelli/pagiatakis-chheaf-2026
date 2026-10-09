#!/usr/bin/env Rscript
## Fig. 1A - the 40 lincRNA modulated by TAC in at least one time point, as row Z-scores of
## log2 CPM averaged within each group, with the mean trajectory of each block.
##
## Rows are split by the direction of the differential calls (TAC vs sham at 1, 2, 4 and 7
## days, thresholded at FDR < 0.05 and |log2FC| > 0.5): up-regulated on top, down-regulated
## below, hierarchically clustered inside each block at the ComplexHeatmap default.
##
## The Z-score is computed once, over the 16 samples. The heatmap shows the mean of those
## Z-scores within each group of two replicates; the trajectories use them at sample level, so
## the ribbon carries the spread across lincRNAs.
##
## Usage: Rscript fig1A_lincRNA_heatmap.R   (run 00_export_mouse_tac.R first)
## Output: results/fig1A_lincRNA_heatmap.pdf, results/fig1A_lincRNA_zscores.tsv

suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

de  <- read.delim(file.path(OUT_DIR, "mouse_tac_DE_calls.tsv"), check.names = FALSE)
cpm <- read.delim(file.path(OUT_DIR, "mouse_tac_log2cpm.tsv"), check.names = FALSE, row.names = 1)

## ---- the 40 lincRNA ---------------------------------------------------------------------
GIORNI <- c("1d", "2d", "4d", "7d")
RES  <- paste0("Results.tac_", GIORNI, "_vs_sham_", GIORNI)
call <- as.matrix(de[, RES])
keep <- de$Genes.gene_biotype == "lincRNA" & rowSums(call != 0) > 0
sel  <- de[keep, ]; call <- call[keep, , drop = FALSE]

## Up-regulated when more time points call it up than down. Up first, then down.
direction <- factor(ifelse(rowSums(call == 1) > rowSums(call == -1), "up-regulated", "down-regulated"),
                    levels = c("up-regulated", "down-regulated"))
cat(sprintf("%d lincRNA: %d up-regulated, %d down-regulated\n", nrow(sel),
            sum(direction == "up-regulated"), sum(direction == "down-regulated")))

## ---- Z-scores over the 16 samples -------------------------------------------------------
smp  <- grep("^(Sham|TAC)[.]", colnames(cpm), value = TRUE)
m    <- as.matrix(cpm[match(sel$Genes.ensembl_gene_id, rownames(cpm)), smp])
rownames(m) <- sel$Genes.external_gene_name
z    <- t(scale(t(m)))

cond <- factor(sub("[.].*$", "", smp), levels = c("Sham", "TAC"))
day  <- factor(sub("_rep[AB]$", "", sub("^.*[.]", "", smp)), levels = GIORNI)
grp  <- interaction(cond, day, drop = TRUE, lex.order = TRUE)

## one column per group of two replicates
zg <- sapply(levels(grp), function(g) rowMeans(z[, grp == g, drop = FALSE]))
gc <- factor(sub("[.].*$", "", levels(grp)), levels = c("Sham", "TAC"))
gd <- factor(sub("^.*[.]", "", levels(grp)), levels = GIORNI)

## ---- trajectory of each block, drawn beside the heatmap ----------------------------------
PAL <- c(Sham = "#7FA8D4", TAC = "#E06A6A")

## Same smoothing as the original panel: a loess through every point, that is every lincRNA in
## every replicate, with the four time points as positions 1 to 4. These are the defaults of
## ggplot2::geom_smooth(se = TRUE) below 1,000 points, and they are what gives the curves their
## shoulder at day 4 and the width of the band. Smoothing the four group means instead would
## flatten both.
traiettoria <- function(righe, titolo) function(index, nm) {
  xi <- seq(1, 4, length.out = 80)
  grid.text(titolo, y = 0.95, just = "top", gp = gpar(fontsize = 7.5))
  pushViewport(viewport(x = 0.22, y = 0.13, width = 0.74, height = 0.72,
                        just = c("left", "bottom"), xscale = c(0.7, 4.3), yscale = c(-2.2, 2.2)))
  grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.8))
  for (cc in levels(cond)) {
    k  <- cond == cc
    df <- data.frame(x = rep(as.integer(day[k]), each = length(righe)),
                     y = as.numeric(z[righe, k]))
    ## With only four distinct time points the local quadratic is rank deficient at the ends
    ## and loess says so. The fit is still the one geom_smooth would produce, so the warning
    ## is silenced rather than the model changed.
    pr <- suppressWarnings({
      fit <- loess(y ~ x, data = df, span = 0.75, degree = 2)
      predict(fit, newdata = data.frame(x = xi), se = TRUE)
    })
    ci  <- qt(0.975, pr$df) * pr$se.fit
    grid.polygon(c(xi, rev(xi)), c(pr$fit - ci, rev(pr$fit + ci)), default.units = "native",
                 gp = gpar(fill = adjustcolor(PAL[[cc]], 0.30), col = NA))
    grid.lines(xi, pr$fit, default.units = "native", gp = gpar(col = PAL[[cc]], lwd = 1.6))
  }
  grid.yaxis(at = -2:2, gp = gpar(fontsize = 6))
  grid.xaxis(at = 1:4, label = GIORNI, gp = gpar(fontsize = 6))
  grid.text("Row Z-Score", x = unit(-2.4, "lines"), rot = 90, gp = gpar(fontsize = 6.5))
  popViewport()
}

blocchi <- split(seq_len(nrow(z)), direction)
etich <- sprintf("%d %s lincRNAs", lengths(blocchi), names(blocchi))

ha <- HeatmapAnnotation(
  Treatment = gc, `Days after surgery/banding` = gd,
  giorno = anno_text(sub("d$", "", as.character(gd)), just = "center",
                     rot = 0, location = unit(0.5, "npc"), gp = gpar(fontsize = 7)),
  col = list(Treatment = c(Sham = "#E8E8E8", TAC = "#9A9A9A"),
             `Days after surgery/banding` = c(`1d` = "#ECE7F2", `2d` = "#C9B3E0",
                                              `4d` = "#9E6FD0", `7d` = "#2222EE")),
  annotation_label = c("Treatment", "Days after surgery/banding", ""),
  annotation_name_gp = gpar(fontsize = 7.5), annotation_name_side = "right",
  show_legend = FALSE, simple_anno_size = unit(3.5, "mm"))

ht <- Heatmap(zg, name = "Row Z-Score",
              col = colorRamp2(c(-2, 0, 2), c("dodgerblue2", "grey20", "yellow2")),
              cluster_columns = FALSE, column_split = gc,
              row_split = direction, cluster_row_slices = FALSE,
              row_gap = unit(0, "mm"),           # the published panel has no gap between blocks
              row_title = NULL, column_title = c("Sham", "TAC"),
              column_title_gp = gpar(fontsize = 9),
              show_row_names = FALSE, show_column_names = FALSE,
              show_row_dend = TRUE, row_dend_width = unit(10, "mm"),
              top_annotation = ha,
              right_annotation = rowAnnotation(
                traiettorie = anno_zoom(align_to = blocchi, which = "row",
                                        panel_fun = function(index, nm)
                                          traiettoria(blocchi[[nm]], etich[match(nm, names(blocchi))])(index, nm),
                                        size = unit(32, "mm"), gap = unit(10, "mm"),
                                        width = unit(55, "mm"), link_width = unit(8, "mm"),
                                        link_gp = gpar(fill = "grey88", col = NA))),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 7.5),
                                          labels_gp = gpar(fontsize = 7)))

lg <- Legend(labels = levels(cond), title = NULL, type = "lines",
             legend_gp = gpar(col = PAL[levels(cond)], lwd = 2),
             labels_gp = gpar(fontsize = 7.5), background = "white")

pdf(file.path(OUT_DIR, "fig1A_lincRNA_heatmap.pdf"), width = 8.6, height = 5.4)
draw(ht, heatmap_legend_side = "left", annotation_legend_list = list(lg),
     annotation_legend_side = "left", padding = unit(c(2, 2, 2, 2), "mm"))
invisible(dev.off())

## written in the order of the panel: up-regulated block first, then down-regulated
o <- order(direction, rownames(zg))
write.table(data.frame(lincRNA = rownames(zg)[o], direction = as.character(direction)[o],
                       zg[o, , drop = FALSE], check.names = FALSE),
            file.path(OUT_DIR, "fig1A_lincRNA_zscores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
