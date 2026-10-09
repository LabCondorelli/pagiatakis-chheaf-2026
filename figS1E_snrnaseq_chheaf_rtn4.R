#!/usr/bin/env Rscript
## Fig. S1E - snRNA-seq UMAP, WGCNA module scores per cell type (z-scored), Nebulosa densities
## of Rtn4, Chheaf and their product, overall and per timepoint.
## Usage: Rscript figS1E_snrnaseq_chheaf_rtn4.R
## Output: results/figS1E_snrnaseq_chheaf_rtn4.{pdf,png}, results/figS1E_module_score_zscore.tsv

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))
options(future.globals.maxSize = 1e16)

PANEL_IN <- c(2.95, 4.5)
GENES <- c(Rtn4 = "Rtn4", Chheaf = "Gm12092")
LIM <- list(rtn4 = c(0, 0.01), chheaf = c(0, 0.045), joint = c(0, 3e-4), joint_tp = c(0, 7e-4))
CT <- c(Cardiomyocyte_1 = "Cardiomyocyte-1", Cardiomyocyte_2 = "Cardiomyocyte-2",
        Cardiomyocyte_3 = "Cardiomyocyte-3", Fibroblast = "Fibroblast",
        Endothelial_Cell = "Endothelial Cell", Myeloid_Cell = "Myeloid Cell",
        Pericyte = "Pericyte", Glia_Cell = "Glia Cell")
## cardiomyocyte subtypes renamed as in the published figure (plotting only)
CM_AS_PUBLISHED <- c(Cardiomyocyte_2 = "Cardiomyocyte_1", Cardiomyocyte_3 = "Cardiomyocyte_2",
                     Cardiomyocyte_1 = "Cardiomyocyte_3")
TIMEPOINTS <- c(Sham = "Sham", day1 = "Day 1", day7 = "Day 7", day14 = "Day 14")
HIGHLIGHT <- c(cell_type = "Cardiomyocyte_2", module = "M3")
PT_UMAP <- 0.15; PT_NEB <- 0.12; FS <- 5
VIR <- viridis::viridis(11)

ref <- readRDS(REFERENCE_RDS)
modules <- read.delim(WGCNA_MODULES, colClasses = "character")
for (m in names(modules)) {
   g <- unique(modules[[m]]); g <- g[!is.na(g) & g != ""]
   ref <- AddModuleScore(ref, features = list(g), name = paste0(m, "_Score"), assay = "RNA")
}
SCORES <- paste0(names(modules), "_Score1")

ct <- ref$cell.type
ct[ct %in% names(CM_AS_PUBLISHED)] <- CM_AS_PUBLISHED[ct[ct %in% names(CM_AS_PUBLISHED)]]
tab <- data.frame(Embeddings(ref, "umap")[, 1:2], cell.type = factor(ct, levels = names(CT)),
                  timepoint = ref$timepoint,
                  t(as.matrix(GetAssayData(ref, assay = "RNA", layer = "data")[GENES, ])),
                  ref@meta.data[, SCORES], check.names = FALSE)
names(tab)[1:2] <- c("umap_1", "umap_2")
rm(ref); gc()
cat(sprintf("%d cells\n", nrow(tab))); print(table(tab$cell.type, tab$timepoint))

th <- theme_classic(base_size = FS) +
   theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black", size = FS - 0.5),
         axis.line = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2),
         axis.ticks.length = unit(0.6, "pt"), plot.title = element_text(size = FS + 1, margin = margin(b = 1)),
         plot.margin = margin(1, 2, 1, 1), legend.position = "bottom", legend.ticks = element_blank(),
         legend.title = element_text(size = FS - 0.5), legend.text = element_text(size = FS - 0.5),
         legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-2, 0, 0, 0))
cbar <- function(title, w = 30) guide_colourbar(title = title, title.position = "bottom", title.hjust = 0.5,
                                                barwidth = unit(w, "pt"), barheight = unit(3.5, "pt"))

## ---- UMAP ----------------------------------------------------------------------------------
xr <- range(tab$umap_1); yr <- range(tab$umap_2)
x0 <- xr[1] - 0.04 * diff(xr); y0 <- yr[1] - 0.04 * diff(yr)
axis_arrow <- function(xend, yend) annotate("segment", x = x0, y = y0, xend = xend, yend = yend, linewidth = 0.25,
                                            arrow = arrow(length = unit(2, "pt"), type = "closed"))
p_umap <- ggplot(tab, aes(umap_1, umap_2, colour = cell.type)) +
   geom_point(size = PT_UMAP, shape = 16, stroke = 0) +
   scale_colour_manual(values = setNames(viridis::viridis(length(CT)), names(CT)), labels = CT, name = NULL) +
   axis_arrow(x0 + 0.22 * diff(xr), y0) + axis_arrow(x0, y0 + 0.22 * diff(yr)) +
   annotate("text", x = x0, y = y0 - 0.02 * diff(yr), label = "UMAP1", hjust = 0, vjust = 1, size = (FS - 0.5) / .pt) +
   annotate("text", x = x0 - 0.02 * diff(xr), y = y0, label = "UMAP2", angle = 90, hjust = 0, vjust = 0, size = (FS - 0.5) / .pt) +
   coord_cartesian(clip = "off") +
   guides(colour = guide_legend(override.aes = list(size = 1.6), keyheight = unit(5.5, "pt"))) +
   theme_void(base_size = FS) +
   theme(legend.position = "right", legend.text = element_text(size = FS), legend.key.width = unit(4, "pt"),
         plot.margin = margin(2, 2, 6, 6))

## ---- module-score heatmap ------------------------------------------------------------------
agg <- aggregate(tab[, SCORES], by = list(cell.type = tab$cell.type), FUN = mean)
z <- scale(as.matrix(agg[match(names(CT), agg$cell.type), SCORES]))
rownames(z) <- names(CT)
write.table(data.frame(cell_type = rownames(z),
                       cell_type_in_reference = c(setNames(names(CM_AS_PUBLISHED), CM_AS_PUBLISHED),
                                                  setNames(names(CT)[-(1:3)], names(CT)[-(1:3)]))[rownames(z)],
                       round(z, 4), check.names = FALSE),
            file.path(OUT_DIR, "figS1E_module_score_zscore.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
hm <- expand.grid(ct = names(CT), mod = SCORES, stringsAsFactors = FALSE)
hm$z <- z[cbind(hm$ct, hm$mod)]
hl <- hm[hm$ct == HIGHLIGHT[["cell_type"]] & hm$mod == paste0(HIGHLIGHT[["module"]], "_Score1"), ]
hm$ct <- factor(hm$ct, levels = rev(names(CT))); hl$ct <- factor(hl$ct, levels = levels(hm$ct))
hm$mod <- factor(hm$mod, levels = SCORES, labels = paste("Module", seq_along(SCORES)))
hl$mod <- factor(hl$mod, levels = SCORES, labels = paste("Module", seq_along(SCORES)))
zmax <- max(abs(hm$z))
p_heat <- ggplot(hm, aes(mod, ct, fill = z)) +
   geom_tile(colour = "white", linewidth = 0.3) +
   geom_tile(data = hl, fill = NA, colour = "black", linewidth = 0.6) +
   scale_fill_gradientn(colours = colorRampPalette(c(viridis::turbo(11)[2], "white", viridis::turbo(11)[9]))(256),
                        limits = c(-zmax, zmax), breaks = c(-2, 0, 2), guide = cbar("Module score (z-score)", w = 40)) +
   scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
   coord_fixed() + theme_void(base_size = FS) +
   theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = FS - 0.5),
         legend.position = "bottom", legend.ticks = element_blank(),
         legend.title = element_text(size = FS - 0.5), legend.text = element_text(size = FS - 0.5),
         legend.margin = margin(1, 0, 0, 0), plot.margin = margin(1, 1, 1, 2))

## ---- Nebulosa ------------------------------------------------------------------------------
dens <- function(d) {
   emb <- as.matrix(d[, c("umap_1", "umap_2")])
   k <- sapply(GENES, function(g) Nebulosa:::calculate_density(d[[g]], emb, "ks", 1))
   cbind(k, joint = apply(k, 1, prod))
}
neb <- function(d, value, lim, breaks, labels, title, legend, y_axis = TRUE, w = 30) {
   d$value <- value
   p <- ggplot(d, aes(umap_1, umap_2, colour = value)) +
      geom_point(size = PT_NEB, shape = 16, stroke = 0) +
      scale_colour_gradientn(colours = VIR, limits = lim, oob = scales::squish, breaks = breaks,
                             labels = labels, guide = cbar(legend, w)) +
      coord_cartesian(xlim = xr, ylim = yr) + labs(x = "umap 1", y = "umap 2", title = title) + th
   if (!y_axis) p <- p + labs(y = NULL) + theme(axis.text.y = element_blank())
   p
}
fix2 <- function(x) sprintf("%.2f", x)
sci <- scales::label_scientific()

all_d <- dens(tab)
row2 <- wrap_plots(neb(tab, all_d[, 1], LIM$rtn4, c(0, 0.01), fix2, "Rtn4", "Density"),
                   neb(tab, all_d[, 2], LIM$chheaf, c(0, 0.04), fix2, "Chheaf", "Density"),
                   neb(tab, all_d[, 3], LIM$joint, c(0, 2e-4), sci, "Rtn4 + Chheaf", "Joint density"),
                   plot_spacer(), nrow = 1, widths = c(1, 1, 1, 0.22))
by_tp <- lapply(names(TIMEPOINTS), function(tp) {
   d <- tab[tab$timepoint == tp, ]
   neb(d, dens(d)[, "joint"], LIM$joint_tp, c(0, 6e-4), sci, TIMEPOINTS[[tp]], NULL,
       y_axis = tp == "Sham", w = 40)
})
p_hdr <- ggplot() + annotate("text", x = 0.5, y = 0.62, label = "Post-MI", size = (FS + 1) / .pt) +
   annotate("segment", x = 0.02, xend = 0.98, y = 0.12, yend = 0.12, linewidth = 0.3) +
   scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
   coord_cartesian(clip = "off") + theme_void()
row3 <- p_hdr + by_tp[[1]] + by_tp[[2]] + by_tp[[3]] + by_tp[[4]] +
   plot_layout(design = "#BBB\nCDEF", heights = c(0.1, 1), widths = c(1.12, 1, 1, 1), guides = "collect") &
   theme(legend.position = "bottom")

panel <- wrap_elements(full = (p_umap | p_heat) + plot_layout(widths = c(1, 1.15))) /
   wrap_elements(full = row2) / wrap_elements(full = row3) +
   plot_layout(heights = c(1, 1.05, 1.1)) +
   plot_annotation(tag_levels = list(c("E", "", "")), theme = theme(plot.margin = margin(2, 2, 2, 2))) &
   theme(plot.tag = element_text(size = FS + 4))

out <- file.path(OUT_DIR, "figS1E_snrnaseq_chheaf_rtn4")
ggsave(paste0(out, ".pdf"), panel, width = PANEL_IN[1], height = PANEL_IN[2], device = cairo_pdf)
ggsave(paste0(out, ".png"), panel, width = PANEL_IN[1], height = PANEL_IN[2], dpi = 600)
cat(sprintf("outlined cell %s x %s: z = %.2f; written %s.{pdf,png}\n",
            HIGHLIGHT[["module"]], HIGHLIGHT[["cell_type"]], hl$z, out))
print(round(z[, , drop = FALSE], 2))
