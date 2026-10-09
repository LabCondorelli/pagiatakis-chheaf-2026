#!/usr/bin/env Rscript
## Fig. S1F - co-expression of each lincRNA of Fig. S1C with its flanking gene, one row per pair.
## Same score and scale as Fig. 1D; lincRNAs in < 20 Visium spots are dimmed.
## Usage: Rscript figS1F_spatial_lincrna_pairs.R
## Output: results/figS1F_spatial_lincrna_pairs.{pdf,png}, results/figS1F_coexpression_per_pair.tsv

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork); library(grid) })
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "00_config.R"))
source(file.path(here, "utils", "spatial_helpers.R"))
options(future.globals.maxSize = 1e16)

PAGE_IN <- c(8.5, 11)
MARGIN_IN <- 0.15; GUTTER_IN <- 0.15; LEGEND_IN <- 0.55
N_BLOCK_COLUMNS <- 2; LABEL_FRAC <- 0.24; BLOCK_GAP <- 0.9; HEADER_H <- 0.30
NOT_EVALUABLE <- 20
NOT_DRAWN <- "Gm12092 Rtn4"
PT <- c(GSM5355663_WT_Sham = 0.212, GSM5943190_WT_MI_day1_2 = 0.222,
        GSM5943192_WT_MI_day7_2 = 0.254, GSM5943195_WT_MI_day14_3 = 0.242)
COLUMNS <- c("Sham", "Day 1", "Day 7", "Day 14")
PAL <- viridis::viridis(11)

## Fig. S1C order; Clhc1 and Rps18-ps6 are not expressed.
S1C <- data.frame(
   linc  = c("Gm12092", "Gm13481", "Gm13481", "Gm5086", "Gm30085", "Gm30085", "A730036I17Rik",
             "A730036I17Rik", "Gm17501", "Gm17501", "Gm42047", "Gm42047"),
   label = c("Chheaf (Gm12092)", "Gm13481", "Gm13481", "Gm5086", "Gm30085", "Gm30085", "A730036I17Rik",
             "A730036I17Rik", "Gm17501", "Gm17501", "Gm42047", "Gm42047"),
   gene  = c("Rtn4", "Zeb2", "Acvr2a", "Enc1", "Pdzrn4", "Gxylt1", "Slc20a1", "Ckap2l",
             "Ccn1", "Ddah1", "Klhdc4", "Slc7a5"),
   side  = c("downstream", "upstream", "downstream", "upstream", "upstream", "downstream",
             "upstream", "downstream", "upstream", "downstream", "upstream", "downstream"),
   stringsAsFactors = FALSE)
stopifnot(setequal(paste(S1C$linc, S1C$gene), paste(PAIRS$linc, PAIRS$gene)))
S1C$key <- paste(S1C$linc, S1C$gene)

vis <- subset(readRDS(VISIUM_RDS), orig.ident %in% PUB_SLIDES)
vis <- NormalizeData(vis, assay = "Spatial", verbose = FALSE)
DefaultAssay(vis) <- "Spatial"
for (sl in PUB_SLIDES) vis <- rotate_section(vis, sl, SLIDE_ROTATION[[sl]])
proj <- read.csv(TANGRAM_GENES, row.names = 1, check.names = FALSE)
if (!all(colnames(vis) %in% rownames(proj))) stop("Tangram spot ids do not match the Visium object")
proj <- proj[colnames(vis), ]
expr <- GetAssayData(vis, assay = "Spatial", layer = "data")
cnt  <- GetAssayData(vis, assay = "Spatial", layer = "counts")

genes <- unique(c(S1C$linc, S1C$gene))
nspot <- setNames(vapply(genes, function(g) sum(cnt[g, ] > 0), numeric(1)), genes)
pct   <- 100 * nspot / ncol(vis)
note  <- function(g, short = FALSE) sprintf(if (short) "%d (%.1f%%)" else "%d spots (%.1f%%)", nspot[[g]], pct[[g]])

COLOC <- setNames(lapply(seq_len(nrow(S1C)), function(i)
   coexpression_score(proj[[tolower(S1C$linc[i])]], expr[S1C$gene[i], ])), S1C$key)
LIMS <- c(min(unlist(COLOC)), COEXPR_CAP)

q <- function(x) c(min = min(x), median = median(x), q95 = unname(quantile(x, 0.95)),
                   q99 = unname(quantile(x, 0.99)), max = max(x), pct_above_cap = 100 * mean(x > COEXPR_CAP))
tab <- cbind(S1C[, c("label", "gene", "side")], drawn = S1C$key != NOT_DRAWN,
             linc_spots = nspot[S1C$linc], gene_spots = nspot[S1C$gene],
             t(vapply(COLOC, q, numeric(6))))
tab[] <- lapply(tab, function(x) if (is.double(x)) signif(x, 4) else x)
write.table(tab, file.path(OUT_DIR, "figS1F_coexpression_per_pair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab, row.names = FALSE)

## ---- section crop: largest connected group of spots, plus one lattice pitch ------------------
crop_limits <- function(p, k = 8L) {
   ld <- ggplot2::layer_data(p, 1)[, c("x", "y")]
   nn <- RANN::nn2(ld, k = k + 1L)
   pitch <- median(nn$nn.dists[, 2])
   keep <- as.vector(nn$nn.dists[, -1]) <= 1.6 * pitch
   g <- igraph::add_edges(igraph::make_empty_graph(nrow(ld), directed = FALSE),
                          rbind(rep(seq_len(nrow(ld)), k)[keep], as.vector(nn$nn.idx[, -1])[keep]))
   cm <- igraph::components(g)
   sel <- cm$membership == which.max(cm$csize)
   list(x = range(ld$x[sel]) + c(-1, 1) * pitch, y = range(ld$y[sel]) + c(-1, 1) * pitch)
}
vis$coloc <- COLOC[[1]]
CROP <- lapply(setNames(PUB_SLIDES, PUB_SLIDES), function(sl)
   crop_limits(SpatialFeaturePlot(vis, features = "coloc", images = sl, pt.size.factor = PT[[sl]], image.alpha = 0)))
ASPECT <- min(1.45, max(0.85, max(vapply(CROP, function(cl) diff(cl$y) / diff(cl$x), numeric(1)))))

map_row <- function(key) {
   vis$coloc <- COLOC[[key]]
   wrap_plots(lapply(PUB_SLIDES, function(sl)
      SpatialFeaturePlot(vis, features = "coloc", images = sl, pt.size.factor = PT[[sl]], image.alpha = 0) +
         scale_fill_gradientn(colours = PAL, limits = LIMS, oob = scales::squish) +
         coord_fixed(xlim = CROP[[sl]]$x, ylim = CROP[[sl]]$y, expand = FALSE) +
         theme(legend.position = "none", plot.margin = margin(1, 1, 1, 1))), nrow = 1)
}

## ---- page geometry: blocks of one lincRNA each, over two columns, never split ---------------
drawn <- S1C[S1C$key != NOT_DRAWN, ]
blocks <- split(drawn, factor(drawn$linc, levels = unique(drawn$linc)))
n_rows <- vapply(blocks, nrow, numeric(1))
col_w <- (PAGE_IN[1] - 2 * MARGIN_IN - GUTTER_IN * (N_BLOCK_COLUMNS - 1)) / N_BLOCK_COLUMNS
map_w <- col_w * (1 - LABEL_FRAC) / length(PUB_SLIDES)
row_h <- map_w * ASPECT
fs <- max(0.62, min(1, col_w * LABEL_FRAC / 1.25))

quota <- ceiling(sum(n_rows) / N_BLOCK_COLUMNS)
cols <- list(integer(0)); used <- 0
for (i in seq_along(blocks)) {
   if (used + n_rows[i] > quota && length(cols) < N_BLOCK_COLUMNS) { cols <- c(cols, list(integer(0))); used <- 0 }
   cols[[length(cols)]] <- c(cols[[length(cols)]], i); used <- used + n_rows[i]
}
content_u <- HEADER_H + max(vapply(cols, function(k) sum(n_rows[k]) + BLOCK_GAP * (length(k) - 1), numeric(1)))
page_u <- (PAGE_IN[2] - 2 * MARGIN_IN) / row_h
if (content_u + LEGEND_IN / row_h > page_u) stop("the blocks do not fit on one page")

block_plot <- function(b, letter = NULL) {
   linc <- b$linc[1]
   lab <- s1c_label_grob(b$label[1], note(linc), data.frame(sym = b$gene, side = b$side,
                         note = vapply(b$gene, note, character(1), short = TRUE)),
                         greyed = nspot[[linc]] < NOT_EVALUABLE, letter = letter, fs = fs)
   wrap_plots(list(wrap_elements(full = lab), wrap_plots(lapply(b$key, map_row), ncol = 1)),
              ncol = 2, widths = c(LABEL_FRAC, 1 - LABEL_FRAC))
}
build_column <- function(idx, first) {
   parts <- list(wrap_plots(list(plot_spacer(), wrap_elements(full = column_header_grob(COLUMNS, fs = fs))),
                            ncol = 2, widths = c(LABEL_FRAC, 1 - LABEL_FRAC)))
   h <- HEADER_H
   for (j in seq_along(idx)) {
      parts <- c(parts, list(block_plot(blocks[[idx[j]]], letter = if (first && j == 1) "A")))
      h <- c(h, n_rows[idx[j]])
      if (j < length(idx)) { parts <- c(parts, list(plot_spacer())); h <- c(h, BLOCK_GAP) }
   }
   if (sum(h) < content_u) { parts <- c(parts, list(plot_spacer())); h <- c(h, content_u - sum(h)) }
   wrap_plots(parts, ncol = 1, heights = h)
}

cbar <- colourbar_grob(PAL, "lincRNA / flanking gene co-expression",
   note = sprintf("cap %s = both genes in the top 5%% of their spots; the %.2f%% of spots above the cap are compressed, not clipped",
                  format(COEXPR_CAP), 100 * mean(unlist(COLOC) > COEXPR_CAP)))
pad <- max(0, (page_u - content_u - LEGEND_IN / row_h) / 2)
fig <- wrap_plots(list(plot_spacer(),
                       wrap_plots(lapply(seq_along(cols), function(k) build_column(cols[[k]], k == 1)), nrow = 1),
                       wrap_elements(full = cbar), plot_spacer()),
                  ncol = 1, heights = c(pad, content_u, LEGEND_IN / row_h, pad)) +
   plot_annotation(theme = theme(plot.margin = margin(MARGIN_IN, MARGIN_IN, MARGIN_IN, MARGIN_IN, unit = "in")))

out <- file.path(OUT_DIR, "figS1F_spatial_lincrna_pairs")
ggsave(paste0(out, ".pdf"), fig, width = PAGE_IN[1], height = PAGE_IN[2])
ggsave(paste0(out, ".png"), fig, width = PAGE_IN[1], height = PAGE_IN[2], dpi = 300)
cat(sprintf("scale %.3f - %.1f over 12 pairs; %d rows drawn; written %s.{pdf,png}\n",
            LIMS[1], COEXPR_CAP, nrow(drawn), out))
