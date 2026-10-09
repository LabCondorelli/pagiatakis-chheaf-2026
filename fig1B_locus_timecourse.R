#!/usr/bin/env Rscript
## Fig. 1B - log2 CPM of Chheaf (Gm12092) and of its two flanking coding genes, Clhc1 and Rtn4,
## over the TAC time course, sham against TAC.
##
## Usage: Rscript fig1B_locus_timecourse.R   (run 00_export_mouse_tac.R first)
## Output: results/fig1B_locus_timecourse.pdf, results/fig1B_locus_timecourse.tsv

suppressPackageStartupMessages({ library(ggplot2); library(ggpubr); library(RColorBrewer) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

cpm <- read.delim(file.path(OUT_DIR, "mouse_tac_log2cpm.tsv"), check.names = FALSE)
smp <- grep("^(Sham|TAC)[.]", colnames(cpm), value = TRUE)

## Chheaf and its expressed flanking coding gene. The other flanking gene, Clhc1, does not pass
## the expression filter of the pipeline and is absent from the matrix.
GENES <- c("Gm12092", "Rtn4")
pal   <- brewer.pal(9, "Set1")[2:1]   # blue = sham, red = TAC

long <- do.call(rbind, lapply(GENES, function(gn) {
  i <- which(cpm$symbol == gn)[1]
  stopifnot(!is.na(i))
  data.frame(gene  = factor(gn, levels = GENES),
             cond  = sub("[.].*$", "", smp),
             time  = sub("_rep[AB]$", "", sub("^.*[.]", "", smp)),
             value = as.numeric(cpm[i, smp]))
}))

## Cubic spline with fixed knots, as in the original analysis
pl <- lapply(GENES, function(gn)
  ggplot(subset(long, gene == gn), aes(time, value, group = cond, colour = cond, fill = cond)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", fx = TRUE, k = 4)) +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                       breaks = scales::pretty_breaks(n = 5)) +
    labs(title = gn, x = NULL, y = "log2 CPM") +
    theme_classic() + theme(legend.position = "none"))

pdf(file.path(OUT_DIR, "fig1B_locus_timecourse.pdf"), width = 6, height = 3)
print(ggarrange(plotlist = pl, nrow = 1, ncol = 2))
invisible(dev.off())

write.table(long, file.path(OUT_DIR, "fig1B_locus_timecourse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d genes x %d samples written\n", length(GENES), length(smp)))
