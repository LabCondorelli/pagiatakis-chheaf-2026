#!/usr/bin/env Rscript
## Fig. 2D - TPM of Chheaf and Rtn4 in WT and Chheaf-KO hearts, sham and 1, 2, 4 weeks after TAC.
## Fig. S2B - Spearman correlation between the two, in TAC hearts at 1 week.
##
## Both use the TPM matrix of the gencode tree, the one that contains the lincRNA.
##
## Usage: Rscript fig2D_S2B_ko_timecourse.R
## Output: results/fig2D_timecourse.{pdf,tsv}, results/figS2B_spearman.{pdf,tsv}

suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

GENES <- c(Chheaf = "ENSMUSG00000084850", Rtn4 = "ENSMUSG00000020458")
pal   <- c(wt = "#377EB8", ko = "#E41A1C")

## The matrix is stored as log2(TPM + 0.5); it is put back on the TPM scale here.
tpm <- read.delim(KO_TPM, check.names = FALSE, stringsAsFactors = FALSE)
smp <- setdiff(colnames(tpm), c("gene_id", "symbol"))
rownames(tpm) <- sub("[.].*$", "", tpm$gene_id)
tpm[, smp] <- 2^as.matrix(tpm[, smp]) - 0.5

## Sample names are Sham.WT.1W_repA and similar; sham exists only at 1 week.
meta <- data.frame(sample = smp,
                   treatment = tolower(sub("^([A-Za-z]+)[.].*$", "\\1", smp)),
                   genotype  = tolower(sub("^[A-Za-z]+[.]([A-Za-z]+)[.].*$", "\\1", smp)),
                   week      = as.numeric(sub("^.*[.]([0-9]+)W_.*$", "\\1", smp)))
meta$x <- ifelse(meta$treatment == "sham", 0, meta$week)   # sham at position 0

long <- do.call(rbind, lapply(names(GENES), function(g)
  data.frame(gene = g, meta, tpm = as.numeric(tpm[GENES[g], smp]))))

## ---- Fig. 2D ----------------------------------------------------------------------------
pl <- lapply(names(GENES), function(g)
  ggplot(subset(long, gene == g), aes(x, tpm, colour = genotype, fill = genotype)) +
    geom_point(size = 2) + geom_smooth(se = TRUE) +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
    scale_x_continuous(breaks = c(0, 1, 2, 4), labels = c("Sham", "1wk", "2wks", "4wks")) +
    labs(title = g, x = NULL, y = "TPM") + theme_classic())

pdf(file.path(OUT_DIR, "fig2D_timecourse.pdf"), width = 8, height = 3)
print(ggarrange(plotlist = pl, nrow = 1, ncol = 2, common.legend = TRUE))
invisible(dev.off())
write.table(long, file.path(OUT_DIR, "fig2D_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- Fig. S2B ---------------------------------------------------------------------------
w <- subset(meta, treatment == "tac" & week == 1)
d <- data.frame(sample = w$sample, genotype = w$genotype,
                Rtn4   = as.numeric(tpm[GENES["Rtn4"],   w$sample]),
                Chheaf = as.numeric(tpm[GENES["Chheaf"], w$sample]))
ct <- suppressWarnings(cor.test(d$Rtn4, d$Chheaf, method = "spearman"))
cat(sprintf("Fig. S2B: n = %d, rho = %.3f, P = %.4g\n", nrow(d), ct$estimate, ct$p.value))

pdf(file.path(OUT_DIR, "figS2B_spearman.pdf"), width = 4, height = 3.5)
print(ggplot(d, aes(Rtn4, Chheaf)) + geom_smooth(method = "lm", formula = y ~ x, colour = "grey40") +
        geom_point(aes(colour = genotype), size = 3) + scale_colour_manual(values = pal) +
        stat_cor(method = "spearman") + labs(x = "Rtn4 (TPM)", y = "Chheaf (TPM)") + theme_classic())
invisible(dev.off())
write.table(d, file.path(OUT_DIR, "figS2B_spearman.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
