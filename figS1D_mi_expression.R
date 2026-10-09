#!/usr/bin/env Rscript
## Fig. S1D - Chheaf and Rtn4 in sham and post-myocardial infarction hearts (4 versus 4).
##
## Expression is the log2 CPM; the statistics come from the limma-voom analysis, the model the
## Methods describe. These are the two FDR printed in the panel, 0.0210 for Chheaf and
## 1.4e-05 for Rtn4.
##
## Usage: Rscript figS1D_mi_expression.R
## Output: results/figS1D_mi_expression.{pdf,tsv}, results/figS1D_mi_stats.tsv

suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

GENES <- c(Chheaf = "ENSMUSG00000084850", Rtn4 = "ENSMUSG00000020458")

cpm <- read.delim(MI_CPM, check.names = FALSE, stringsAsFactors = FALSE)
smp <- setdiff(colnames(cpm), c("ensgID", "symbol"))
long <- do.call(rbind, lapply(names(GENES), function(g) {
  i <- which(cpm$ensgID == GENES[g])[1]
  data.frame(gene = factor(g, levels = names(GENES)),
             sample = smp,
             group = factor(sub("_rep[A-D]$", "", smp), levels = c("sham", "infarct")),
             log2cpm = as.numeric(cpm[i, smp]))
}))

d <- read.delim(MI_DE, check.names = FALSE, stringsAsFactors = FALSE)
stats <- do.call(rbind, lapply(names(GENES), function(g) {
  i <- which(sub("[.].*$", "", d$Genes.ensembl_gene_id) == GENES[g])[1]
  ## one contrast only, so the columns carry no suffix in our tables and the contrast name
  ## in those of 00_de_and_cpm.R
  co <- grep("^Coef", names(d), value = TRUE)[1]
  fd <- grep("^P[.]value[.]adj", names(d), value = TRUE)[1]
  data.frame(gene = g, log2FC = d[[co]][i], FDR = d[[fd]][i])
}))
print(stats, row.names = FALSE, digits = 3)

pdf(file.path(OUT_DIR, "figS1D_mi_expression.pdf"), width = 5, height = 3)
print(ggplot(long, aes(group, log2cpm, fill = group)) +
        stat_summary(fun = mean, geom = "col", width = 0.6) +
        geom_point(size = 2) + facet_wrap(~ gene, scales = "free_y") +
        scale_fill_manual(values = c(sham = "#377EB8", infarct = "#E41A1C")) +
        labs(x = NULL, y = "log2 CPM") + theme_classic() + theme(legend.position = "none"))
invisible(dev.off())

write.table(long,  file.path(OUT_DIR, "figS1D_mi_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(stats, file.path(OUT_DIR, "figS1D_mi_stats.tsv"),      sep = "\t", quote = FALSE, row.names = FALSE)
