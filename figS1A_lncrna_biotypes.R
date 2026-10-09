#!/usr/bin/env Rscript
## Fig. S1A - GENCODE biotype composition of the lncRNA modulated in the TAC time course.
##
## Modulated means called significant (FDR < 0.05 and |log2FC| > 0.5) in at least one TAC vs
## sham comparison. The biotype names are those of GENCODE vM25, which keeps lincRNA and
## antisense separate; vM26 and later merge them into lncRNA.
##
## Usage: Rscript figS1A_lncrna_biotypes.R   (run 00_export_mouse_tac.R first)
## Output: results/figS1A_lncrna_biotypes.{pdf,tsv}

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

BIOTYPES <- c("lincRNA", "antisense", "sense_intronic",
              "bidirectional_promoter_lncRNA", "sense_overlapping")

de   <- read.delim(file.path(OUT_DIR, "mouse_tac_DE_calls.tsv"), check.names = FALSE)
RES  <- paste0("Results.tac_", c("1d", "2d", "4d", "7d"), "_vs_sham_", c("1d", "2d", "4d", "7d"))
keep <- de$Genes.gene_biotype %in% BIOTYPES & rowSums(as.matrix(de[, RES]) != 0) > 0

tab <- table(factor(de$Genes.gene_biotype[keep], levels = BIOTYPES))
tab <- tab[tab > 0]
print(tab)
cat(sprintf("total modulated lncRNA: %d\n", sum(tab)))

pdf(file.path(OUT_DIR, "figS1A_lncrna_biotypes.pdf"), width = 5, height = 5)
pie(as.numeric(tab), labels = sprintf("%s (%d)", names(tab), as.numeric(tab)),
    col = grey.colors(length(tab)), main = "Modulated lncRNA by biotype")
invisible(dev.off())

write.table(data.frame(biotype = names(tab), n = as.numeric(tab)),
            file.path(OUT_DIR, "figS1A_lncrna_biotypes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
