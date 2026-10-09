#!/usr/bin/env Rscript
## Fig. S2C - differential expression summary of the genes uniquely down-regulated at 1 week
## after TAC in KO versus WT hearts.
##
## "Uniquely down-regulated" means that at 1 week the TAC KO versus TAC WT contrast calls genes
## down and none up: three genes, Rtn4, Serpinb1a and Rbms3. The script derives them from the
## data, then reports each one across all six contrasts, with the log2FC where the call is
## significant (FDR < 0.05 and |log2FC| > 0.5) and NS elsewhere.
##
## Usage: Rscript figS2C_de_summary.R
## Output: results/figS2C_de_summary.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

WEEKS <- c("1wk", "2wk", "4wk")
CONTRASTS <- c(paste0("tac_wt_vs_sham_wt_", WEEKS),
               paste0("tac_ko_vs_tac_wt_", WEEKS), "tac_ko_vs_tac_wt")

pull <- function(tab, genes, tree) {
  d <- read.delim(tab, check.names = FALSE, stringsAsFactors = FALSE)
  do.call(rbind, lapply(genes, function(g) {
    i <- which(d$Genes.gene_name == g)[1]
    data.frame(gene = g, tree = tree, contrast = CONTRASTS,
               log2FC = sapply(CONTRASTS, function(cc) d[[paste0("Coef.", cc)]][i]),
               FDR    = sapply(CONTRASTS, function(cc) d[[paste0("P.value.adj.", cc)]][i]),
               call   = sapply(CONTRASTS, function(cc) d[[paste0("Results.", cc)]][i]))
  }))
}

## Genes that define the panel: down and nothing up, at 1 week, in the protein-coding tree
de  <- read.delim(KO_DE_ALL, check.names = FALSE, stringsAsFactors = FALSE)
sel <- de[["Results.tac_ko_vs_tac_wt_1wk"]]
cat(sprintf("TAC KO vs TAC WT at 1 week: %d down, %d up\n", sum(sel == -1), sum(sel == 1)))
genes <- sort(de$Genes.gene_name[sel == -1])
cat("  ", paste(genes, collapse = ", "), "\n")

out <- pull(KO_DE_ALL, genes, "pc")
## Chheaf itself is a lincRNA and lives in the separate lncRNA tree
if (file.exists(KO_DE_LNC)) out <- rbind(out, pull(KO_DE_LNC, "Gm12092", "lncRNA"))
out$shown <- ifelse(out$call == 0, "NS", sprintf("%.2f", out$log2FC))

print(out[, c("gene", "tree", "contrast", "log2FC", "FDR", "shown")], row.names = FALSE, digits = 3)
write.table(out, file.path(OUT_DIR, "figS2C_de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
