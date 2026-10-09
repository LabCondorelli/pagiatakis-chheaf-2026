#!/usr/bin/env Rscript
## Data S1 - differential expression of TAC KO versus TAC WT at 1, 2 and 4 weeks.
##
## One table per time point, from the protein-coding tree of the Chheaf-KO RNA-seq.
##
## Usage: Rscript dataS1_ko_vs_wt_de.R
## Output: results/dataS1_tac_ko_vs_tac_wt_{1,2,4}wk.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

WEEKS <- c("1wk", "2wk", "4wk")
de <- read.delim(KO_DE_ALL, check.names = FALSE, stringsAsFactors = FALSE)

for (w in WEEKS) {
  cc <- paste0("tac_ko_vs_tac_wt_", w)
  d <- data.frame(GeneID   = de$Genes.gene_id,
                  GeneName = de$Genes.gene_name,
                  logFC    = de[[paste0("Coef.", cc)]],
                  t        = de[[paste0("t.", cc)]],
                  P.Value  = de[[paste0("P.value.", cc)]],
                  adj.P.Val = de[[paste0("P.value.adj.", cc)]],
                  Significant = ifelse(de[[paste0("Results.", cc)]] == 0, "no", "yes"),
                  Biotype  = de$Genes.gene_biotype)
  d <- d[order(d$adj.P.Val, d$P.Value), ]
  f <- file.path(OUT_DIR, sprintf("dataS1_tac_ko_vs_tac_wt_%s.tsv", w))
  write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  %-4s %d genes, %d significant -> %s\n", w, nrow(d),
              sum(d$Significant == "yes"), basename(f)))
}
