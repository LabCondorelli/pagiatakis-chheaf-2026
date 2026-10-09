#!/usr/bin/env Rscript
## Fig. S2E - GSEA of the DisGeNET "Hypertrophic Cardiomyopathy" gene set on TAC KO versus
## TAC WT hearts, ranked by the moderated t statistic.
##
## Three things are pinned so that a rerun gives the same numbers: the gene set library is read
## from the frozen copy in data/, the orthologue mapping from the frozen table in data/, and the
## seed and the number of processes are fixed. The mapping and the ranking are those of the run
## that produced the panel.
##
## Usage: Rscript figS2E_gsea_disgenet.R
## Output: results/figS2E_gsea.tsv, results/figS2E_leading_edge.tsv, results/figS2E_enrichment.pdf

suppressPackageStartupMessages({ library(fgsea); library(ggplot2) })
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(HERE, "00_config.R"))

GMT     <- p("CHH_GMT",     file.path(HERE, "data", "disgenet_harmonizome.gmt.gz"))
ORTHO   <- p("CHH_ORTHO",   file.path(HERE, "data", "mouse_human_orthologues.tsv"))
PATHWAY <- "Hypertrophic Cardiomyopathy"
CONTRASTS <- c("1wk", "2wk", "4wk")
NPERM   <- 10000
LEADING_EDGE_CONTRAST <- p("CHH_LE_CONTRAST", "2wk")   # the contrast shown in the panel

de <- read.delim(KO_DE_ALL, check.names = FALSE, stringsAsFactors = FALSE)

## Mouse to human orthologues, from the Ensembl BioMart table frozen in data/. The join is on
## the Ensembl gene id, not on the symbol, and where a mouse gene has several human homologues
## the first one is kept: both as in the analysis that produced the panel. Rebuilding the table
## against a newer Ensembl release changes the enrichment, see data/README.md.
bm <- read.delim(ORTHO, stringsAsFactors = FALSE)
de$human <- bm$human_symbol[match(sub("[.].*$", "", de$Genes.gene_id), bm$ensembl_gene_id)]

if (!file.exists(GMT))
  stop("gene set library not found. Harmonizome is CC BY-NC-SA 4.0 and is not redistributed\n",
       "  here; get it with:  bash data/fetch_disgenet.sh")
pathways <- gmtPathways(GMT)

rank_by <- function(col) {
  s <- de[[col]]; g <- de$human
  ok <- !is.na(s) & !is.na(g) & g != ""
  s <- s[ok]; names(s) <- g[ok]
  sort(s[!duplicated(names(s))], decreasing = TRUE)
}

set.seed(123)
res <- do.call(rbind, lapply(CONTRASTS, function(w) {
  st <- rank_by(paste0("t.tac_ko_vs_tac_wt_", w))
  r  <- fgsea(pathways, st, nPermSimple = NPERM, nproc = 1,
              minSize = 1, maxSize = length(st) - 1)
  r  <- r[r$pathway == PATHWAY, ]
  data.frame(contrast = w, NES = r$NES, pval = r$pval, padj = r$padj,
             size = r$size, leading_edge = lengths(r$leadingEdge))
}))
print(res, row.names = FALSE, digits = 3)
write.table(res, file.path(OUT_DIR, "figS2E_gsea.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Leading edge and enrichment curve of the contrast shown in the panel
set.seed(123)
st <- rank_by(paste0("t.tac_ko_vs_tac_wt_", LEADING_EDGE_CONTRAST))
r  <- fgsea(pathways, st, nPermSimple = NPERM, nproc = 1, minSize = 1, maxSize = length(st) - 1)
le <- r$leadingEdge[[which(r$pathway == PATHWAY)]]
write.table(data.frame(gene = le), file.path(OUT_DIR, "figS2E_leading_edge.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("leading edge at %s: %d genes\n", LEADING_EDGE_CONTRAST, length(le)))

pdf(file.path(OUT_DIR, "figS2E_enrichment.pdf"), width = 5, height = 3.5)
print(plotEnrichment(pathways[[PATHWAY]], st) +
        labs(title = sprintf("%s, TAC KO vs TAC WT %s", PATHWAY, LEADING_EDGE_CONTRAST)))
invisible(dev.off())
