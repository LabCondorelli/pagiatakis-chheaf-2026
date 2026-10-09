#!/usr/bin/env Rscript
## Exports the two plain tables that Fig. 1A, Fig. 1B and Data S2 read, from the output of
## the lab edgeR/voom pipeline on the mouse TAC time course (GENCODE vM25, 8 groups of 2:
## sham and TAC at 1, 2, 4 and 7 days).
##
## The matrix is the voom log2 CPM, log2((counts + 0.5) / (TMM library size + 1) * 1e6),
## on the 15,565 genes that pass filterByExpr.
##
## Usage: Rscript 00_export_mouse_tac.R
## Output: results/mouse_tac_log2cpm.tsv, results/mouse_tac_DE_calls.tsv

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

## log2 CPM, one column per sample
m <- read.delim(MOUSE_TAC_CPM, check.names = FALSE, stringsAsFactors = FALSE)
write.table(m, file.path(OUT_DIR, "mouse_tac_log2cpm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Differential calls already thresholded by the pipeline at FDR < 0.05 and |log2FC| > 0.5,
## coded -1 / 0 / +1 in the Results.* columns, one per time point (TAC vs sham).
de <- read.delim(MOUSE_TAC_DE_SIGN, check.names = FALSE)
write.table(de, file.path(OUT_DIR, "mouse_tac_DE_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("log2 CPM: %d genes x %d samples\nDE calls: %d genes\n",
            nrow(m), sum(grepl("^(Sham|TAC)[.]", colnames(m))), nrow(de)))
