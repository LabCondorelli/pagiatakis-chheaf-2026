#!/usr/bin/env Rscript
## Fig. S3B, expression values - Chheaf and Rtn4 across mouse cardiomyocyte development
## (E14.5, neonatal, adult) and after TAC, from the developmental RNA-seq (GENCODE vM25,
## 2 replicates per stage).
##
## Reports three quantities, which differ for a transcript this low: raw CPM on total assigned
## reads, linear CPM after filtering and TMM, and the log2 CPM used by the pipeline.
##
## Usage: Rscript figS3B_development_cpm.R
## Output: results/figS3B_development_cpm.tsv

suppressPackageStartupMessages({ library(edgeR) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

x <- read.delim(DEV_COUNTS, skip = 1, check.names = FALSE)
cols <- grep("[.]bam$", colnames(x), value = TRUE)
m <- as.matrix(x[, cols])
rownames(m) <- sub("[.].*$", "", x[[1]])
colnames(m) <- basename(dirname(cols))          # one directory per sample

stage <- sub("_rep[AB]$", "", colnames(m))
d <- DGEList(counts = m, group = stage)
d <- calcNormFactors(d[filterByExpr(d), , keep.lib.sizes = FALSE])

GENES <- c(Chheaf = "ENSMUSG00000084850", Rtn4 = "ENSMUSG00000020458")
raw <- t(t(m) / colSums(m)) * 1e6                # before filtering and TMM
lin <- cpm(d); lg <- cpm(d, log = TRUE)

out <- do.call(rbind, lapply(names(GENES), function(g) {
  i <- which(rownames(d) == GENES[g]); j <- which(rownames(m) == GENES[g])
  data.frame(gene = g, stage = unique(stage),
             reads    = as.numeric(tapply(m[j, ],   stage, sum)[unique(stage)]),
             cpm_raw  = as.numeric(tapply(raw[j, ], stage, mean)[unique(stage)]),
             cpm_tmm  = as.numeric(tapply(lin[i, ], stage, mean)[unique(stage)]),
             log2_cpm = as.numeric(tapply(lg[i, ],  stage, mean)[unique(stage)]))
}))

print(out, row.names = FALSE, digits = 4)
write.table(out, file.path(OUT_DIR, "figS3B_development_cpm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
