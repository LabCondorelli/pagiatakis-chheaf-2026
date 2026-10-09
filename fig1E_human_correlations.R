#!/usr/bin/env Rscript
## Fig. 1E - Spearman correlations between RTN4 (NOGO-A) and the other transcripts of the
## locus, computed separately in non-failing (NF) and hypertrophic (HCM) human left ventricles.
## The same model also gives the differential expression of syntCHHEAF shown in Fig. R1D.
##
## Design: group + age + sex. Correlations are computed on log2 CPM
## with the age and sex effects subtracted, so that only the group term is left in the matrix.
##
## Usage: Rscript fig1E_human_correlations.R
## Output: results/fig1E_correlations.tsv, results/fig1E_differential_expression.tsv

suppressPackageStartupMessages({ library(edgeR); library(limma) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_config.R"))

## Five unstranded non-failing libraries. A strand specific count is not defined for them and
## syntCHHEAF is antisense to RTN4, so they are excluded.
UNSTRANDED <- c("nf_rep18", "nf_rep37", "nf_rep139", "nf_rep142", "nf_rep143")

## Transcripts of the locus. syntCHHEAF_CLS is the GENCODE CLS model anchTM_000000915718.
TX <- c(syntCHHEAF  = "syntCHHEAF_CLS",
        NogoA       = "ENST00000337526",
        NPPA        = "ENST00000376480",
        NogoA_trunc = "ENST00000461004",
        lncRTN4     = "ENST00000834600")

## ---- counts -----------------------------------------------------------------------------
## featureCounts writes its command line on line 1 and names the columns after the BAM paths.
read_fc <- function(f) {
  x <- read.delim(f, skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
  cols <- grep("[.]bam$", colnames(x), value = TRUE)
  m <- as.matrix(x[, cols, drop = FALSE])
  rownames(m) <- x[[1]]
  colnames(m) <- basename(dirname(cols))   # one directory per sample
  m
}
g  <- read_fc(file.path(COUNTS_DIR, "RNA_genes_counts"))
tx <- read_fc(file.path(COUNTS_DIR, "RNA_transcripts_counts"))

smp <- setdiff(colnames(tx), UNSTRANDED)
g <- g[, smp, drop = FALSE]; tx <- tx[, smp, drop = FALSE]
grp <- ifelse(grepl("^hcm", smp), "HCM", "NF")

## ---- design -----------------------------------------------------------------------------
meta <- read.delim(SAMPLE_TABLE, stringsAsFactors = FALSE)
meta <- meta[match(smp, meta$sample), ]
stopifnot(!anyNA(meta$age), !anyNA(meta$sex))
design <- cbind(hcm     = as.numeric(grp == "HCM"),
                nf      = as.numeric(grp == "NF"),
                age     = as.numeric(meta$age),
                sexmale = as.numeric(meta$sex == "male"))
rownames(design) <- smp

## ---- normalisation ----------------------------------------------------------------------
## TMM factors and library sizes come from the gene level matrix and are carried over to the
## transcript level one, so that both are on the same scale.
dg  <- calcNormFactors(DGEList(counts = g)[filterByExpr(DGEList(counts = g), group = grp,
                                                        min.count = 5), , keep.lib.sizes = FALSE])
dt  <- DGEList(counts = tx)
dt$samples$lib.size     <- dg$samples$lib.size
dt$samples$norm.factors <- dg$samples$norm.factors
d2  <- dt[filterByExpr(dt, group = grp, min.count = 2), , keep.lib.sizes = TRUE]

## ---- fit --------------------------------------------------------------------------------
v   <- voom(d2, design = design)
fit <- lmFit(v, design)
cov <- design[, c("age", "sexmale"), drop = FALSE]
b   <- fit$coefficients[, c("age", "sexmale"), drop = FALSE]; b[is.na(b)] <- 0
M   <- as.matrix(v$E) - b %*% t(cov)          # log2 CPM without the age and sex effects

tt <- topTable(eBayes(contrasts.fit(fit, makeContrasts(hcm - nf, levels = design)), robust = TRUE),
               number = Inf, sort.by = "none")

## ---- output -----------------------------------------------------------------------------
id  <- function(nm) which(sub("[.].*$", "", rownames(M)) == TX[nm])[1]
val <- function(nm) as.numeric(M[id(nm), ])
sp  <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman"))
## Fisher z test for the difference between two independent correlations
fz  <- function(r1, n1, r2, n2) 2 * pnorm(-abs((atanh(r1) - atanh(r2)) / sqrt(1/(n1-3) + 1/(n2-3))))

pairs <- list(c("NogoA", "NPPA"), c("NogoA", "syntCHHEAF"),
              c("NogoA", "lncRTN4"), c("NogoA", "NogoA_trunc"))
co <- do.call(rbind, lapply(pairs, function(p) {
  x <- val(p[1]); y <- val(p[2])
  h <- sp(x[grp == "HCM"], y[grp == "HCM"]); n <- sp(x[grp == "NF"], y[grp == "NF"])
  data.frame(pair = paste(p, collapse = "_vs_"),
             rho_HCM = unname(h$estimate), p_HCM = h$p.value,
             rho_NF  = unname(n$estimate), p_NF  = n$p.value,
             p_diff  = fz(h$estimate, sum(grp == "HCM"), n$estimate, sum(grp == "NF")))
}))

rt <- sub("[.].*$", "", rownames(tt))
de <- do.call(rbind, lapply(names(TX), function(nm) {
  j <- which(rt == TX[nm])[1]
  data.frame(name = nm, logFC = tt$logFC[j], P = tt$P.Value[j], FDR = tt$adj.P.Val[j])
}))

cat(sprintf("samples: %d (%d HCM, %d NF); transcripts tested: %d\n",
            length(smp), sum(grp == "HCM"), sum(grp == "NF"), nrow(d2)))
print(co, row.names = FALSE, digits = 3)
print(de, row.names = FALSE, digits = 3)

write.table(co, file.path(OUT_DIR, "fig1E_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(de, file.path(OUT_DIR, "fig1E_differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(M = M, design = design, counts = d2$counts, group = grp, samples = smp),
        file.path(OUT_DIR, "human_matrix.rds"))
