#!/usr/bin/env Rscript
## From a featureCounts matrix to the differential expression table, the log2 CPM matrix and the
## TPM matrix that the panel scripts read. This is the differential expression step of the
## Methods, for the four mouse datasets.
##
## Model: edgeR filterByExpr and TMM, then limma voomWithQualityWeights, lmFit and
## eBayes(robust = TRUE), which is what produced every mouse table of the paper. A gene is called
## at FDR < 0.05 and |log2FC| > 0.5, coded -1 / 0 / +1 in the Results columns.
##
## The log2 CPM is limma's, log2((count + 0.5) / (TMM library size + 1) * 1e6). The prior count
## of 0.5 is a parameter here because it changes the values plotted for a transcript as low as
## Chheaf.
##
## The expression filter sets min.count to half the median library size in millions of reads,
## rounded to the nearest million, a cutoff of 0.5 counts per million: 32 for the TAC time
## course, 29.5 for development and hypertrophy, 30.5 for the knockout, 132 for the infarction
## dataset. CHH_MIN_COUNT = auto computes it, a number overrides it. The value is printed with
## the results because every FDR depends on how many genes were tested.
##
## Run it with the prefix results/mouse_tac and Fig. 1A, Fig. 1B, Fig. S1A and Data S2 find
## their inputs with no further step.
##
## The panels look for their contrasts by name, so the contrast file has to use these:
##   TAC time course  tac_1d_vs_sham_1d, tac_2d_vs_sham_2d, tac_4d_vs_sham_4d, tac_7d_vs_sham_7d
##   Chheaf knockout  tac_wt_vs_sham_wt_1wk, _2wk, _4wk; tac_ko_vs_tac_wt_1wk, _2wk, _4wk; and
##                    tac_ko_vs_tac_wt, the same comparison pooled over the three weeks
## The infarction dataset has one contrast and its name is free.
##
## The knockout was quantified and analysed twice, against a protein coding annotation and
## against a long non-coding one, so the filter, the voom fit and the Benjamini-Hochberg
## universe are those of one group at a time; Fig. S2C, Fig. S2E and Data S1 read the protein
## coding table. Reproduce it by quantifying against that annotation and running this script on
## the result. CHH_BIOTYPES subsets a combined matrix by biotype, which gives the same gene
## universe on counts of its own, since the coding annotation alone assigns a shared read to
## the coding gene and moves 2,854 of the expressed coding genes by more than 5 per cent:
##   protein coding  CHH_BIOTYPES=protein_coding,polymorphic_pseudogene
##   lncRNA          CHH_BIOTYPES=lincRNA,antisense,sense_intronic,sense_overlapping,
##                   bidirectional_promoter_lncRNA,processed_transcript,macro_lncRNA,
##                   3prime_overlapping_ncRNA,TEC
##
## Usage: Rscript 00_de_and_cpm.R <counts> <samples> <contrasts> <prefix>
##   counts     featureCounts output, the one written with --primary and -g gene_id
##   samples    tab separated, columns: sample, group. One row per count column, in any order
##   contrasts  tab separated, columns: name, contrast. For example
##                  tac_1d_vs_sham_1d <TAB> tac.1d - sham.1d
##   prefix     output prefix
## Output: <prefix>_log2cpm.tsv, <prefix>_DE_calls.tsv, <prefix>_TPM.tsv

suppressPackageStartupMessages({ library(edgeR); library(limma) })
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(HERE, "00_config.R"))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 4) stop("usage: Rscript 00_de_and_cpm.R <counts> <samples> <contrasts> <prefix>")
COUNTS <- a[1]; SAMPLES <- a[2]; CONTRASTS <- a[3]; PREFIX <- a[4]

MIN_COUNT   <- p("CHH_MIN_COUNT", "auto")   # "auto" = 0.5 CPM, or an explicit number
PRIOR_COUNT <- as.numeric(p("CHH_PRIOR_COUNT", "0.5"))   # log2 CPM
FDR_MAX     <- as.numeric(p("CHH_FDR",         "0.05"))
LFC_MIN     <- as.numeric(p("CHH_LFC",         "0.5"))
BIOTYPES    <- p("CHH_BIOTYPES", "")        # "" = every gene; see the note on the knockout above

## ---- counts ------------------------------------------------------------------------------
## featureCounts writes its command line on line 1 and names the columns after the BAM paths.
x <- read.delim(COUNTS, skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
cols <- grep("[.]bam$", colnames(x), value = TRUE)
len  <- x$Length                                # union of the exons, for the TPM
m <- as.matrix(x[, cols, drop = FALSE])
rownames(m) <- x[[1]]; names(len) <- x[[1]]
colnames(m) <- basename(dirname(cols))          # one directory per sample

smp <- read.delim(SAMPLES, stringsAsFactors = FALSE)
if (!all(c("sample", "group") %in% names(smp))) stop("the sample table needs columns: sample, group")
## One row per count column, no repetitions: a repeated row would be fitted as an extra replicate.
if (anyDuplicated(smp$sample))    stop("repeated sample in the sample table")
if (anyDuplicated(colnames(m)))   stop("repeated sample in the count matrix")
if (!setequal(colnames(m), smp$sample))
  stop("the sample table and the count matrix do not hold the same samples: ",
       paste(c(setdiff(colnames(m), smp$sample), setdiff(smp$sample, colnames(m))), collapse = ", "))
m <- m[, smp$sample, drop = FALSE]
lev <- unique(make.names(smp$group))
if (length(lev) != length(unique(smp$group)))
  stop("two group names collapse onto one once made syntactic: ", paste(unique(smp$group), collapse = ", "))
grp <- factor(make.names(smp$group), levels = lev)

## ---- annotation, from the GENCODE GTF ------------------------------------------------------
gtf <- read.delim(GTF_MM10, header = FALSE, comment.char = "#", quote = "",
                  colClasses = c("NULL", "NULL", "character", rep("NULL", 5), "character"))
gtf <- gtf[gtf[[1]] == "gene", 2]
get <- function(k) sub(paste0('.*', k, ' "([^"]*)".*'), "\\1", gtf)
ann <- data.frame(gene_id = get("gene_id"), gene_name = get("gene_name"),
                  gene_biotype = get("gene_type"), stringsAsFactors = FALSE)
## Both sides are stripped of the GENCODE version before matching, so a count matrix quantified
## against another patch of the same release still joins.
ann$key <- sub("[.].*$", "", ann$gene_id)
gid_all <- sub("[.].*$", "", rownames(m))
if (anyDuplicated(gid_all)) stop("the count matrix holds the same gene id twice once the version is dropped")
ann <- ann[match(gid_all, ann$key), ]
nomatch <- sum(is.na(ann$gene_id))
if (nomatch > 0.01 * nrow(m))
  stop(sprintf("%s of %s genes are absent from the GTF: wrong annotation release?",
               format(nomatch, big.mark = ","), format(nrow(m), big.mark = ",")))

if (nzchar(BIOTYPES)) {
  sel <- ann$gene_biotype %in% strsplit(BIOTYPES, ",")[[1]]
  cat(sprintf("biotypes %s: %s of %s genes\n", BIOTYPES,
              format(sum(sel, na.rm = TRUE), big.mark = ","), format(nrow(m), big.mark = ",")))
  m <- m[which(sel), , drop = FALSE]; ann <- ann[which(sel), ]; gid_all <- gid_all[which(sel)]
}

## ---- model -------------------------------------------------------------------------------
d <- DGEList(m, group = grp)
if (MIN_COUNT == "auto") MIN_COUNT <- round(median(colSums(m)) / 1e6) * 0.5
MIN_COUNT <- as.numeric(MIN_COUNT)
if (!is.finite(MIN_COUNT) || MIN_COUNT <= 0)
  stop("min.count came out as ", MIN_COUNT, ": libraries too shallow for the rule, set CHH_MIN_COUNT")
keep <- filterByExpr(d, group = grp, min.count = MIN_COUNT)
d <- calcNormFactors(d[keep, , keep.lib.sizes = FALSE])
## filterByExpr also needs a gene to reach min.total.count over all samples, and to clear
## min.count in at least min.prop of the smallest group once that group exceeds large.n. The
## three are printed with min.count because they move the gene universe too.
cat(sprintf("genes: %s quantified, %s kept\n", format(nrow(m), big.mark = ","),
            format(sum(keep), big.mark = ",")))
cat(sprintf(paste("filter: min.count = %s (0.5 CPM of a median library of %.1f M),",
                  "min.total.count = 15, large.n = 10, min.prop = 0.7; smallest group %d\n"),
            format(MIN_COUNT), median(colSums(m))/1e6, min(table(grp))))

X <- model.matrix(~ 0 + grp); colnames(X) <- levels(grp)
v <- voomWithQualityWeights(d, X, plot = FALSE)
fit <- eBayes(lmFit(v, X), robust = TRUE)

cn <- read.delim(CONTRASTS, stringsAsFactors = FALSE)
if (!all(c("name", "contrast") %in% names(cn))) stop("the contrast file needs columns: name, contrast")
cm <- makeContrasts(contrasts = cn$contrast, levels = X)
colnames(cm) <- cn$name
g <- eBayes(contrasts.fit(fit, cm), robust = TRUE)

## ---- output ------------------------------------------------------------------------------
## Both spellings of the gene columns are written because the laboratory tables use
## Genes.ensembl_gene_id / Genes.external_gene_name in one tree and Genes.gene_id /
## Genes.gene_name in the other, and the panel scripts read whichever their tree has.
out <- data.frame(AveExpr = fit$Amean, check.names = FALSE)
for (k in cn$name) {
  tt <- topTable(g, coef = k, number = Inf, sort.by = "none")
  out[[paste0("Coef.", k)]]          <- tt$logFC
  out[[paste0("t.", k)]]             <- tt$t
  out[[paste0("P.value.", k)]]       <- tt$P.Value
  out[[paste0("P.value.adj.", k)]]   <- tt$adj.P.Val
  out[[paste0("Results.", k)]]       <- ifelse(tt$adj.P.Val < FDR_MAX & abs(tt$logFC) > LFC_MIN,
                                               sign(tt$logFC), 0)
}
## The identifiers are written as the bare accession, the form the panel scripts match on. ann
## is aligned to the count matrix row by row, so the genes that passed the filter are its rows
## under keep.
ann_keep <- ann[keep, ]
gid      <- gid_all[keep]
out$Genes.gene_id            <- gid
out$Genes.ensembl_gene_id    <- gid
out$Genes.gene_name          <- ann_keep$gene_name
out$Genes.external_gene_name <- ann_keep$gene_name
out$Genes.gene_biotype       <- ann_keep$gene_biotype
write.table(out, paste0(PREFIX, "_DE_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cpm_log <- cpm(d, log = TRUE, prior.count = PRIOR_COUNT)
write.table(data.frame(ensgID = gid, symbol = ann_keep$gene_name, cpm_log, check.names = FALSE),
            paste0(PREFIX, "_log2cpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## TPM over every quantified gene, stored as log2(TPM + 0.5), the scale Fig. 2D reads.
rpk <- m / (len[rownames(m)] / 1000)
tpm <- log2(sweep(rpk, 2, colSums(rpk) / 1e6, "/") + 0.5)
write.table(data.frame(gene_id = gid_all, symbol = ann$gene_name, tpm, check.names = FALSE),
            paste0(PREFIX, "_TPM.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (k in cn$name)
  cat(sprintf("  %-28s %5s called\n", k, format(sum(out[[paste0("Results.", k)]] != 0), big.mark = ",")))
cat(sprintf("\nwritten: %s_DE_calls.tsv, %s_log2cpm.tsv and %s_TPM.tsv\n", PREFIX, PREFIX, PREFIX))
