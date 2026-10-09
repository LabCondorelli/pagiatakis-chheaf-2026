#!/usr/bin/env Rscript
## snRNA-seq reference from the raw GEO matrices: decontX, QC filters of the original analysis,
## then doublets, cell types and UMAP from MANUSCRIPT_CELLS (original doublet calls had no seed).
## MANUSCRIPT_CELLS = "none": scDblFinder and re-annotation, close to but not the published one.
## Usage: Rscript prep2_snrnaseq_reference.R
## Output: work/snrnaseq_reference.rds, results/prep2_cell_counts.tsv
##         (results/prep2_cluster_markers.tsv with MANUSCRIPT_CELLS = "none")

suppressPackageStartupMessages({
   library(Seurat); library(SingleCellExperiment); library(scDblFinder); library(celda)
})
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "..", "00_config.R"))
source(file.path(here, "..", "utils", "sc_helpers.R"))
options(future.globals.maxSize = 1e16)

LIBRARIES <- c(
   "GSM5355657_nonCM_Sham_1", "GSM5355658_nonCM_Sham_2", "GSM5355659_nonCM_MI_day1",
   "GSM5355660_nonCM_MI_day1", "GSM5355661_nonCM_MI_day7", "GSM5355662_nonCM_MI_day7",
   "GSM5943175_CM_Sham_1", "GSM5943176_CM_Sham_2", "GSM5943177_CM_MI_day1_IZ+BZ_1",
   "GSM5943178_CM_MI_day1_IZ+BZ_2", "GSM5943179_CM_MI_day1_RZ_1", "GSM5943180_CM_MI_day1_RZ_2",
   "GSM5943181_CM_MI_day7_IZ+BZ_1", "GSM5943182_CM_day7_IZ+BZ_2", "GSM5943183_CM_MI_day7_RZ_1",
   "GSM5943184_CM_MI_day7_RZ_2", "GSM5943185_CM_MI_day14_IZ+BZ_1", "GSM5943186_CM_MI_day14_IZ+BZ_2",
   "GSM5943187_CM_MI_day14_RZ_1", "GSM5943188_CM_MI_day14_RZ_2")

## harmony.snn_res.0.4 cluster -> cell type (MANUSCRIPT_CELLS = "none").
CELL_TYPE_MAP <- c("0" = "Cardiomyocyte", "1" = "Cardiomyocyte", "2" = "Cardiomyocyte",
                   "3" = "Endothelial_Cell", "4" = "Fibroblast", "5" = "Myeloid_Cell",
                   "6" = "Myeloid_Cell", "7" = "Pericyte", "8" = "Fibroblast",
                   "9" = "Endothelial_Cell", "10" = "Fibroblast", "11" = "Fibroblast",
                   "12" = "Oxphos_Cell", "13" = "Endothelial_Cell", "14" = "Glia_Cell")
CM_SUBTYPE_MAP <- c("0" = "Cardiomyocyte_1", "1" = "Cardiomyocyte_2", "2" = "Cardiomyocyte_3")
FROM_TABLE <- MANUSCRIPT_CELLS != "none"

## ---- loading and ambient RNA ---------------------------------------------------------------
srt <- lapply(LIBRARIES, function(l)
   srt_from_featuremtx_barcode_folder(file.path(RAW_DIR, l), sample.name = sub("^GSM\\d+_", "", l)))
srt <- merge(srt[[1]], srt[-1])
srt[["RNA"]] <- as(srt[["RNA"]], Class = "Assay")
srt <- lapply(SplitObject(srt, split.by = "orig.ident"), decontX_remove)
srt <- merge(srt[[1]], srt[-1])
n_loaded <- ncol(srt)

## ---- QC and doublets -------------------------------------------------------------------------
srt <- filter_zero_expr_genes(srt, assay.oi = "Decont")
srt <- filter_zero_expr_genes(srt, assay.oi = "RNA")
srt <- add_all_qc_metrics(srt, assay.oi = "RNA")
srt <- add_all_qc_metrics(srt, assay.oi = "Decont")
srt <- subset(srt, percent.mito < 60)
srt <- subset(srt, log10UMI > 0.82)
## merge() suffixes barcodes by position: key on library|barcode.
key <- paste(srt$orig.ident, sub("_[0-9]+$", "", colnames(srt)), sep = "|")
if (FROM_TABLE) {
   pub <- read.delim(MANUSCRIPT_CELLS)
   pub_key <- paste(pub$library, pub$barcode, sep = "|")
   if (!all(pub_key %in% key)) stop(sum(!pub_key %in% key), " published cells fail the QC above")
   srt <- srt[, key %in% pub_key]
} else {
   doublets <- scDblFinder_clusters(srt, nrep = 3, samples = "orig.ident", assay.oi = "Decont")
   srt <- srt[, !colnames(srt) %in% doublets]
}
srt <- subset(srt, nCount_Decont > 500)
n_qc <- ncol(srt)

srt$timepoint <- ifelse(grepl("day\\d+_", srt$orig.ident),
                        sub(".*(day\\d+)_.*", "\\1", srt$orig.ident), "Sham")

## ---- cardiomyocyte-enriched libraries: annotation ---------------------------------------------
ref <- subset(srt, cells = colnames(srt)[!grepl("nonCM", srt$orig.ident)]); rm(srt); gc()
ref <- subset(ref, percent.ribo < 10)
if (FROM_TABLE) {
   if (ncol(ref) != nrow(pub)) stop(nrow(pub) - ncol(ref), " published cells removed by the UMI or ribosomal filter")
   i <- match(paste(ref$orig.ident, sub("_[0-9]+$", "", colnames(ref)), sep = "|"), pub_key)
   ref$cell.type <- pub$cell_type[i]
   umap <- as.matrix(pub[i, c("umap_1", "umap_2")]); dimnames(umap) <- list(colnames(ref), c("umap_1", "umap_2"))
   ref[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "umap_", assay = "Decont")
} else {
   ref <- preprocess_Normalize(ref, assay = "Decont", nfeatures = 2000)
   ref <- harmonize_Umap(ref, harmonize = "orig.ident", min.dist = 0.3, n.neighbors = 50)
   DefaultAssay(ref) <- "Decont"
   ref <- SCTransform(ref, vst.flavor = "v2", method = "glmGamPoi", assay = "Decont", verbose = FALSE)
   ref <- correct_sct(ref)

   Idents(ref) <- ref$harmony.snn_res.0.4
   markers <- FindAllMarkers(ref, assay = "SCT", slot = "data", min.pct = 0.4, min.diff.pct = 0.3, only.pos = TRUE)
   write.table(markers, file.path(OUT_DIR, "prep2_cluster_markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
   if (!setequal(levels(Idents(ref)), names(CELL_TYPE_MAP)))
      stop("clusters found (", paste(levels(Idents(ref)), collapse = ","), ") differ from CELL_TYPE_MAP: ",
           "redo the map from results/prep2_cluster_markers.tsv")
   ref$cell.type <- unname(CELL_TYPE_MAP[as.character(ref$harmony.snn_res.0.4)])
   ref <- subset(ref, cell.type != "Oxphos_Cell")
   ref <- harmonize_Umap(ref, harmonize = "orig.ident", clustering = FALSE, do_harmony = FALSE,
                         min.dist = 0.55, n.neighbors = 40)
   cm <- as.character(ref$harmony.snn_res.0.4) %in% names(CM_SUBTYPE_MAP)
   ref$cell.type[cm] <- CM_SUBTYPE_MAP[as.character(ref$harmony.snn_res.0.4[cm])]
}

DefaultAssay(ref) <- "RNA"
if (inherits(ref[["RNA"]], "Assay5")) {
   ref <- JoinLayers(ref, assay = "RNA")
   ref[["RNA"]] <- as(ref[["RNA"]], Class = "Assay")
}
ref <- NormalizeData(ref, assay = "RNA")
saveRDS(ref, REFERENCE_RDS)

counts <- as.data.frame(table(cell_type = ref$cell.type))
write.table(rbind(data.frame(cell_type = c("loaded (20 libraries)", "after QC and doublets"),
                             Freq = c(n_loaded, n_qc)), counts),
            file.path(OUT_DIR, "prep2_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("reference: %d cells (%d loaded, %d after QC), cells from %s -> %s\n", ncol(ref), n_loaded, n_qc,
            if (FROM_TABLE) "the published table" else "scDblFinder", REFERENCE_RDS))
print(counts, row.names = FALSE)
