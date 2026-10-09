## Helpers of prep2, unchanged from the original analysis.

srt_from_featuremtx_barcode_folder <- function(path.to.folder, sample.name) {
   srt <- CreateSeuratObject(counts = Read10X(path.to.folder))
   srt$orig.ident <- sample.name
   srt
}

decontX_remove <- function(seurat_obj, rounded = TRUE) {
   decontX_results <- celda::decontX(seurat_obj@assays$RNA@counts)
   m <- decontX_results$decontXcounts
   if (isTRUE(rounded)) m <- round(m)
   seurat_obj[["Decont"]] <- CreateAssayObject(m)
   seurat_obj
}

## Threshold scales with the number of genes, not cells.
filter_zero_expr_genes <- function(srt, assay.oi = "RNA") {
   counts <- GetAssayData(object = srt, slot = "counts", assay = assay.oi)
   keep_genes <- Matrix::rowSums(counts > 0) >= round(nrow(counts) / 2000)
   srt[[assay.oi]] <- CreateAssay5Object(counts = counts[keep_genes, ])
   srt
}

add_all_qc_metrics <- function(srt, assay.oi = "RNA") {
   if (class(srt[[assay.oi]]) != "Assay") srt[[assay.oi]] <- as(srt[[assay.oi]], Class = "Assay")
   DefaultAssay(srt) <- assay.oi
   try.err <- try(JoinLayers(srt), silent = TRUE)
   if (class(try.err) != "try-error") srt <- JoinLayers(srt)
   if (assay.oi == "RNA") {
      srt$log10UMI <- log10(srt$nFeature_RNA) / log10(srt$nCount_RNA)
      srt[["percent.mito"]] <- PercentageFeatureSet(srt, features = rownames(srt)[grep(rownames(srt), pattern = "^mt-")])
      srt[["percent.ribo"]] <- PercentageFeatureSet(srt, features = rownames(srt)[grep(rownames(srt), pattern = "^Rp[sl]")])
      srt[["percent_hb"]] <- PercentageFeatureSet(object = srt, pattern = "^Hb[^(p)]")
   }
   srt[[assay.oi]] <- as(srt[[assay.oi]], Class = "Assay5")
   srt
}

scDblFinder_clusters <- function(seurat_obj, nrep = 3, samples = "orig.ident", assay.oi = "RNA") {
   DefaultAssay(seurat_obj) <- assay.oi
   seurat_obj <- DietSeurat(seurat_obj, assays = assay.oi)
   sce <- as.SingleCellExperiment(seurat_obj)
   runs <- replicate(nrep, scDblFinder::scDblFinder(sce, samples = samples))
   calls <- lapply(lapply(runs, as.Seurat, data = NULL), function(x)
      rownames(x@meta.data)[x@meta.data$scDblFinder.class == "doublet"])
   Reduce(intersect, calls)
}

## Human mt/ribo patterns: removes nothing on mouse genes.
preprocess_Normalize <- function(srt, assay, nfeatures = 3000) {
   DefaultAssay(srt) <- assay
   srt[[assay]] <- as(srt[[assay]], Class = "Assay")
   srt[[assay]] <- split(srt[[assay]], f = srt$orig.ident)
   srt <- NormalizeData(srt, assay = assay)
   srt <- FindVariableFeatures(srt, nfeatures = nfeatures, assay = assay)
   mito_ribo <- grep("^MT-|^RP[SL]", rownames(srt), value = TRUE)
   VariableFeatures(srt) <- setdiff(VariableFeatures(srt), mito_ribo)
   srt <- ScaleData(srt, VariableFeatures(srt), assay = assay)
   RunPCA(srt, assay = assay)
}

pc_selector <- function(seurat_object, pc.name) {
   emb <- seurat_object[[pc.name]]@cell.embeddings
   sds <- apply(emb, 2, sd)
   pct <- sds / sum(sds) * 100
   sort(which((pct[-length(pct)] - pct[-1]) > 0.1), decreasing = TRUE)[1] + 1
}

harmonize_Umap <- function(srt, harmonize = "orig.ident", reduction_save = "harmony",
                           umap_name = "umap", min.dist = 0.1, n.neighbors = 50,
                           clustering = TRUE, do_harmony = TRUE) {
   if (do_harmony) {
      srt <- harmony::RunHarmony(srt, harmonize, reduction.save = reduction_save)
      srt <- FindNeighbors(srt, reduction = reduction_save,
                           dims = 1:pc_selector(srt, reduction_save),
                           graph.name = c(paste0(reduction_save, ".nn"), paste0(reduction_save, ".snn")))
   }
   if (clustering)
      srt <- FindClusters(srt, resolution = seq(0.1, 1, 0.1), graph.name = paste0(reduction_save, ".snn"))
   RunUMAP(srt, dims = 1:pc_selector(srt, reduction_save), reduction = reduction_save,
           min.dist = min.dist, n.neighbors = n.neighbors, reduction.name = umap_name)
}

correct_sct <- function(srt, assay.oi = "SCT") {
   vec <- unlist(unname(SCTResults(object = srt, slot = "umi.assay", assay = assay.oi)))
   if (!all(vec == vec[1]))
      for (i in seq_along(vec))
         slot(object = srt@assays[[assay.oi]]@SCTModel.list[[i]], name = "umi.assay") <- assay.oi
   PrepSCTFindMarkers(srt, verbose = FALSE, assay = assay.oi)
}
