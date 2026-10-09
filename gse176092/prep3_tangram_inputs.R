#!/usr/bin/env Rscript
## Tangram inputs (h5ad): Visium counts and decontaminated snRNA-seq counts, all genes.
## Usage: Rscript prep3_tangram_inputs.R
## Output: work/tangram/spatial.h5ad, work/tangram/scrna.h5ad

suppressPackageStartupMessages({ library(Seurat) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "..", "00_config.R"))
options(future.globals.maxSize = 1e16)
dir.create(TANGRAM_DIR, showWarnings = FALSE, recursive = TRUE)

to_h5ad <- function(obj, assay, outfile) {
   counts <- GetAssayData(obj, assay = assay, layer = "counts")
   x <- CreateSeuratObject(counts = counts, meta.data = obj@meta.data)
   x[["RNA"]] <- as(x[["RNA"]], Class = "Assay")
   sceasy::convertFormat(x, from = "seurat", to = "anndata", outFile = outfile)
   cat(sprintf("%s: %d x %d genes\n", outfile, ncol(x), nrow(x)))
}

to_h5ad(readRDS(VISIUM_RDS), "Spatial", file.path(TANGRAM_DIR, "spatial.h5ad"))
gc()
ref <- readRDS(REFERENCE_RDS)
if (inherits(ref[["Decont"]], "Assay5")) ref[["Decont"]] <- JoinLayers(ref[["Decont"]])
to_h5ad(ref, "Decont", file.path(TANGRAM_DIR, "scrna.h5ad"))
