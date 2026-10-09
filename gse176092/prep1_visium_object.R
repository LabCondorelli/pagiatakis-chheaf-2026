#!/usr/bin/env Rscript
## Visium object of the 14 sections from the raw GEO files. Scale factors set by hand (GEO has
## none); spots with <= 500 UMIs removed.
## Usage: Rscript prep1_visium_object.R
## Output: work/visium_object.rds, results/prep1_spot_qc.tsv

suppressPackageStartupMessages({ library(Seurat); library(magick) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "..", "00_config.R"))
options(future.globals.maxSize = 1e16)

SLIDES <- c("GSM5355663_WT_Sham", "GSM5355666_WT_MI_day7", "GSM5355668_WT_MI_day14",
            "GSM5943189_WT_MI_day1_1", "GSM5943190_WT_MI_day1_2", "GSM5943191_WT_MI_day1_3",
            "GSM5943192_WT_MI_day7_2", "GSM5943193_WT_MI_day7_3", "GSM5943194_WT_MI_day14_2",
            "GSM5943195_WT_MI_day14_3", "GSM5943196_shCsrp3_MI_day1", "GSM5943197_shCsrp3_MI_day14",
            "GSM5943198_Csrp30E_MI_day1", "GSM5943199_Csrp30E_MI_day14")
MIN_UMI <- 500

load_section <- function(slide) {
   raw <- file.path(RAW_DIR, slide)
   spatial <- file.path(WORK_DIR, "visium_spatial", slide, "spatial")
   dir.create(spatial, recursive = TRUE, showWarnings = FALSE)
   gunzip_one <- function(pattern, dest) {
      gz <- list.files(raw, pattern, full.names = TRUE)
      if (length(gz) != 1) stop(slide, ": expected one ", pattern, " in ", raw)
      R.utils::gunzip(gz, destname = dest, overwrite = TRUE, remove = FALSE)
   }
   gunzip_one("tissue_positions_list\\.csv\\.gz$", file.path(spatial, "tissue_positions_list.csv"))
   gunzip_one("\\.tif\\.gz$", file.path(spatial, "tissue_hires_image.tif"))
   jsonlite::write_json(list(tissue_hires_scalef = 1.0, tissue_lowres_scalef = 0.5,
                             fiducial_diameter_fullres = 87, spot_diameter_fullres = 55),
                        file.path(spatial, "scalefactors_json.json"), auto_unbox = TRUE)
   image_read(file.path(spatial, "tissue_hires_image.tif")) |> image_scale("5%") |>
      image_write(file.path(spatial, "tissue_lowres_image.png"), quality = 10, format = "png")

   obj <- CreateSeuratObject(counts = Read10X(data.dir = raw), assay = "Spatial")
   img <- Read10X_Image(image.dir = spatial, filter.matrix = FALSE)[Cells(obj)]
   DefaultAssay(img) <- "Spatial"
   obj[[slide]] <- img
   obj$orig.ident <- slide
   obj
}

vis <- lapply(SLIDES, load_section)
visium <- merge(vis[[1]], vis[-1]); rm(vis)

qc <- data.frame(section = visium$orig.ident, umi = visium$nCount_Spatial)
qc_tab <- do.call(rbind, lapply(split(qc, qc$section), function(d)
   data.frame(section = d$section[1], spots = nrow(d), removed = sum(d$umi <= MIN_UMI),
              kept = sum(d$umi > MIN_UMI))))
write.table(qc_tab[match(SLIDES, qc_tab$section), ], file.path(OUT_DIR, "prep1_spot_qc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

visium <- subset(visium, nCount_Spatial > MIN_UMI)
visium[["Spatial"]] <- JoinLayers(visium[["Spatial"]])
saveRDS(visium, VISIUM_RDS)
cat(sprintf("Visium object: %d spots x %d genes over %d sections -> %s\n",
            ncol(visium), nrow(visium), length(SLIDES), VISIUM_RDS))
print(qc_tab[match(SLIDES, qc_tab$section), ], row.names = FALSE)
