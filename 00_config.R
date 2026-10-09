## Shared paths for the R scripts. Replace the placeholders with your own files.
## Every path can also be set from the environment, e.g. CHH_COUNTS=/data/counts Rscript fig1E_*.R
## Values that belong to one script only are declared in that script, not here.
##
## The placeholders name each file by what it has to contain, not by where ours happens to sit.

p <- function(var, default) { v <- Sys.getenv(var); if (nzchar(v)) v else default }

## GENCODE vM25, mouse. Same variable name as in 00_config.sh, so one setting serves both.
GTF_MM10 <- p("GTF_MM10", "<PATH_TO>/gencode.vM25.annotation.gtf")

## Human left ventricle, GENCODE v47 plus the CLS model: the directory 00_quantify_human.sh
## writes, holding RNA_genes_counts and RNA_transcripts_counts
COUNTS_DIR   <- p("CHH_COUNTS",  "<PATH_TO>/human_counts")

## Human sample table, tab separated, columns: sample, age, race, sex
SAMPLE_TABLE <- p("CHH_SAMPLES", "<PATH_TO>/human_samples.tsv")

## Mouse TAC time course (sham and TAC at 1, 2, 4 and 7 days), GENCODE vM25.
## log2 CPM, TMM normalised. Columns: ensgID,
## symbol, one per sample. Whether the table is filtered by expression changes nothing here:
## the filter only drops rows, the values of the genes that stay are the same.
MOUSE_TAC_CPM     <- p("CHH_TAC_CPM", "<PATH_TO>/mouse_tac_log2cpm.tsv")

## Same time course, differential calls coded -1 / 0 / +1 in the Results.* columns, one per
## time point, already thresholded at FDR < 0.05 and |log2FC| > 0.5
MOUSE_TAC_DE_SIGN <- p("CHH_TAC_DE",  "<PATH_TO>/mouse_tac_DE_sign.tsv")

## Mouse development and hypertrophy (E14.5, neonatal, adult, TAC), featureCounts gene matrix
DEV_COUNTS <- p("CHH_DEV_COUNTS", "<PATH_TO>/mouse_development_counts.txt")

## Mouse Chheaf-KO (sham and TAC at 1, 2 and 4 weeks), limma-voom tables with all genes and
## their statistics, protein coding and lncRNA trees analysed separately
KO_DE_ALL <- p("CHH_KO_DE_ALL", "<PATH_TO>/ko_DE_protein_coding.txt")
KO_DE_LNC <- p("CHH_KO_DE_LNC", "<PATH_TO>/ko_DE_lncRNA.txt")

## Same experiment, TPM matrix, lncRNA included
KO_TPM    <- p("CHH_KO_TPM",    "<PATH_TO>/ko_TPM.txt")

## Mouse myocardial infarction (4 sham, 4 infarct): log2 CPM matrix and the limma-voom table
MI_CPM <- p("CHH_MI_CPM", "<PATH_TO>/mi_log2cpm.txt")
MI_DE  <- p("CHH_MI_DE",  "<PATH_TO>/mi_DE.txt")

## ---- single nuclei and spatial transcriptomics of the infarcted heart (GSE176092) ----
## GEO GSE176092, one folder per sample as written by 00_download_geo.sh.
RAW_DIR <- p("CHH_RAW", "<PATH_TO>/GSE176092")

## WGCNA modules M1-M12: one column per module, gene symbols, tab-separated.
WGCNA_MODULES <- p("CHH_WGCNA", "<PATH_TO>/wgcna_modules.txt")

## Published snRNA-seq cells: library, barcode, cell type, UMAP ("none" = scDblFinder).
MANUSCRIPT_CELLS <- p("CHH_CELLS", "<PATH_TO>/manuscript_cells.tsv")

## Intermediate objects (~8 GB).
WORK_DIR <- p("CHH_WORK", "work")
VISIUM_RDS    <- file.path(WORK_DIR, "visium_object.rds")
REFERENCE_RDS <- file.path(WORK_DIR, "snrnaseq_reference.rds")
TANGRAM_DIR   <- file.path(WORK_DIR, "tangram")
TANGRAM_GENES <- file.path(TANGRAM_DIR, "tangram_projection_19genes.csv")

## ---- shared by the spatial panels (Fig. 1D, Fig. S1F) -------------------------------------

## The four sections of Fig. 1D, in column order.
PUB_SLIDES <- c("GSM5355663_WT_Sham", "GSM5943190_WT_MI_day1_2",
                "GSM5943192_WT_MI_day7_2", "GSM5943195_WT_MI_day14_3")

## Day 1 is stored upside down.
SLIDE_ROTATION <- c(GSM5355663_WT_Sham = 0, GSM5943190_WT_MI_day1_2 = 180,
                    GSM5943192_WT_MI_day7_2 = 0, GSM5943195_WT_MI_day14_3 = 0)

## Pairs of Fig. S1C; all twelve set the shared colour scale.
PAIRS <- data.frame(
   linc = c("Gm12092", "Gm13481", "Gm13481", "Gm5086", "Gm30085", "Gm30085",
            "A730036I17Rik", "A730036I17Rik", "Gm17501", "Gm17501", "Gm42047", "Gm42047"),
   gene = c("Rtn4", "Zeb2", "Acvr2a", "Enc1", "Pdzrn4", "Gxylt1",
            "Slc20a1", "Ckap2l", "Ccn1", "Ddah1", "Klhdc4", "Slc7a5"),
   stringsAsFactors = FALSE)

## Co-expression cap: 2.6 = both genes in their top 5%; higher values are squished.
COEXPR_CAP <- 2.6

OUT_DIR <- p("CHH_OUT", "results")
for (d in c(OUT_DIR, WORK_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
