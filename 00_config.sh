## Shared paths for the shell scripts. Replace the placeholders with your own copies.
## Every path can also be set from the environment.

: "${GTF_MM10:=<PATH_TO>/gencode.vM25.annotation.gtf}"        # GENCODE vM25, mouse
: "${BEDTOOLS:=bedtools}"                                      # bedtools v2.31.1
: "${OUT_DIR:=${CHH_OUT:-results}}"           # CHH_OUT is the name the R scripts use

## ---- single nuclei and spatial transcriptomics of the infarcted heart (GSE176092) ----
: "${CHH_RAW:=<PATH_TO>/GSE176092}"        # raw GEO files, one folder per sample
: "${CHH_WORK:=work}"                      # intermediate objects, about 8 GB
: "${CHH_WGCNA:=<PATH_TO>/wgcna_modules.txt}"
: "${CHH_CELLS:=<PATH_TO>/manuscript_cells.tsv}"                            # "none" = scDblFinder
: "${CHH_TANGRAM_MAP:=<PATH_TO>/tangram_map_manuscript.tsv}"                # "none" = fit again
: "${CHH_TANGRAM_COUNTS:=<PATH_TO>/tangram_celltype_counts_manuscript.tsv}"
: "${R_ENV:=scrnaseq}"                     # conda env of the R steps (sessionInfo_gse176092_R.txt)
: "${TANGRAM_ENV:=tangram}"                # conda env of prep4 (sessionInfo_gse176092_python.txt)
export CHH_RAW CHH_WORK CHH_OUT CHH_WGCNA CHH_CELLS CHH_TANGRAM_MAP CHH_TANGRAM_COUNTS

mkdir -p "$OUT_DIR"
