# Analysis code for Pagiatakis et al., *Chheaf*

One script per figure panel that rests on a computational analysis.

The code is under the MIT licence, in `LICENSE`. The licence covers the code. Two inputs belong
to other people and are left out of the repository, with instructions to fetch them in
`data/README.md`: the DisGeNET gene set library and the WGCNA modules of Yamada et al.

## Authorship

The single nuclei and spatial transcriptomics of GSE176092 were written by Lorenzo Marcone: the
chain in `gse176092/`, the three panels it prepares (Fig. 1D, Fig. S1E, Fig. S1F) and the
helpers in `utils/`. The rest is by Simone Serio.

## How to run

Paths are placeholders. Set them either by editing `00_config.R` and `00_config.sh`, or from
the environment:

```bash
export CHH_COUNTS=/data/magnet/counts/featureCounts
export CHH_SAMPLES=/data/magnet/batchesV2.txt
Rscript fig1E_human_correlations.R
```

Every script writes to `results/` and prints the numbers it produced, so a run can be checked
against the paper without opening the output files. `expected/` holds the same tables as we
obtained them, for a direct diff. Recomputing the expression matrices from the counts moves the
last decimals of the panels that plot them; the gene sets and the calls do not move.

The gene set library used by the GSEA is not in this repository: Harmonizome is CC BY-NC-SA 4.0,
which does not allow redistribution from an MIT licensed repository. Fetch it once with
`bash data/fetch_disgenet.sh`, which verifies the checksum of the version used in the paper.
The mouse to human orthologue table it is scored against is in the repository and can be
rebuilt with `bash data/fetch_orthologues.sh`; a newer Ensembl release changes the enrichment,
see `data/README.md`.

Requirements: R with `edgeR`, `limma`, `ComplexHeatmap`, `ggplot2`, `ggpubr`, `fgsea`;
Python with `pandas`, `numpy` and `pyBigWig`; `featureCounts`, `bedtools`, `liftOver`.
These panels need two further environments, one for R with `Seurat`, `SingleCellExperiment`,
`scDblFinder`, `celda`, `harmony`, `Nebulosa`, `sceasy`, `magick` and `patchwork`, one for Python
with `scanpy`, `tangram-sc` and `torch`. The versions of all four are in `sessionInfo_R.txt`
and `sessionInfo_python.txt` here and in `gse176092/`.

## Inputs

The expression matrices are not in this repository. What is deposited is the raw sequencing
data, at the accessions below; the matrices are obtained from it with the quantification
parameters stated in the Methods, and for the human cohort with `00_quantify_human.sh`. Every
script states the file it needs at the top.

| script | dataset | accession |
|---|---|---|
| `00_quantify_human.sh`, `fig1E_human_correlations.R`, `fig1E_coverage_tracks.py` | human left ventricle, 118 non-failing and 26 hypertrophic | PRJNA595151 |
| `fig1E_human_correlations.R`, transcript model | GENCODE CLS, adult human heart | anchTM_000000915718, Zenodo 15004659 |
| `00_de_and_cpm.R`, `00_export_mouse_tac.R`, `fig1A_*`, `fig1B_*`, `figS1A_*`, `dataS2_*` | mouse cardiomyocytes, 1 to 7 days after TAC | PRJNA1372455 |
| `fig2D_S2B_*`, `figS2C_*`, `figS2E_*`, `dataS1_*` | Chheaf wild-type and knockout cardiomyocytes | PRJNA1372682 |
| `figS1D_mi_expression.R` | mouse hearts after myocardial infarction | PRJNA227375 |
| `figS3B_development_cpm.R` | cardiomyocytes at E14.5, neonatal, adult and TAC | PRJNA278078 |
| `fig3B_tf_cobound_liftover.sh` | Gata4, Nkx2-5 and Tead1 in fetal and adult heart | PRJNA510592 |
| `fig1D_*`, `figS1E_*`, `figS1F_*`, `gse176092/` | single nuclei and spatial transcriptomics of the infarcted heart | GSE176092 |
| Fig. 3B, promoter capture Hi-C arcs | iPSC and iPSC-derived cardiomyocytes, Montefiori et al. | E-MTAB-6014 |
| Fig. 3B, human H3K27ac | failing and non-failing human left ventricle, Gilsbach et al. | PRJNA353755 |

`00_de_and_cpm.R` turns any of these featureCounts matrices into the tables the panels read.
The knockout was quantified and analysed twice, against a protein coding and a long non-coding
annotation, so reproducing its tables means running the script once per annotation; the script
header says why subsetting one combined matrix is not equivalent.

`fig1E_synteny_blocks.sh` needs no deposited data: it queries the UCSC chain alignment track.

## Order

Four scripts prepare inputs that the panel scripts read. Run them first.

| script | produces |
|---|---|
| `00_quantify_human.sh` | featureCounts matrices for the human left ventricles |
| `00_de_and_cpm.R` | log2 CPM matrix and differential expression table of a mouse dataset |
| `00_export_mouse_tac.R` | the same two tables for the TAC time course, from our pipeline output |
| `dataS2_00_neighbours.sh` | nearest coding gene on each side of every lncRNA |

`00_de_and_cpm.R` is the differential expression step of the Methods: from a featureCounts
matrix, a sample table and a list of contrasts it writes the log2 CPM matrix, the differential
expression table and the TPM matrix that the mouse panels read. Run it on the TAC time course
with the prefix `results/mouse_tac` and Fig. 1A, Fig. 1B and Fig. S1A find their inputs
directly; `CHH_TAC_DE_ALL=results/mouse_tac_DE_calls.tsv` points Data S2 at the same table.
On that dataset it reproduces our own pipeline: the same 15,565 genes pass the filter, the log2
CPM values differ by a median of 0.0002, and the differential calls and the 40 lincRNA of
Fig. 1A are identical. `00_export_mouse_tac.R` is the shortcut for whoever already has the
output of that pipeline.

The three panels that come from GSE176092 start from GEO instead of from our matrices, and have
their own chain, in `gse176092/`. Run it with `bash gse176092/run_all.sh [--skip-download]`, or
step by step:

| script | produces | reads |
|---|---|---|
| `00_download_geo.sh` | the raw files, 20 snRNA-seq and 14 Visium | GEO |
| `prep1_visium_object.R` | the Visium object | the raw Visium |
| `prep2_snrnaseq_reference.R` | the snRNA-seq reference | the raw snRNA-seq, `manuscript_cells.tsv` |
| `prep3_tangram_inputs.R` | the Tangram inputs | prep1, prep2 |
| `prep4_tangram_projection.py` | the projection of the 19 genes | prep3, `tangram_*_manuscript.tsv` |

prep1 and prep2 need about 200 GB and several hours, prep4 up to a day when the Tangram map is
refitted. The R steps and prep4 need two different environments, named in `00_config.sh`;
their versions are in `gse176092/sessionInfo_R.txt` and `gse176092/sessionInfo_python.txt`.

## Panels

| panel | script | what it computes |
|---|---|---|
| Fig. 1A | `fig1A_lincRNA_heatmap.R` | Z-scores of the 40 lincRNA modulated by TAC, split by direction |
| Fig. 1B | `fig1B_locus_timecourse.R` | log2 CPM of Chheaf and Rtn4 over the time course |
| Fig. 1D | `fig1D_spatial_chheaf_rtn4.R` | Rtn4 and Chheaf on the sections, Tangram projection |
| Fig. 1E | `fig1E_synteny_blocks.sh` | human/mouse syntenic DNA blocks, the orange bars |
| Fig. 1E | `fig1E_coverage_tracks.py` | strand resolved coverage in non-failing and HCM hearts |
| Fig. 1E | `fig1E_human_correlations.R` | the four Spearman correlations, per group |
| Fig. 2D | `fig2D_S2B_ko_timecourse.R` | TPM of Chheaf and Rtn4 in WT and KO hearts |
| Fig. 3B | `fig3B_tf_cobound_liftover.sh` | TF co-bound regions projected onto hg38 |
| Fig. S1A | `figS1A_lncrna_biotypes.R` | biotype composition of the modulated lncRNA |
| Fig. S1C | `dataS2_lincRNA_selection.py` | the lincRNA and flanking gene pairs the panel plots, from the chain that writes Data S2 |
| Fig. S1D | `figS1D_mi_expression.R` | Chheaf and Rtn4 after myocardial infarction |
| Fig. S1E | `figS1E_snrnaseq_chheaf_rtn4.R` | snRNA-seq UMAP, module scores, Nebulosa densities |
| Fig. S1F | `figS1F_spatial_lincrna_pairs.R` | co-expression of the other lincRNA pairs |
| Fig. S2B | `fig2D_S2B_ko_timecourse.R` | Spearman between Rtn4 and Chheaf at 1 week |
| Fig. S2C | `figS2C_de_summary.R` | the genes uniquely down-regulated at 1 week |
| Fig. S2E | `figS2E_gsea_disgenet.R` | GSEA of the DisGeNET hypertrophic cardiomyopathy set |
| Fig. S3B | `figS3B_development_cpm.R` | Chheaf and Rtn4 across development |
| Data S1 | `dataS1_ko_vs_wt_de.R` | TAC KO versus TAC WT at 1, 2 and 4 weeks |
| Data S2 | `dataS2_lincRNA_selection.py` | the 40 lincRNA, their neighbours and the selection rule |

The four files the GSE176092 chain reads are in `data/gse176092/`: the WGCNA modules of Yamada
et al., and the published cell identities, Tangram map and per-cell-type counts. The original
doublet calls and the Tangram fit had no seed, so those three come frozen from the published
analysis; everything else is recomputed. Setting `CHH_CELLS=none` and `CHH_TANGRAM_MAP=none` recomputes them too,
giving figures that are close but not identical.
