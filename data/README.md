# Frozen inputs

The GSEA of Fig. S2E needs two inputs that come from sources carrying no version number, so a
rerun a year later would silently score a different gene set. Both are pinned here by sha256.

The orthologue table is in the repository. The gene set library is **not**: Harmonizome is
distributed under CC BY-NC-SA 4.0, which does not allow redistribution from an MIT licensed
repository. Download it with `bash data/fetch_disgenet.sh`, which checks the checksum and warns
if the upstream file has changed.

## `disgenet_harmonizome.gmt.gz`  (not in the repository, see above)

DisGeNET disease gene sets, Harmonizome "crisp" library, downloaded on 30 September 2026 from
`https://maayanlab.cloud/static/hdfs/harmonizome/data/disgenetdisease/gene_set_library_crisp.gmt.gz`.

sha256 `05fcf5dcd8df7b7b79a8e16e23f50dedd67c482cb9a86a80b580faa05da06256`.

7,786 gene sets. "Hypertrophic Cardiomyopathy" has 516 genes, the same count as the run that
produced the published panel.

## `mouse_human_orthologues.tsv`  (in the repository)

Mouse to human mapping used by the GSEA of Fig. S2E: Ensembl gene id of the mouse gene and the
associated human gene symbol, 20,137 rows. It is the same source and the same rule as the analysis
that produced the panel, which queried Ensembl BioMart through `biomaRt` on the day it ran: where
a mouse gene has several human homologues the first one returned is kept, and the join onto the
differential expression table is on the gene id, not on the symbol.

sha256 `b7c16853e0abf0c218ff036f0437fa4bd74eed8cf8dc26e92040cdc678fe652d`.

Rebuild it with `bash fetch_orthologues.sh`. Two things to know before doing so.

**`www.ensembl.org` refuses the BioMart endpoint**, it answers 308 with an empty body, which is
why the `biomaRt` package reports that it cannot query any Ensembl site. The mirrors answer, so
the script uses `useast.ensembl.org`; `asia.ensembl.org` works too and `BIOMART` overrides it.

**A newer Ensembl release changes the result.** The frozen copy was built from release 116 and
reproduces the published panel closely but not exactly, because the original run used the release
current in March 2026:

| | published run | this table |
|---|---|---|
| genes in the ranked list | 11,603 | 11,580 |
| genes of the set found in it | 428 | 426 |
| NES | 1.50 | 1.52 |
| nominal P | 1.1e-04 | 5.9e-05 |
| FDR over all 7,786 sets | 0.030 | 0.019 |
| leading edge | 108 genes | 109 genes, 107 shared |

The published run itself is not reproducible to the digit: `fgsea` was called on 111 processes
with no seed, and three runs of it on the same day gave FDR 0.030, 0.044 and 0.045. This script
fixes the seed and uses one process, so its own number does not move.

## `gse176092/`

Four tables that the GSE176092 chain reads. They are not code and the MIT licence of this
repository does not extend to them: each keeps the terms of the source it comes from.

| file | what it holds | source |
|---|---|---|
| `wgcna_modules.txt` | the twelve WGCNA modules M1 to M12, gene symbols | **not in the repository**, see below |
| `manuscript_cells.tsv` | library, barcode, cell type and UMAP of the snRNA-seq cells | the analysis of this paper |
| `tangram_map_manuscript.tsv` | the Tangram map, spots by cell type | the analysis of this paper |
| `tangram_celltype_counts_manuscript.tsv` | counts summed per cell type, genes by cell type | the analysis of this paper |

### `wgcna_modules.txt`  (not in the repository)

The twelve WGCNA modules are published as supplementary tables of the paper that deposited
GSE176092: S. Yamada et al., *Spatiotemporal transcriptome analysis reveals critical roles for
mechano-sensing genes at the border zone in remodeling after myocardial infarction*, Nature
Cardiovascular Research 1, 1072 to 1083 (2022). They belong to that paper and are left for the
reader to download from it.

Build the file as twelve tab separated columns named `M1` to `M12`, one gene symbol per cell,
padded with empty cells to the longest module. The copy used here has 1,572 rows and sha256
`9b8b79b842743badd58b97b62054583de7a616e3c58a486e2650d32e0ada7243`. `CHH_WGCNA` points at it.

The doublet calls and the Tangram fit were run without a seed, so the last three are kept as
they were rather than refitted. `CHH_CELLS=none` and `CHH_TANGRAM_MAP=none` compute them again,
giving figures that are close to the published ones.
