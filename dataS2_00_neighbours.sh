#!/bin/bash
## Nearest protein-coding neighbour on each side of every lncRNA in GENCODE vM25.
## Overlapping genes are excluded (-io), so a lincRNA inside a read-through model still gets
## its two real flanking genes.
##
## Usage: bash dataS2_00_neighbours.sh
## Output: results/neighbours_left.txt, results/neighbours_right.txt

set -eo pipefail
source "$(dirname "$0")/00_config.sh"

## GTF to BED: chrom, start-1, end, gene_id, 0, strand, gene_name, gene_type
awk -F'\t' '$3=="gene"{
  match($9,/gene_id "[^"]+"/);   id=substr($9,RSTART+9,RLENGTH-10)
  match($9,/gene_name "[^"]+"/); nm=substr($9,RSTART+11,RLENGTH-12)
  match($9,/gene_type "[^"]+"/); bt=substr($9,RSTART+11,RLENGTH-12)
  print $1"\t"$4-1"\t"$5"\t"id"\t0\t"$7"\t"nm"\t"bt }' "$GTF_MM10" > "$OUT_DIR/genes.bed"

## The lncRNA biotypes present in the differential expression table
awk -F'\t' '$8=="lincRNA"||$8=="antisense"||$8=="processed_transcript"||
            $8=="bidirectional_promoter_lncRNA"||$8=="macro_lncRNA"' "$OUT_DIR/genes.bed" \
  | sort -k1,1 -k2,2n > "$OUT_DIR/lncrna.bed"
awk -F'\t' '$8=="protein_coding"' "$OUT_DIR/genes.bed" | sort -k1,1 -k2,2n > "$OUT_DIR/pc.bed"

## -D ref: distance signed on genome coordinates, so "left" and "right" do not depend on strand
"$BEDTOOLS" closest -D ref -id -io -t first -a "$OUT_DIR/lncrna.bed" -b "$OUT_DIR/pc.bed" \
  > "$OUT_DIR/neighbours_left.txt"
"$BEDTOOLS" closest -D ref -iu -io -t first -a "$OUT_DIR/lncrna.bed" -b "$OUT_DIR/pc.bed" \
  > "$OUT_DIR/neighbours_right.txt"

wc -l "$OUT_DIR/lncrna.bed" "$OUT_DIR/pc.bed" "$OUT_DIR/neighbours_left.txt"
