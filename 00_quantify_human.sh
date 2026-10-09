#!/bin/bash
## Human left ventricle quantification (MAGNet, BioProject PRJNA595151): GENCODE v47 plus the
## GENCODE CLS long-read model of syntCHHEAF, counted in a single featureCounts run.
##
## Gene level runs without -O, transcript level with -O, as in the original quantification.
## Strandedness is per sample: 139 libraries reverse stranded (-s 2), 5 unstranded (-s 0).
##
## Usage: bash 00_quantify_human.sh
## Output: results/RNA_genes_counts, results/RNA_transcripts_counts

set -eo pipefail
source "$(dirname "$0")/00_config.sh"

: "${GTF_HG38:=<PATH_TO>/gencode.v47.annotation.gtf}"   # GENCODE v47, human
: "${BAM_LIST:=<PATH_TO>/bams.txt}"                     # one BAM path per line, one directory per sample
: "${STRAND:=2}"                                        # featureCounts -s
: "${FEATURECOUNTS:=featureCounts}"                     # subread v2.0.6
: "${THREADS:=20}"

GTF="$OUT_DIR/gencode.v47.plus_CLS.gtf"

## syntCHHEAF is the CLS model anchTM_000000915718: chr2:55,112,678-55,114,107, plus strand,
## one exon of 1,430 bp.
if [ ! -s "$GTF" ]; then
  cp "$GTF_HG38" "$GTF"
  {
    printf 'chr2\tCLS\tgene\t55112678\t55114107\t.\t+\t.\tgene_id "syntCHHEAF_CLS"; gene_name "syntCHHEAF_CLS"; gene_type "lncRNA";\n'
    printf 'chr2\tCLS\ttranscript\t55112678\t55114107\t.\t+\t.\tgene_id "syntCHHEAF_CLS"; transcript_id "syntCHHEAF_CLS"; transcript_name "anchTM_000000915718"; transcript_type "lncRNA";\n'
    printf 'chr2\tCLS\texon\t55112678\t55114107\t.\t+\t.\tgene_id "syntCHHEAF_CLS"; transcript_id "syntCHHEAF_CLS"; exon_number "1";\n'
  } >> "$GTF"
fi

BAMS=$(tr '\n' ' ' < "$BAM_LIST")

"$FEATURECOUNTS" -p --countReadPairs -B -s "$STRAND" -T "$THREADS" --primary -g gene_id \
    -a "$GTF" -o "$OUT_DIR/RNA_genes_counts" $BAMS

"$FEATURECOUNTS" -O -p --countReadPairs -B -s "$STRAND" -T "$THREADS" --primary -g transcript_id \
    -a "$GTF" -o "$OUT_DIR/RNA_transcripts_counts" $BAMS

grep "^syntCHHEAF_CLS" "$OUT_DIR/RNA_transcripts_counts" | cut -c1-120
