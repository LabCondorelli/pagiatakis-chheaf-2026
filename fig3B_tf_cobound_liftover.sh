#!/bin/bash
## Fig. 3B, transcription factor co-occupancy track - cardiac TF co-bound regions of the mouse
## Rtn4 locus (Akerberg et al., fetal E12.5 and adult heart) projected onto hg38.
##
## Reciprocal liftOver: mm10 to hg38 and back, keeping only the intervals that map uniquely in
## both directions. TF code in the name: G GATA4, A MEF2A, C MEF2C, N NKX2-5, S SRF, T TBX5,
## E TEAD1.
##
## Usage: bash fig3B_tf_cobound_liftover.sh
## Output: results/fig3B_tf_cobound_hg38.bed

set -eo pipefail
source "$(dirname "$0")/00_config.sh"

: "${COBOUND_ADULT:=<PATH_TO>/adult_cobound_mm10.bed}"   # from Akerberg et al. Supplementary Data 3, mm9 to mm10
: "${COBOUND_FETAL:=<PATH_TO>/fetal_cobound_mm10.bed}"
: "${CHAIN_MM10_HG38:=<PATH_TO>/mm10ToHg38.over.chain}"
: "${CHAIN_HG38_MM10:=<PATH_TO>/hg38ToMm10.over.chain}"
: "${LIFTOVER:=liftOver}"
: "${BEDTOOLS:=bedtools}"
: "${MIN_MATCH:=0.5}"                                     # the value stated in the Methods

LOCUS="chr11:29519909-29752900"                           # Rtn4 locus, mm10
W="$OUT_DIR/fig3B"; mkdir -p "$W"

printf '%s\t%s\t%s\n' "${LOCUS%%:*}" "$(echo "${LOCUS#*:}" | cut -d- -f1)" "$(echo "${LOCUS#*:}" | cut -d- -f2)" > "$W/locus.bed"

## Restrict to the locus and label each interval by developmental stage
for s in adult fetal; do
  src=$([ "$s" = adult ] && echo "$COBOUND_ADULT" || echo "$COBOUND_FETAL")
  "$BEDTOOLS" intersect -a "$src" -b "$W/locus.bed" -wa | sort -k1,1 -k2,2n \
    | awk -v S="$s" 'BEGIN{OFS="\t"; s=toupper(substr(S,1,1)) substr(S,2)} {print $1,$2,$3,s"("$4")"}'
done | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,"id_"NR}' > "$W/merged_mm10.bed"

"$LIFTOVER" -minMatch=$MIN_MATCH "$W/merged_mm10.bed"  "$CHAIN_MM10_HG38" "$W/hg38.bed" "$W/hg38.unmapped"
"$LIFTOVER" -minMatch=$MIN_MATCH "$W/hg38.bed"         "$CHAIN_HG38_MM10" "$W/back.bed" "$W/back.unmapped"

## Keep only intervals that are unique in both directions
uniq_ids () { cut -f5 "$1" | sort | uniq -c | awk '$1==1{print $2}'; }
comm -12 <(uniq_ids "$W/hg38.bed") <(uniq_ids "$W/back.bed") > "$W/reciprocal.ids"

awk 'BEGIN{OFS="\t"} NR==FNR{k[$1]=1; next} ($5 in k){print $1,$2,$3,$4}' \
  "$W/reciprocal.ids" "$W/hg38.bed" | sort -k1,1 -k2,2n > "$OUT_DIR/fig3B_tf_cobound_hg38.bed"

wc -l "$W/merged_mm10.bed" "$OUT_DIR/fig3B_tf_cobound_hg38.bed"
