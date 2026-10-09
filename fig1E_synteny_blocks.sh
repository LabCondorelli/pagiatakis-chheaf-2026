#!/bin/bash
## Fig. 1E, orange bars - human/mouse syntenic DNA blocks, i.e. the aligned blocks of the
## UCSC mouse chain alignment track (chainMm10Link) over the displayed hg38 interval.
##
## Usage: bash fig1E_synteny_blocks.sh
## Output: results/fig1E_synteny_blocks.bed

set -eo pipefail
source "$(dirname "$0")/00_config.sh"

CHROM=chr2; START=55000000; END=55220000

curl -s "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=chainMm10Link;chrom=${CHROM};start=${START};end=${END}" \
  | python3 -c '
import json, sys
d = json.load(sys.stdin)["chainMm10Link"]
for b in sorted(d, key=lambda x: x["tStart"]):
    print(b["tName"], b["tStart"], b["tEnd"], "chain%d" % b["chainId"], sep="\t")
' > "$OUT_DIR/fig1E_synteny_blocks.bed"

awk '{n++; c+=$3-$2; len[n]=$3-$2} END{
  asort(len); printf "%d blocks, %d bp aligned, median block %d bp\n", n, c, len[int(n/2)+1] }' \
  "$OUT_DIR/fig1E_synteny_blocks.bed"
