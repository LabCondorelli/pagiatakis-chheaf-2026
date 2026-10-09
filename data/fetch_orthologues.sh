#!/bin/bash
## Rebuilds data/mouse_human_orthologues.tsv from Ensembl BioMart.
##
## The frozen copy in this directory is what the GSEA of Fig. S2E used, so this script is only
## needed to refresh it against a newer Ensembl release. Doing so changes the mapping and
## therefore the enrichment: see data/README.md.
##
## www.ensembl.org refuses the BioMart endpoint, the mirrors answer, which is also why the
## biomaRt package fails here. One mouse gene can have several human homologues; the first one
## returned is kept, as in the analysis that produced the panel.
##
## Usage: bash fetch_orthologues.sh
## Output: mouse_human_orthologues.tsv

set -eo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${BIOMART:=https://useast.ensembl.org/biomart/martservice}"

Q=$(mktemp); trap 'rm -f "$Q"' EXIT
cat > "$Q" <<'XML'
<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
<Dataset name="mmusculus_gene_ensembl" interface="default">
<Filter name="with_hsapiens_homolog" excluded="0"/>
<Attribute name="ensembl_gene_id"/><Attribute name="hsapiens_homolog_associated_gene_name"/>
</Dataset></Query>
XML

curl -sL --max-time 600 --data-urlencode "query@$Q" "$BIOMART" \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "ensembl_gene_id","human_symbol"}
                $2!="" && !($1 in seen) {seen[$1]=1; print $1,$2}' \
  | { read -r h; echo "$h"; sort -k1,1; } > "$HERE/mouse_human_orthologues.tsv"

wc -l < "$HERE/mouse_human_orthologues.tsv"
