#!/bin/bash
# Downloads the DisGeNET gene set library used by figS2E_gsea_disgenet.R and checks that it is
# the same file the published analysis used.
#
# The file is not kept in this repository: Harmonizome is distributed under CC BY-NC-SA 4.0,
# which does not allow redistribution from an MIT licensed repository. The sha256 below pins
# the version, so the script stops when the upstream file has changed.
set -eo pipefail
cd "$(dirname "$0")"
URL=https://maayanlab.cloud/static/hdfs/harmonizome/data/disgenetdisease/gene_set_library_crisp.gmt.gz
SHA=05fcf5dcd8df7b7b79a8e16e23f50dedd67c482cb9a86a80b580faa05da06256

curl -sSL -o disgenet_harmonizome.gmt.gz "$URL"
got=$(sha256sum disgenet_harmonizome.gmt.gz | cut -d' ' -f1)
if [ "$got" = "$SHA" ]; then
  echo "  downloaded, sha256 matches the version used in the paper"
else
  echo "  WARNING: sha256 is $got, expected $SHA"
  echo "  the upstream library has changed; results will differ from the published ones"
fi
