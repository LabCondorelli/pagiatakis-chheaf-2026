#!/bin/bash
## Whole chain from the raw GEO files to the three panels, in order.
##
## Two environments are needed and they are not interchangeable: R_ENV for the R steps and
## TANGRAM_ENV for prep4, both named in 00_config.sh. Each step is run in its own, through
## `conda run`. Without conda, set RSCRIPT and PYTHON to the interpreters of the two
## environments and they are used instead. Both are checked before anything starts, so a
## missing interpreter fails in seconds rather than after the first preparation step.
##
## The preparation steps are heavy: prep1 and prep2 need about 200 GB and several hours, prep4
## up to a day when the Tangram map is refitted (CHH_TANGRAM_MAP=none). The three figure steps
## need about 50 to 100 GB.
##
## Usage: bash run_all.sh [--skip-download]
## Output: $CHH_OUT

set -eo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"      # the panel scripts and the configuration
source "${ROOT}/00_config.sh"

: "${RSCRIPT:=}"      # Rscript of R_ENV, to bypass conda
: "${PYTHON:=}"       # python3 of TANGRAM_ENV, to bypass conda

run_r ()  { if [ -n "$RSCRIPT" ]; then "$RSCRIPT" "$@"
            else conda run --no-capture-output -n "$R_ENV" Rscript "$@"; fi; }
run_py () { if [ -n "$PYTHON" ];  then "$PYTHON" "$@"
            else conda run --no-capture-output -n "$TANGRAM_ENV" python3 "$@"; fi; }

## Preflight: both interpreters must answer before the first long step.
run_r  --version >/dev/null 2>&1 || { echo "no Rscript: set RSCRIPT, or create the '$R_ENV' environment" >&2; exit 1; }
run_py --version >/dev/null 2>&1 || { echo "no python3: set PYTHON, or create the '$TANGRAM_ENV' environment" >&2; exit 1; }

[[ "$1" == "--skip-download" ]] || bash "${HERE}/00_download_geo.sh"

## prep1 and prep2 are independent; prep3 needs both, prep4 needs prep3.
run_r "${HERE}/prep1_visium_object.R"
run_r "${HERE}/prep2_snrnaseq_reference.R"
run_r "${HERE}/prep3_tangram_inputs.R"
run_py -u "${HERE}/prep4_tangram_projection.py" --dir "${CHH_WORK}/tangram" \
       --map "${CHH_TANGRAM_MAP}" --counts "${CHH_TANGRAM_COUNTS}"

## Fig. 1D and Fig. S1F read prep1 and prep4; Fig. S1E reads prep2.
run_r "${ROOT}/fig1D_spatial_chheaf_rtn4.R"
run_r "${ROOT}/figS1F_spatial_lincrna_pairs.R"
run_r "${ROOT}/figS1E_snrnaseq_chheaf_rtn4.R"
