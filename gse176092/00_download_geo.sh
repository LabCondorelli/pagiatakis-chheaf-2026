#!/bin/bash
## Raw files of GEO GSE176092: 20 snRNA-seq libraries, 14 Visium sections. Rerun to resume
## if NCBI answers with a captcha.
## Usage: bash 00_download_geo.sh
## Output: ${CHH_RAW}/<sample>/

set -eo pipefail
source "$(dirname "$0")/../00_config.sh"

ACC="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
get() { curl -fsS --retry 8 --retry-delay 5 --retry-all-errors "$@"; }

# GSM id -> folder name
SAMPLES=(
  GSM5355657:GSM5355657_nonCM_Sham_1        GSM5355658:GSM5355658_nonCM_Sham_2
  GSM5355659:GSM5355659_nonCM_MI_day1       GSM5355660:GSM5355660_nonCM_MI_day1
  GSM5355661:GSM5355661_nonCM_MI_day7       GSM5355662:GSM5355662_nonCM_MI_day7
  GSM5943175:GSM5943175_CM_Sham_1           GSM5943176:GSM5943176_CM_Sham_2
  GSM5943177:GSM5943177_CM_MI_day1_IZ+BZ_1  GSM5943178:GSM5943178_CM_MI_day1_IZ+BZ_2
  GSM5943179:GSM5943179_CM_MI_day1_RZ_1     GSM5943180:GSM5943180_CM_MI_day1_RZ_2
  GSM5943181:GSM5943181_CM_MI_day7_IZ+BZ_1  GSM5943182:GSM5943182_CM_day7_IZ+BZ_2
  GSM5943183:GSM5943183_CM_MI_day7_RZ_1     GSM5943184:GSM5943184_CM_MI_day7_RZ_2
  GSM5943185:GSM5943185_CM_MI_day14_IZ+BZ_1 GSM5943186:GSM5943186_CM_MI_day14_IZ+BZ_2
  GSM5943187:GSM5943187_CM_MI_day14_RZ_1    GSM5943188:GSM5943188_CM_MI_day14_RZ_2
  GSM5355663:GSM5355663_WT_Sham             GSM5355666:GSM5355666_WT_MI_day7
  GSM5355668:GSM5355668_WT_MI_day14         GSM5943189:GSM5943189_WT_MI_day1_1
  GSM5943190:GSM5943190_WT_MI_day1_2        GSM5943191:GSM5943191_WT_MI_day1_3
  GSM5943192:GSM5943192_WT_MI_day7_2        GSM5943193:GSM5943193_WT_MI_day7_3
  GSM5943194:GSM5943194_WT_MI_day14_2       GSM5943195:GSM5943195_WT_MI_day14_3
  GSM5943196:GSM5943196_shCsrp3_MI_day1     GSM5943197:GSM5943197_shCsrp3_MI_day14
  GSM5943198:GSM5943198_Csrp30E_MI_day1     GSM5943199:GSM5943199_Csrp30E_MI_day14
)

for s in "${SAMPLES[@]}"; do
  gsm="${s%%:*}"; dir="${CHH_RAW}/${s#*:}"
  mkdir -p "${dir}"
  urls=$(get "${ACC}?acc=${gsm}&targ=self&form=text&view=brief" |
         sed -n 's/^!Sample_supplementary_file_[0-9]* = ftp:/https:/p' | tr -d '\r')
  [[ -n "${urls}" ]] || { echo "no supplementary files listed for ${gsm}" >&2; exit 1; }
  for u in ${urls}; do
    f="${u##*/}"
    case "${f}" in
      *_barcodes.tsv.gz) out=barcodes.tsv.gz ;;
      *_features.tsv.gz) out=features.tsv.gz ;;
      *_matrix.mtx.gz)   out=matrix.mtx.gz ;;
      *)                 out="${f}" ;;
    esac
    [[ -s "${dir}/${out}" ]] && continue
    get -C - -o "${dir}/${out}.part" "${u}"
    mv "${dir}/${out}.part" "${dir}/${out}"
  done
  echo "${gsm} -> ${dir} ($(ls "${dir}" | wc -l) files)"
  sleep 1
done
