#!/usr/bin/env bash
# Download RNA-seq delivery F26A430001919_MUStijhxR and check MD5 sums.
# Run inside tmux with the 094 conda environment active.
# Usage: bash download_rnaseq.sh /path/to/destination

set -euo pipefail

DEST="${1:-/hpcdata/Mimir/adrian/research/094/data/raw}"

aws s3 sync s3://mustijhxr-598731762349/F26A430001919_MUStijhxR/ "${DEST}" --profile seqdata

find "${DEST}" -iname '*md5*' -execdir md5sum -c {} \;
