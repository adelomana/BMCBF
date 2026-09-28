#!/usr/bin/env bash
set -euo pipefail
cd /hpcdata/Mimir/adrian/research/099/data/X208SC26065695-Z01-F001_02/01.RawData

for f in *_1.fq.gz; do
  s="${f%_1.fq.gz}"                        # e.g. LumA_cntr_1
  mkdir -p "$s"
  mv -n "${s}_1.fq.gz" "${s}_2.fq.gz" "$s"/
done

echo "Sample folders: $(ls -d */ | wc -l)  |  fastq files inside: $(ls */*.fq.gz | wc -l)"