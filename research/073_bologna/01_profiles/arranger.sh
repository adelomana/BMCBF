#!/usr/bin/env bash
cd /hpcdata/Mimir/adrian/research/073_bologna/data || exit 1
for r1 in *_1.fq.gz; do
    sample="${r1%_1.fq.gz}"        # D0BioP1_1.fq.gz -> D0BioP1
    mkdir -p "$sample"
    mv "$r1" "${sample}_2.fq.gz" "$sample/"
done