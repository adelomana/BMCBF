#!/usr/bin/env bash
cd /hpcdata/Mimir/adrian/research/079/data
for r1 in *_1.fq.gz; do
    sample="${r1%_1.fq.gz}"        # C1_1.fq.gz -> C1
    mkdir -p "$sample"
    mv "$r1" "${sample}_2.fq.gz" "$sample/"
done