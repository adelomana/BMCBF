#!/usr/bin/env bash
set -euo pipefail

PARENT="/Users/adrian/research/bmcbf/004_keilir/results/001_bam"

for d in "$PARENT"/MITF_{A,M}_Untreated_*_[123]; do
    sample=$(basename "$d")

    in="$d/${sample}.cleaned.mapq30.bam"
    out="$d/fragment_lengths.cleaned.txt"

    echo "Processing $sample"

    cmd="samtools view -F 0x04 $in | awk 'function abs(x){return (x < 0 ? -x : x)}{if (\$9 != 0) sizes[abs(\$9)]++}END{for (s in sizes) print s, sizes[s]/2}' OFS=\"\t\" | sort -n > $out"

    echo "$cmd"
    eval "$cmd"
done
