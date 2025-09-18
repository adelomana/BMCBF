#!/usr/bin/env bash
set -euo pipefail

PARENT="/Users/adrian/research/bmcbf/004_keilir/results/001_bam"

for d in "$PARENT"/MITF_A_Untreated_*_[123]; do
    sample=$(basename "$d")

    if [[ "$sample" == *IgG* ]]; then
        in="$d/${sample}.removed_dup.bam"
    else
        in="$d/${sample}.marked_dup.bam"
    fi

    out="$d/${sample}.cleaned.mapq30.bam"

    echo "Cleaning $in -> $out"

    # Get total and kept read counts
    total=$(samtools view --threads 8 -c "$in")
    kept=$(samtools view --threads 8 -c -q 30 "$in")

    # Compute percentage in bash (floating point with awk)
    perc=$(awk -v k="$kept" -v t="$total" 'BEGIN{printf "%.2f", (k/t)*100}')

    echo "Total reads: $total | Retained: $kept | ${perc}%"

    # Write cleaned BAM
    samtools view --threads 8 -b -q 30 "$in" > "$out"
done
