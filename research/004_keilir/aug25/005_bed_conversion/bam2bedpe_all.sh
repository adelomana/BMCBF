#!/usr/bin/env bash
set -euo pipefail

PARENT="/Users/adrian/research/bmcbf/004_keilir/results/001_bam"
OUTDIR="/Users/adrian/research/bmcbf/004_keilir/results/002_bed"

for d in "$PARENT"/MITF_*_Untreated_*_[123]; do
    sample=$(basename "$d")
    in="$d/${sample}.cleaned.mapq30.bam"

    mkdir -p "$OUTDIR/$sample"

    qname_bam="$OUTDIR/$sample/${sample}.qnamesort.bam"
    bedpe="$OUTDIR/$sample/${sample}.bedpe"
    clean="$OUTDIR/$sample/${sample}.clean.bed"
    fragments="$OUTDIR/$sample/${sample}.fragments.bed"

    echo "Processing $sample"

    # 1. Sort BAM by queryname
    cmd1="samtools sort -n -@ 8 -o $qname_bam $in"
    echo "$cmd1"
    eval "$cmd1"

    # 2. BAM to BEDPE
    cmd2="bedtools bamtobed -bedpe -i $qname_bam > $bedpe"
    echo "$cmd2"
    eval "$cmd2"

    # 3. Filter pairs: same chromosome, fragment length < 1000
    cmd3="awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' $bedpe > $clean"
    echo "$cmd3"
    eval "$cmd3"

    # 4. Collapse to fragment intervals
    cmd4="cut -f 1,2,6 $clean | sort -k1,1 -k2,2n -k3,3n > $fragments"
    echo "$cmd4"
    eval "$cmd4"

done
