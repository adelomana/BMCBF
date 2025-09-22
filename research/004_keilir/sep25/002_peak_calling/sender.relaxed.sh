#!/usr/bin/env bash
set -euo pipefail

# Paths
BEDDIR="/Users/adrian/research/bmcbf/004_keilir/results/002_bed"
OUTDIR="/Users/adrian/research/bmcbf/004_keilir/results/003_seacr"
SEACR="/Users/adrian/software/SEACR/SEACR_1.3.sh"
GENOME_SIZES="/Users/adrian/software/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as.genome.sizes"   
mkdir -p "$OUTDIR"

for frag in "$BEDDIR"/MITF_*FLAG_*/*.fragments.bed; do
    sample=$(basename "$(dirname "$frag")")
    ctrl=${sample/FLAG/IgG}
    ctrl_bed="$BEDDIR/$ctrl/$ctrl.fragments.bed"

    if [[ ! -f "$ctrl_bed" ]]; then
        echo "Warning: no control for $sample, skipping."
        continue
    fi

    exp_bg="$BEDDIR/$sample/${sample}.bedgraph"
    ctrl_bg="$BEDDIR/$ctrl/${ctrl}.bedgraph"

    # Run SEACR
    out_prefix="$OUTDIR/${sample}_vs_${ctrl}"
    cmd3="bash $SEACR $exp_bg $ctrl_bg norm relaxed $out_prefix"
    echo "$cmd3"; eval "$cmd3"

    out_prefix="$OUTDIR/${sample}_vs_${ctrl}_auc"
    cmd3="bash $SEACR $exp_bg 0.01 norm relaxed $out_prefix"
    echo "$cmd3"; eval "$cmd3"

    echo ""
done
