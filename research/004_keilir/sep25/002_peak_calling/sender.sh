#!/usr/bin/env bash
set -euo pipefail

# Paths
BEDDIR="/Users/adrian/research/bmcbf/004_keilir/results/000_bedgraph"
OUTDIR="/Users/adrian/research/bmcbf/004_keilir/results/002_seacr"
SEACR="/Users/adrian/software/SEACR/SEACR_1.3.sh"
mkdir -p "$OUTDIR"

for frag in "$BEDDIR"/MITF_*FLAG_*/human.bedgraph; do
    sample=$(basename "$(dirname "$frag")")
    ctrl=${sample/FLAG/IgG}

    exp_bg="$BEDDIR/$sample/human.bedgraph"
    ctrl_bg="$BEDDIR/$ctrl/human.bedgraph"

    #echo $exp_bg
    #echo $ctrl_bg
    

    if [[ ! -f "$ctrl_bg" ]]; then
        echo "Warning: no control for $sample, skipping."
        continue
    fi

    # Run SEACR
    out_prefix="$OUTDIR/${sample}_vs_${ctrl}"
    cmd3="bash $SEACR $exp_bg $ctrl_bg norm stringent $out_prefix"
    echo "$cmd3"
    eval "$cmd3"

    out_prefix="$OUTDIR/${sample}_vs_${ctrl}_auc"
    cmd3="bash $SEACR $exp_bg 0.01 norm stringent $out_prefix"
    echo "$cmd3"
    eval "$cmd3"

    out_prefix="$OUTDIR/${sample}_vs_${ctrl}"
    cmd3="bash $SEACR $exp_bg $ctrl_bg norm relaxed $out_prefix"
    echo "$cmd3"
    eval "$cmd3"

    out_prefix="$OUTDIR/${sample}_vs_${ctrl}_auc"
    cmd3="bash $SEACR $exp_bg 0.01 norm relaxed $out_prefix"
    echo "$cmd3"
    eval "$cmd3"

    echo ""
done
