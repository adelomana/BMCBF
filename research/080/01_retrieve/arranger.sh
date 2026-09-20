#!/usr/bin/env bash
set -euo pipefail

DATADIR="/hpcdata/Mimir/adrian/research/080/data"

echo "=== Step 1: Fixing naming inconsistency ==="
for f in "$DATADIR"/MITF-X6-rep-[0-9]*; do
    [[ -e "$f" ]] || continue
    filename=$(basename "$f")
    newname="$DATADIR/${filename/rep-/rep}"
    echo "Renaming: $f -> $newname"
    mv "$f" "$newname"
done

echo ""
echo "=== Step 2: Moving pairs into folders ==="
for r1 in "$DATADIR"/*-R1.f*q.gz; do
    [[ -e "$r1" ]] || continue
    filename=$(basename "$r1")
    base="${filename%-R1.*}"
    r2="$DATADIR/${filename/R1/R2}"
    outdir="$DATADIR/$base"

    if [[ ! -f "$r2" ]]; then
        echo "WARNING: R2 pair not found for $r1 — skipping"
        continue
    fi

    echo "Creating folder: $outdir/"
    mkdir -p "$outdir"
    mv "$r1" "$r2" "$outdir/"
    echo "  Moved: $(basename $r1) + $(basename $r2)"
done

echo ""
echo "=== Done ==="