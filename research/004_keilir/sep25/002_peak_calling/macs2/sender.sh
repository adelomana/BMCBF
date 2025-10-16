#!/usr/bin/env bash
set -euo pipefail

# Base directory containing your BEDPE folders
BASE="/Users/adrian/research/bmcbf/004_keilir/results/results001"
OUTDIR="/Users/adrian/research/bmcbf/004_keilir/results/macs3"

# Create the output directory if missing
mkdir -p "$OUTDIR"

# Loop over both MITF_A and MITF_M
for TF in MITF_A MITF_M; do
  # Loop over replicates 1–3
  for REP in 1 2 3; do
    FLAG_DIR="${BASE}/${TF}_Untreated_FLAG_${REP}"
    IGG_DIR="${BASE}/${TF}_Untreated_IgG_${REP}"

    FLAG_FILE="${FLAG_DIR}/human.bedpe"
    IGG_FILE="${IGG_DIR}/human.bedpe"

    if [[ ! -f "$FLAG_FILE" || ! -f "$IGG_FILE" ]]; then
      echo "Skipping ${TF} replicate ${REP} — missing files."
      continue
    fi

    OUT_SUBDIR="${OUTDIR}/${TF}_FLAG${REP}_vs_IgG${REP}"
    mkdir -p "$OUT_SUBDIR"

    CMD="macs3 callpeak -t ${FLAG_FILE} -c ${IGG_FILE} -f BEDPE -g hs -n ${TF}_FLAG${REP}_vs_IgG${REP} --outdir ${OUT_SUBDIR} --qvalue 0.01 --call-summits --cutoff-analysis"

    echo "--------------------------------------------"
    echo "Running: $CMD"
    echo "--------------------------------------------"
    time eval "$CMD"
    echo ""
  done
done
