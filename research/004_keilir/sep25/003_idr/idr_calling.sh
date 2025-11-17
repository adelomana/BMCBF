#!/usr/bin/env bash
set -euo pipefail

# ------------ CONFIG ------------
MACS_DIR="/hpcdata/Mimir/adrian/research/004_keilir/results/run001_macs3"               # where macs3 results live
OUT_DIR="/hpcdata/Mimir/adrian/research/004_keilir/results/idr"                         # IDR outputs
BLACKLIST="/hpcdata/Mimir/adrian/research/004_keilir/data/black/hg38-blacklist.v2.bed"  # wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
MERGE_DIST=50                                                                           # bp to cluster peaks across pairs
IDR_THRESH=0.05                                                                         # standard TF threshold
# --------------------------------

mkdir -p "$OUT_DIR"

run_cmd () {
  echo ">>> $*"
  time eval "$*"
}

process_tf () {
  local TF="$1"   # MITF_A or MITF_M
  echo "========================================"
  echo " TF: $TF"
  echo "========================================"

  # Build expected narrowPeak paths from your MACS3 naming scheme:
  # outdir: 004_macs3/${TF}_FLAG{rep}_vs_IgG{rep}/
  # file:   ${TF}_FLAG{rep}_vs_IgG{rep}_peaks.narrowPeak
  declare -a REPS=("1" "2" "3")
  declare -a NP
  for R in "${REPS[@]}"; do
    local SUBDIR="${MACS_DIR}/${TF}_FLAG${R}_vs_IgG${R}"
    local NAME="${TF}_FLAG${R}_vs_IgG${R}"
    local FILE="${SUBDIR}/${NAME}_peaks.narrowPeak"
    if [[ ! -f "$FILE" ]]; then
      echo "!! Missing narrowPeak for ${TF} rep${R}: $FILE"
      return 1
    fi
    NP+=("$FILE")
  done

  local TF_OUT="${OUT_DIR}/${TF}"
  mkdir -p "$TF_OUT"

  # Step 0: optional blacklist removal
  declare -a CLEAN
  for i in "${!NP[@]}"; do
    local base="rep$((i+1))"
    local out_clean="${TF_OUT}/${base}.clean.narrowPeak"
    if [[ -n "${BLACKLIST}" && -f "${BLACKLIST}" ]]; then
      run_cmd "bedtools intersect -v -a ${NP[$i]} -b ${BLACKLIST} > ${out_clean}"
    else
      run_cmd "cp ${NP[$i]} ${out_clean}"
    fi
    CLEAN+=("$out_clean")
  done

  # Step 1: sort for IDR by signalValue (col 7), then -log10(p) col8, -log10(q) col9
  declare -a SORTED
  for i in "${!CLEAN[@]}"; do
    local out_sorted="${TF_OUT}/rep$((i+1)).idr.sorted.narrowPeak"
    run_cmd "sort -k7,7gr -k8,8gr -k9,9gr ${CLEAN[$i]} > ${out_sorted}"
    SORTED+=("$out_sorted")
  done

  # Step 2: pairwise IDR (r1–r2, r1–r3, r2–r3)
  run_cmd "idr --samples ${SORTED[0]} ${SORTED[1]} --input-file-type narrowPeak --rank signal.value --idr-threshold ${IDR_THRESH} --output-file ${TF_OUT}/idr_r1_r2.txt --plot"
  run_cmd "idr --samples ${SORTED[0]} ${SORTED[2]} --input-file-type narrowPeak --rank signal.value --idr-threshold ${IDR_THRESH} --output-file ${TF_OUT}/idr_r1_r3.txt --plot"
  run_cmd "idr --samples ${SORTED[1]} ${SORTED[2]} --input-file-type narrowPeak --rank signal.value --idr-threshold ${IDR_THRESH} --output-file ${TF_OUT}/idr_r2_r3.txt --plot"

  # Extract passing peaks from each pair (col 12 = IDR score)
  run_cmd "awk '\$12>=1.30103{print \$1\"\\t\"\$2\"\\t\"\$3}' ${TF_OUT}/idr_r1_r2.txt > ${TF_OUT}/pass_r1_r2.bed"
  run_cmd "awk '\$12>=1.30103{print \$1\"\\t\"\$2\"\\t\"\$3}' ${TF_OUT}/idr_r1_r3.txt > ${TF_OUT}/pass_r1_r3.bed"
  run_cmd "awk '\$12>=1.30103{print \$1\"\\t\"\$2\"\\t\"\$3}' ${TF_OUT}/idr_r2_r3.txt > ${TF_OUT}/pass_r2_r3.bed"

  # Step 3: consensus = present in ≥2 of 3 pairs (merge nearby peaks within MERGE_DIST)
  run_cmd "cat ${TF_OUT}/pass_r1_r2.bed ${TF_OUT}/pass_r1_r3.bed ${TF_OUT}/pass_r2_r3.bed | sort -k1,1V -k2,2n -k3,3n | bedtools merge -d ${MERGE_DIST} -c 1 -o count > ${TF_OUT}/consensus.tmp.bed"
  run_cmd "awk '\$4>=2{print \$1\"\\t\"\$2\"\\t\"\$3}' ${TF_OUT}/consensus.tmp.bed > ${TF_OUT}/${TF}.IDR_consensus.bed"
  run_cmd "rm -f ${TF_OUT}/consensus.tmp.bed"

  echo ">> Done: ${TF_OUT}/${TF}.IDR_consensus.bed"
  echo
}

# Run for both TF flavors
process_tf "MITF_A"
process_tf "MITF_M"

echo "All done. Results in: ${OUT_DIR}"
