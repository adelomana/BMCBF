#!/usr/bin/env bash
set -euo pipefail
#! adrian conda activate my_env_3_10

# ==========================================================
# Aggregate FE and q-value across replicates for MITF_A IDR consensus peaks
# Author: Adrian's helper
# ==========================================================

# -------- CONFIG --------
BASE="/hpcdata/Mimir/adrian/research/004_keilir/results/idr/MITF_A"
CONS="${BASE}/MITF_A.IDR_consensus.bed"

REP1="${BASE}/rep1.idr.sorted.narrowPeak"
REP2="${BASE}/rep2.idr.sorted.narrowPeak"
REP3="${BASE}/rep3.idr.sorted.narrowPeak"

GENOME="/users/home/adrian/software/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as.genome.sizes"

# Aggregation modes
WITHIN_REP_AGG="max"       # aggregate across overlapping peaks within replicate: max or mean
ACROSS_REP_AGG="median"    # aggregate across replicates: median or mean

# Output file
OUT="${BASE}/MITF_A.IDR_consensus.with_FE_q.${WITHIN_REP_AGG}_within.${ACROSS_REP_AGG}_across.tsv"
# -------------------------

run() { echo ">>> $*"; eval "$*"; }

# ---- Sanity checks ----
for f in "$CONS" "$REP1" "$REP2" "$REP3" "$GENOME"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done

# ---- 1) Sort all inputs by genome order ----
run "bedtools sort -i \"$CONS\" -g \"$GENOME\" > tmp.cons.sorted.bed"
run "bedtools sort -i \"$REP1\" -g \"$GENOME\" > tmp.rep1.sorted.narrowPeak"
run "bedtools sort -i \"$REP2\" -g \"$GENOME\" > tmp.rep2.sorted.narrowPeak"
run "bedtools sort -i \"$REP3\" -g \"$GENOME\" > tmp.rep3.sorted.narrowPeak"

# ---- 1b) Filter replicate peaks to chromosomes that appear in the consensus ----
# This removes chrM, random contigs, alt scaffolds, etc., if they are absent from the consensus.
run "awk 'NR==FNR{c[\$1]; next} (\$1 in c)' tmp.cons.sorted.bed tmp.rep1.sorted.narrowPeak > tmp.rep1.filtered.narrowPeak"
run "awk 'NR==FNR{c[\$1]; next} (\$1 in c)' tmp.cons.sorted.bed tmp.rep2.sorted.narrowPeak > tmp.rep2.filtered.narrowPeak"
run "awk 'NR==FNR{c[\$1]; next} (\$1 in c)' tmp.cons.sorted.bed tmp.rep3.sorted.narrowPeak > tmp.rep3.filtered.narrowPeak"

run "mv tmp.rep1.filtered.narrowPeak tmp.rep1.sorted.narrowPeak"
run "mv tmp.rep2.filtered.narrowPeak tmp.rep2.sorted.narrowPeak"
run "mv tmp.rep3.filtered.narrowPeak tmp.rep3.sorted.narrowPeak"

# ---- 2) Map FE (col7) and q (col9) per replicate ----
run "cp tmp.cons.sorted.bed tmp.step0.bed"
run "bedtools map -a tmp.step0.bed -b tmp.rep1.sorted.narrowPeak -c 7,9 -o ${WITHIN_REP_AGG},${WITHIN_REP_AGG} > tmp.step1.bed"
run "bedtools map -a tmp.step1.bed -b tmp.rep2.sorted.narrowPeak -c 7,9 -o ${WITHIN_REP_AGG},${WITHIN_REP_AGG} > tmp.step2.bed"
run "bedtools map -a tmp.step2.bed -b tmp.rep3.sorted.narrowPeak -c 7,9 -o ${WITHIN_REP_AGG},${WITHIN_REP_AGG} > tmp.step3.bed"

# After mapping, the last 6 fields are: FE_r1 q_r1 FE_r2 q_r2 FE_r3 q_r3

# ---- 3) Aggregate across replicates ----
if [[ "$ACROSS_REP_AGG" == "median" ]]; then
  run "awk 'BEGIN{OFS=\"\t\"}
    function push(v, A,    n){ if(v!=\".\" && v!=\"NA\" && v!=\"\"){ n=length(A)+1; A[n]=v+0 } }
    function med3(A,    n,a,b,c,t){
      n=length(A); if(n==0)return \"NA\"; if(n==1)return A[1]; if(n==2)return (A[1]+A[2])/2;
      a=A[1]; b=A[2]; c=A[3]; if(a>b){t=a;a=b;b=t} if(b>c){t=b;b=c;c=t} if(a>b){t=a;a=b;b=t}; return b
    }
    {
      fe1=\$(NF-5); q1=\$(NF-4); fe2=\$(NF-3); q2=\$(NF-2); fe3=\$(NF-1); q3=\$NF;
      delete F; delete Q; push(fe1,F); push(fe2,F); push(fe3,F); push(q1,Q); push(q2,Q); push(q3,Q);
      print \$1,\$2,\$3, med3(F), med3(Q)
    }' tmp.step3.bed > \"$OUT\""
elif [[ "$ACROSS_REP_AGG" == "mean" ]]; then
  run "awk 'BEGIN{OFS=\"\t\"}
    function add(x,s,c,   y){ if(x!=\".\" && x!=\"NA\" && x!=\"\"){ y=x+0; s+=y; c++ } ; return s \"\\t\" c }
    {
      fe1=\$(NF-5); q1=\$(NF-4); fe2=\$(NF-3); q2=\$(NF-2); fe3=\$(NF-1); q3=\$NF;
      s=0; c=0; split(add(fe1,s,c),a,\"\\t\"); s=a[1]; c=a[2];
               split(add(fe2,s,c),a,\"\\t\"); s=a[1]; c=a[2];
               split(add(fe3,s,c),a,\"\\t\"); s=a[1]; c=a[2];
      fe_mean = (c>0 ? s/c : \"NA\");
      s=0; c=0; split(add(q1,s,c),a,\"\\t\"); s=a[1]; c=a[2];
               split(add(q2,s,c),a,\"\\t\"); s=a[1]; c=a[2];
               split(add(q3,s,c),a,\"\\t\"); s=a[1]; c=a[2];
      q_mean = (c>0 ? s/c : \"NA\");
      print \$1,\$2,\$3, fe_mean, q_mean
    }' tmp.step3.bed > \"$OUT\""
else
  echo "ERROR: ACROSS_REP_AGG must be median or mean" >&2; exit 1
fi

echo
echo "------------------------------------------"
echo "Wrote: $OUT"
echo "------------------------------------------"
echo "Temporary files kept for inspection (tmp.*.bed)"
