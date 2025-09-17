#!/usr/bin/env bash

PARENT_DIR="/Users/adrian/research/bmcbf/004_keilir/results/001_bam"
PICARD="/Users/adrian/software/picard/picard.jar"

for bam in "$PARENT_DIR"/*/*.bam; do
  sample_dir=$(dirname "$bam")
  sample=$(basename "$sample_dir")          # e.g., MITF_A_Untreated_FLAG_1
  unit=${sample: -1}                        # last char (replicate number)

  rg_bam="$sample_dir/${sample}.RG.bam"

  echo "working with $sample"

  # Step 1: Add read groups
  add_rg_cmd="time java -jar $PICARD AddOrReplaceReadGroups I=$bam O=$rg_bam RGID=$sample RGLB=lib_$sample RGPL=ILLUMINA RGPU=unit$unit RGSM=$sample VERBOSITY=WARNING"
  echo "$add_rg_cmd"
  eval "$add_rg_cmd"
  echo ""

  # Step 2: Mark or remove duplicates depending on IgG/FLAG
  if [[ "$sample" == *IgG* ]]; then
    dup_bam="$sample_dir/${sample}.removed_dup.bam"
    metrics="$sample_dir/${sample}.removed_dup.info.txt"
    dup_cmd="time java -jar $PICARD MarkDuplicates I=$rg_bam O=$dup_bam REMOVE_DUPLICATES=true M=$metrics VERBOSITY=WARNING"
    echo "$dup_cmd"
    eval "$dup_cmd"
  else
    dup_bam="$sample_dir/${sample}.marked_dup.bam"
    metrics="$sample_dir/${sample}.marked_dup.info.txt"
    dup_cmd="time java -jar $PICARD MarkDuplicates I=$rg_bam O=$dup_bam REMOVE_DUPLICATES=false M=$metrics VERBOSITY=WARNING"
    echo "$dup_cmd"
    eval "$dup_cmd"
  fi
  echo ""
done