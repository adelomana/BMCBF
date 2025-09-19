#!/usr/bin/env bash
set -euo pipefail

PARENT_DIR="/Users/adrian/research/bmcbf/004_keilir/results/001_bam"
PICARD="/Users/adrian/software/picard/picard.jar"
OUTDIR="/Users/adrian/research/bmcbf/004_keilir/results/002_bed"

for bam in "$PARENT_DIR"/*/*.bam; do
  sample_dir=$(dirname "$bam")
  sample=$(basename "$sample_dir")
  unit=${sample: -1}

  echo "=== Working with $sample ==="

  rg_bam="$sample_dir/${sample}.RG.bam"

  # Step 1: Add read groups
  add_rg="time java -jar $PICARD AddOrReplaceReadGroups I=$bam O=$rg_bam RGID=$sample RGLB=lib_$sample RGPL=ILLUMINA RGPU=unit$unit RGSM=$sample VERBOSITY=WARNING"
  echo "$add_rg"
  eval "$add_rg"

  # Step 2: Mark or remove duplicates
  if [[ "$sample" == *IgG* ]]; then
    dup_bam="$sample_dir/${sample}.removed_dup.bam"
    metrics="$sample_dir/${sample}.removed_dup.info.txt"
    remove_dup="true"
  else
    dup_bam="$sample_dir/${sample}.marked_dup.bam"
    metrics="$sample_dir/${sample}.marked_dup.info.txt"
    remove_dup="false"
  fi

  markdup="time java -jar $PICARD MarkDuplicates I=$rg_bam O=$dup_bam REMOVE_DUPLICATES=$remove_dup M=$metrics VERBOSITY=WARNING"
  echo "$markdup"
  eval "$markdup"

  # Step 3: Fragment size distribution
  frag_len="$sample_dir/fragment_lengths.txt"
  frag_cmd="time samtools view -F 0x04 $dup_bam | awk 'function abs(x){return (x<0?-x:x)}{if(\$9!=0) sizes[abs(\$9)]++}END{for(s in sizes) print s, sizes[s]/2}' OFS=\"\t\" | sort -n > $frag_len"
  echo "$frag_cmd"
  eval "$frag_cmd"

  # Step 4: MAPQ filtering
  out="$sample_dir/${sample}.cleaned.mapq30.bam"
  echo "Cleaning $dup_bam -> $out"
  filter_cmd="samtools view --threads 8 -b -q 30 $dup_bam > $out"
  echo "$filter_cmd"
  eval "$filter_cmd"

  frag_len="$sample_dir/fragment_lengths.mapq30.txt"
  frag_cmd="time samtools view -F 0x04 $out | awk 'function abs(x){return (x<0?-x:x)}{if(\$9!=0) sizes[abs(\$9)]++}END{for(s in sizes) print s, sizes[s]/2}' OFS=\"\t\" | sort -n > $frag_len"
  echo "$frag_cmd"
  eval "$filter_cmd"

  # Step 5: Convert to BEDPE fragments
  mkdir -p "$OUTDIR/$sample"
  qname_bam="$OUTDIR/$sample/${sample}.qnamesort.bam"
  bedpe="$OUTDIR/$sample/${sample}.bedpe"
  clean="$OUTDIR/$sample/${sample}.clean.bed"
  fragments="$OUTDIR/$sample/${sample}.fragments.bed"

  cmd1="samtools sort -n -@ 8 -o $qname_bam $out"
  cmd2="bedtools bamtobed -bedpe -i $qname_bam > $bedpe"
  cmd3="awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' $bedpe > $clean"
  cmd4="cut -f 1,2,6 $clean | sort -k1,1 -k2,2n -k3,3n > $fragments"

  echo "$cmd1"
  echo "$cmd2"
  echo "$cmd3"
  echo "$cmd4"
  eval "$cmd1"; eval "$cmd2"; eval "$cmd3"; eval "$cmd4"

  echo ""
done
