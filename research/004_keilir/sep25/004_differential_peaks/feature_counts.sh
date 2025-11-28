featureCounts \
  -p -B -C \           # paired-end, require proper pairs, avoid chimeras
  -T 8 \               # threads
  -a union_strict_peaks.bed \
  -F SAF \             # or GTF-like; would need conversion to SAF
  -o union_counts_featureCounts.txt \
  ../bams/A1.bam ../bams/A2.bam ../bams/A3.bam \
  ../bams/B1.bam ../bams/B2.bam ../bams/B3.bam
