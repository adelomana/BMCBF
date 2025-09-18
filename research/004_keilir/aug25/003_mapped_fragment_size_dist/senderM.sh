#!/usr/bin/env bash

time samtools view -F 0x04 /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_FLAG_1/MITF_M_Untreated_FLAG_1.marked_dup.bam | awk 'function abs(x){return (x < 0 ? -x : x)}{if ($9 != 0) sizes[abs($9)]++}END{for (s in sizes) print s, sizes[s]/2}' OFS="\t" | sort -n > /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_FLAG_1/fragment_lengths.txt

time samtools view -F 0x04 /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_FLAG_2/MITF_M_Untreated_FLAG_2.marked_dup.bam | awk 'function abs(x){return (x < 0 ? -x : x)}{if ($9 != 0) sizes[abs($9)]++}END{for (s in sizes) print s, sizes[s]/2}' OFS="\t" | sort -n > /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_FLAG_2/fragment_lengths.txt

time samtools view -F 0x04 /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_FLAG_3/MITF_M_Untreated_FLAG_3.marked_dup.bam | awk 'function abs(x){return (x < 0 ? -x : x)}{if ($9 != 0) sizes[abs($9)]++}END{for (s in sizes) print s, sizes[s]/2}' OFS="\t" | sort -n > /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_FLAG_3/fragment_lengths.txt

time samtools view -F 0x04 /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_IgG_1/MITF_M_Untreated_IgG_1.removed_dup.bam | awk 'function abs(x){return (x < 0 ? -x : x)}{if ($9 != 0) sizes[abs($9)]++}END{for (s in sizes) print s, sizes[s]/2}' OFS="\t" | sort -n > /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_IgG_1/fragment_lengths.txt

time samtools view -F 0x04 /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_IgG_2/MITF_M_Untreated_IgG_2.removed_dup.bam | awk 'function abs(x){return (x < 0 ? -x : x)}{if ($9 != 0) sizes[abs($9)]++}END{for (s in sizes) print s, sizes[s]/2}' OFS="\t" | sort -n > /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_IgG_2/fragment_lengths.txt

time samtools view -F 0x04 /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_IgG_3/MITF_M_Untreated_IgG_3.removed_dup.bam | awk 'function abs(x){return (x < 0 ? -x : x)}{if ($9 != 0) sizes[abs($9)]++}END{for (s in sizes) print s, sizes[s]/2}' OFS="\t" | sort -n > /Users/adrian/research/bmcbf/004_keilir/results/001_bam/MITF_M_Untreated_IgG_3/fragment_lengths.txt
