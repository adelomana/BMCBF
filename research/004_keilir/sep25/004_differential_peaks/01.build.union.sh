#!/bin/bash

wc -l /Users/adrian/research/bmcbf/004_keilir/results/differential_peaks/*.bed

cat /Users/adrian/research/bmcbf/004_keilir/results/differential_peaks/purple.a.bed /Users/adrian/research/bmcbf/004_keilir/results/differential_peaks/purple.m.bed \
  | sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge -d 50 > /Users/adrian/research/bmcbf/004_keilir/results/differential_peaks/union_peaks.bed
