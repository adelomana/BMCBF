#!/bin/bash

cd /Users/adrian/research/bmcbf/004_keilir/results/differential_peaks

cat ../post_idr_peaks/tmp.a.step3.bed ../post_idr_peaks/tmp.m.step3.bed \
  | sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge -d 50 > union_peaks.bed
