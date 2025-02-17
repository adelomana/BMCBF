#!/bin/bash


#
# step 01. quality control on 53
#

# this step took around 13 hours in necio5 on 53 samples
#greenPipes --inputdir /Users/adrian/research/keilir/data/raw_fastq/ --inputfile metadata.txt --libraryType pair --outputdir /Users/adrian/research/keilir/results/53/ --modes qc

#
# step 02. alignment
#
greenPipes --mode alignment --inputdir /Users/adrian/research/keilir/data/raw_fastq/ --inputfile metadata.txt --libraryType pair --outputdir /Users/adrian/research/keilir/results/53/ --refgenome /Users/adrian/research/keilir/data/osfstorage-archive/Reference_genomes/hg38_human/GRCh38.p13 --blackListedRegions /Users/adrian/research/keilir/data/osfstorage-archive/EncodeBlackListRegions/hg38-blacklist.v2.bed --spikein /Users/adrian/research/keilir/data/osfstorage-archive/Reference_genomes/ecoli/ecoli
