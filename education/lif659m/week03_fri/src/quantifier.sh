#!/usr/bin/env bash
set -euo pipefail # a command to fail if any step in a pipe fails

#
# 0. preliminaries
#

# install kallisto
cd ..
mkdir software

wget https://github.com/pachterlab/kallisto/releases/download/v0.51.1/kallisto_linux-v0.51.1.tar.gz
tar xvf kallisto_linux-v0.51.1.tar.gz

# retrive the reference transcriptome index for zebra fish
wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/zebrafish_index_standard.tar.xz
tar xvf zebrafish_index_standard.tar.xz

# inspect the data using head. They are FASTQ files, paired FASTQ files, so two paired files describe one sample. Note the paired naming on the reads


#
# 1. quantify a transcriptome using only one core.
#
your code

#
# 2. quantify a transcriptome using two threads
#
your code

#
# 3. capture messages into a log file
#
your code
