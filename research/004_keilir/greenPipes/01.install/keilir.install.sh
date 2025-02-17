#!/bin/bash

#
# before you start in necio5
# brew install wget
# install homer following web instructions
# brew install brewsci/bio/sambamba
# meme with macports
#

# install fast screen and tune so you have config in your bin folder


# cd /Users/adrian/hub/keilir/03_green/01.install
# rm -rf greenPipes
# rm -rf Pair-end
# conda env create --name keilir014 -f environment.yaml
# conda activate keilir014

# idr
wget https://github.com/kundajelab/idr/archive/2.0.4.zip
pip install 2.0.4.zip
idr
rm 2.0.4.zip

# final installation
rm -rf greenPipes
git clone https://github.com/snizam001/greenPipes
cd greenPipes
chmod -R +x $(pwd)/greenPipes/*
pip install $(pwd)/
cd ..
