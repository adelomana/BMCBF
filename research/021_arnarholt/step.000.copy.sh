#!/bin/bash

#SBATCH --job-name=arnarholt
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4       
#SBATCH --hint=multithread        
#SBATCH --output=messages.out.txt   
#SBATCH --error=messages.err.txt 

date
pwd

#
# 2. copy files 
#
echo ""
echo "about to copy files"
date
time cp -rf /hpcdata/Mimir/adrian/entrance2nextcloud/021_arnarholt_data /hpcdata/Mimir/adrian/research/021_arnarholt/.
date