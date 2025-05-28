#!/bin/bash

#SBATCH --job-name=003
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32      
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.003.out.txt   
#SBATCH --error=messages/messages.003.err.txt 

date
pwd
hostname
conda env list

time nanovar --threads 32 /hpcdata/Mimir/adrian/research/021_arnarholt/results/aligned_bam/sample01.bam /hpcdata/Mimir/adrian/research/021_arnarholt/data/reference/release_113/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_rm.toplevel.fa /hpcdata/Mimir/adrian/research/021_arnarholt/results/variants

time nanovar --threads 32 /hpcdata/Mimir/adrian/research/021_arnarholt/results/aligned_bam/sample02.bam /hpcdata/Mimir/adrian/research/021_arnarholt/data/reference/release_113/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_rm.toplevel.fa /hpcdata/Mimir/adrian/research/021_arnarholt/results/variants

date