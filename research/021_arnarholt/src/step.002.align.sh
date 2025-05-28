#!/bin/bash

#SBATCH --job-name=002
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32      
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.002.out.txt   
#SBATCH --error=messages/messages.002.err.txt 

date
pwd

time /users/home/adrian/software/dorado/dorado-0.9.1-linux-x64/bin/dorado aligner /hpcdata/Mimir/adrian/research/021_arnarholt/data/reference/release_113/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_rm.toplevel.fa.gz /hpcdata/Mimir/adrian/research/021_arnarholt/results/bam/sample01.bam --output-dir /hpcdata/Mimir/adrian/research/021_arnarholt/results/aligned_bam -t 32 -v --emit-summary

time /users/home/adrian/software/dorado/dorado-0.9.1-linux-x64/bin/dorado aligner /hpcdata/Mimir/adrian/research/021_arnarholt/data/reference/release_113/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_rm.toplevel.fa.gz /hpcdata/Mimir/adrian/research/021_arnarholt/results/bam/sample02.bam --output-dir /hpcdata/Mimir/adrian/research/021_arnarholt/results/aligned_bam -t 32 -v --emit-summary