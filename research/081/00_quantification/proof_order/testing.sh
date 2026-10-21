#!/bin/bash
    
#SBATCH --job-name=testing
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread  
#SBATCH --output=testing.out.txt   
#SBATCH --error=testing.err.txt 

cd /hpcdata/Mimir/adrian/research/081/src/msATDC5_1


echo 'about to 111222'
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/mouse/index.idx -o pepe -t 16 -b 16 --verbose \
    msATDC5_1_part1.clean.R1.fq.gz msATDC5_1_part2.clean.R1.fq.gz msATDC5_1_part3.clean.R1.fq.gz msATDC5_1_part1.clean.R2.fq.gz msATDC5_1_part2.clean.R2.fq.gz msATDC5_1_part3.clean.R2.fq.gz


echo 'about to 121212'
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/mouse/index.idx -o joan -t 16 -b 16 --verbose \
    msATDC5_1_part1.clean.R1.fq.gz msATDC5_1_part1.clean.R2.fq.gz msATDC5_1_part2.clean.R1.fq.gz msATDC5_1_part2.clean.R2.fq.gz msATDC5_1_part3.clean.R1.fq.gz msATDC5_1_part3.clean.R2.fq.gz


