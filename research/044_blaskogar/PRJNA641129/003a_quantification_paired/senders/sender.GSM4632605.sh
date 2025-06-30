#!/bin/bash

#SBATCH --job-name=GSM4632605
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.GSM4632605.out.txt   
#SBATCH --error=messages/messages.GSM4632605.err.txt 

date
pwd

#
# 1. Create a temporary directory with a unique identifier associated with your jobid
# 
scratchlocation=/scratch/users
if [ ! -d $scratchlocation/$USER ]; then
mkdir -p $scratchlocation/$USER
fi
tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo "the scratch dir is:"
echo $tdir

#
# 2. copy files into scratch
#
echo ""
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/research/044_vala_GSE/data/PRJNA641129/GSM4632605 $tdir/.
date


#
# 3. call Trimmomatic
#
echo ""
echo "about to call trimmomatic"
date
cd $tdir
mkdir $tdir/clean_fastq
time java -jar /users/home/adrian/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 GSM4632605/PSQ5_Pre_236_S62_L003_R1_001.fastq.gz GSM4632605/PSQ5_Pre_236_S62_L003_R2_001.fastq.gz clean_fastq/PSQ5_Pre_236_S62_L003__R1_clean.fastq.gz clean_fastq/PSQ5_Pre_236_S62_L003__R1_garbage.fastq.gz clean_fastq/PSQ5_Pre_236_S62_L003__R2_clean.fastq.gz clean_fastq/PSQ5_Pre_236_S62_L003__R2_garbage.fastq.gz ILLUMINACLIP:/users/home/adrian/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_a -t 16 -b 100 --rf-stranded --verbose clean_fastq/PSQ5_Pre_236_S62_L003__R1_clean.fastq.gz clean_fastq/PSQ5_Pre_236_S62_L003__R2_clean.fastq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_b -t 16 -b 100 --fr-stranded --verbose clean_fastq/PSQ5_Pre_236_S62_L003__R1_clean.fastq.gz clean_fastq/PSQ5_Pre_236_S62_L003__R2_clean.fastq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_c -t 16 -b 100 --verbose clean_fastq/PSQ5_Pre_236_S62_L003__R1_clean.fastq.gz clean_fastq/PSQ5_Pre_236_S62_L003__R2_clean.fastq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
date
mkdir /hpcdata/Mimir/adrian/research/044_vala_GSE/results/GSM4632605_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/044_vala_GSE/results/GSM4632605_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    