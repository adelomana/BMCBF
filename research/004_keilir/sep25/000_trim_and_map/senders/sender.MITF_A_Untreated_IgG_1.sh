#!/bin/bash

#SBATCH --job-name=MITF_A_Untreated_IgG_1
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.MITF_A_Untreated_IgG_1.out.txt   
#SBATCH --error=messages/messages.MITF_A_Untreated_IgG_1.err.txt 

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

date
pwd
conda info --envs

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
cp -rf /hpcdata/Mimir/adrian/door/MITF_A_Untreated_IgG_1 $tdir/.
date
cd $tdir
cd MITF_A_Untreated_IgG_1

#
# 3. call fastp
#
echo ""
echo "about to call fastp"

date
time fastp -i novaseqxplus4_20240513_LH00203_0044_A22F577LT3.s_1_1711_1503.001.R1.fastq.gz -I novaseqxplus4_20240513_LH00203_0044_A22F577LT3.s_1_1711_1503.001.R2.fastq.gz -o clean.R1.fq.gz -O clean.R2.fq.gz --thread 32 --trim_poly_g --detect_adapter_for_pe --allow_gap_overlap_trimming --length_required 25 -h report.html -j report.json
date

#
# 4. call bowtie2
#
echo ""
echo "about to call bowtie2"

date
time bowtie2 -x /users/home/adrian/software/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as -1 clean.R1.fq.gz -2 clean.R2.fq.gz --end-to-end --very-sensitive --no-mixed --no-discordant  -p 32 -I 10 -X 700 | samtools view -bS -@ 32 - | samtools sort -@ 32 -o human.bam
date

#
# 5. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
mkdir /hpcdata/Mimir/adrian/research/keilir/results/MITF_A_Untreated_IgG_1

date
cp -rf report* /hpcdata/Mimir/adrian/research/keilir/results/MITF_A_Untreated_IgG_1/.
cp -rf human.bam /hpcdata/Mimir/adrian/research/keilir/results/MITF_A_Untreated_IgG_1/.
date

#
# 6. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
    