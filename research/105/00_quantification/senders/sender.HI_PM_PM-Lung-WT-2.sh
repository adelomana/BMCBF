#!/bin/bash
    
#SBATCH --job-name=HI_PM_PM-Lung-WT-2
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread  
#SBATCH --output=messages/messages.HI_PM_PM-Lung-WT-2.out.txt   
#SBATCH --error=messages/messages.HI_PM_PM-Lung-WT-2.err.txt 

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
cp -rf /hpcdata/Mimir/adrian/door2next/20261118_105_parinaz_cmtr2/decode/HI_PM_PM-Lung-WT-2 $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/HI_PM_PM-Lung-WT-2
ls
date
time /users/home/adrian/software/fastp/fastp -i /hpcdata/Mimir/adrian/door2next/20261118_105_parinaz_cmtr2/decode/HI_PM_PM-Lung-WT-2/22LH5YLT4_7_GTCCTGGATA_CGGTGGCGAA_S78_L007_R1_001.fastq.gz -I /hpcdata/Mimir/adrian/door2next/20261118_105_parinaz_cmtr2/decode/HI_PM_PM-Lung-WT-2/22LH5YLT4_7_GTCCTGGATA_CGGTGGCGAA_S78_L007_R2_001.fastq.gz -o HI_PM_PM-Lung-WT-2_part1.clean.R1.fq.gz -O HI_PM_PM-Lung-WT-2_part1.clean.R2.fq.gz --thread 16 --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 -h HI_PM_PM-Lung-WT-2_part1.fastp.html -j HI_PM_PM-Lung-WT-2_part1.fastp.json
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_HI_PM_PM-Lung-WT-2_reverse -t 16 -b 100 --rf-stranded --verbose HI_PM_PM-Lung-WT-2_part1.clean.R1.fq.gz HI_PM_PM-Lung-WT-2_part1.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_HI_PM_PM-Lung-WT-2_forward -t 16 -b 100 --fr-stranded --verbose HI_PM_PM-Lung-WT-2_part1.clean.R1.fq.gz HI_PM_PM-Lung-WT-2_part1.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_HI_PM_PM-Lung-WT-2_unstranded -t 16 -b 100 --verbose HI_PM_PM-Lung-WT-2_part1.clean.R1.fq.gz HI_PM_PM-Lung-WT-2_part1.clean.R2.fq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir /hpcdata/Mimir/adrian/research/105/results/HI_PM_PM-Lung-WT-2_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/105/results/HI_PM_PM-Lung-WT-2_processed/.
cp *fastp.html /hpcdata/Mimir/adrian/research/105/results/HI_PM_PM-Lung-WT-2_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    