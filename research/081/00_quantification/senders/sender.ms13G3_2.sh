#!/bin/bash
    
#SBATCH --job-name=ms13G3_2
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread  
#SBATCH --output=messages/messages.ms13G3_2.out.txt   
#SBATCH --error=messages/messages.ms13G3_2.err.txt 

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
cp -rf /hpcdata/Mimir/adrian/door2next/081/data/ms13G3_2 $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/ms13G3_2
ls
date
time /users/home/adrian/software/fastp/fastp -i /hpcdata/Mimir/adrian/door2next/081/data/ms13G3_2/ms13G3_2_MKRN250009502-1A_22MH3HLT4_L3_1.fq.gz -I /hpcdata/Mimir/adrian/door2next/081/data/ms13G3_2/ms13G3_2_MKRN250009502-1A_22MH3HLT4_L3_2.fq.gz -o ms13G3_2_part1.clean.R1.fq.gz -O ms13G3_2_part1.clean.R2.fq.gz --thread 16 --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 -h ms13G3_2_part1.fastp.html -j ms13G3_2_part1.fastp.json
time /users/home/adrian/software/fastp/fastp -i /hpcdata/Mimir/adrian/door2next/081/data/ms13G3_2/ms13G3_2_MKRN250009502-1A_22MH55LT4_L4_1.fq.gz -I /hpcdata/Mimir/adrian/door2next/081/data/ms13G3_2/ms13G3_2_MKRN250009502-1A_22MH55LT4_L4_2.fq.gz -o ms13G3_2_part2.clean.R1.fq.gz -O ms13G3_2_part2.clean.R2.fq.gz --thread 16 --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 -h ms13G3_2_part2.fastp.html -j ms13G3_2_part2.fastp.json
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/mouse/index.idx -o kallisto_output_ms13G3_2_reverse -t 16 -b 100 --rf-stranded --verbose ms13G3_2_part1.clean.R1.fq.gz ms13G3_2_part1.clean.R2.fq.gz ms13G3_2_part2.clean.R1.fq.gz ms13G3_2_part2.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/mouse/index.idx -o kallisto_output_ms13G3_2_forward -t 16 -b 100 --fr-stranded --verbose ms13G3_2_part1.clean.R1.fq.gz ms13G3_2_part1.clean.R2.fq.gz ms13G3_2_part2.clean.R1.fq.gz ms13G3_2_part2.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/mouse/index.idx -o kallisto_output_ms13G3_2_unstranded -t 16 -b 100 --verbose ms13G3_2_part1.clean.R1.fq.gz ms13G3_2_part1.clean.R2.fq.gz ms13G3_2_part2.clean.R1.fq.gz ms13G3_2_part2.clean.R2.fq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir /hpcdata/Mimir/adrian/research/081/results/ms13G3_2_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/081/results/ms13G3_2_processed/.
cp *fastp.html /hpcdata/Mimir/adrian/research/081/results/ms13G3_2_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    