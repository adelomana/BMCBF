#!/bin/bash
    
#SBATCH --job-name=T251
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread  
#SBATCH --output=messages/messages.T251.out.txt   
#SBATCH --error=messages/messages.T251.err.txt 

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
cp -rf /hpcdata/Mimir/adrian/door2next/bioprint/X204SC25118826-Z01-F001/01.RawData/T251 $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/T251
ls
date
time /users/home/adrian/software/fastp/fastp -i /hpcdata/Mimir/adrian/door2next/bioprint/X204SC25118826-Z01-F001/01.RawData/T251/T251_1.fq.gz -I /hpcdata/Mimir/adrian/door2next/bioprint/X204SC25118826-Z01-F001/01.RawData/T251/T251_2.fq.gz -o T251_part1.clean.R1.fq.gz -O T251_part1.clean.R2.fq.gz --thread 16 --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 -h T251_part1.fastp.html -j T251_part1.fastp.json
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_T251_reverse -t 16 -b 100 --rf-stranded --verbose T251_part1.clean.R1.fq.gz T251_part1.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_T251_forward -t 16 -b 100 --fr-stranded --verbose T251_part1.clean.R1.fq.gz T251_part1.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_T251_unstranded -t 16 -b 100 --verbose T251_part1.clean.R1.fq.gz T251_part1.clean.R2.fq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir /hpcdata/Mimir/adrian/research/067_urmia/results/bioprint/T251_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/067_urmia/results/bioprint/T251_processed/.
cp *fastp.html /hpcdata/Mimir/adrian/research/067_urmia/results/bioprint/T251_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    