#!/bin/bash
    
#SBATCH --job-name=RPCI-WM1-1
#SBATCH --partition=mimir-interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.RPCI-WM1-1.out.txt   
#SBATCH --error=messages/messages.RPCI-WM1-1.err.txt 

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
cp -rf /hpcdata/Mimir/shared/adrian/2025.08.26_adrian2erna/rnaseq/F25A910000195_ANIytzuT/RPCI-WM1-1 $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/RPCI-WM1-1
ls
date
time /users/home/adrian/software/fastp/fastp -i RPCI-WM1-1_1.fq.gz -I RPCI-WM1-1_2.fq.gz -o clean.R1.fq.gz -O clean.R2.fq.gz --thread 16 --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 --thread 16 -h report.html -j report.json
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_a -t 16 -b 100 --rf-stranded --verbose clean.R1.fq.gz clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_b -t 16 -b 100 --fr-stranded --verbose clean.R1.fq.gz clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_c -t 16 -b 100 --verbose clean.R1.fq.gz clean.R2.fq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir /hpcdata/Mimir/shared/adrian/2025.09.12_adrian2julia/RPCI-WM1-1_processed
cp -rf kallisto_output_* /hpcdata/Mimir/shared/adrian/2025.09.12_adrian2julia/RPCI-WM1-1_processed/.
cp -rf report* /hpcdata/Mimir/shared/adrian/2025.09.12_adrian2julia/RPCI-WM1-1_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    