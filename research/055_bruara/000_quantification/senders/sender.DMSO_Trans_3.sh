#!/bin/bash
    
#SBATCH --job-name=DMSO_Trans_3
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread  
#SBATCH --output=messages/messages.DMSO_Trans_3.out.txt   
#SBATCH --error=messages/messages.DMSO_Trans_3.err.txt 

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
cp -rf /hpcdata/Mimir/adrian/research/055_bruara/data/rnaseq/DMSO_Trans_3 $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/DMSO_Trans_3
ls
date
time /users/home/adrian/software/fastp/fastp -i DMSO_Trans_3_1.fq.gz -I DMSO_Trans_3_2.fq.gz -o clean.R1.fq.gz -O clean.R2.fq.gz --thread 16 --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 --thread 16 -h report.html -j report.json
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/zebra/index.idx -o kallisto_output_a -t 16 -b 100 --rf-stranded --verbose clean.R1.fq.gz clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/zebra/index.idx -o kallisto_output_b -t 16 -b 100 --fr-stranded --verbose clean.R1.fq.gz clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/zebra/index.idx -o kallisto_output_c -t 16 -b 100 --verbose clean.R1.fq.gz clean.R2.fq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir /hpcdata/Mimir/adrian/research/055_bruara/results/rnaseq/DMSO_Trans_3_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/055_bruara/results/rnaseq/DMSO_Trans_3_processed/.
cp -rf report* /hpcdata/Mimir/adrian/research/055_bruara/results/rnaseq/DMSO_Trans_3_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    