#!/bin/bash
    
#SBATCH --job-name=D0BioP2
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread  
#SBATCH --output=messages/messages.D0BioP2.out.txt   
#SBATCH --error=messages/messages.D0BioP2.err.txt 

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
cp -rf /hpcdata/Mimir/adrian/research/073_bologna/data/D0BioP2 $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/D0BioP2
ls
date
time /users/home/adrian/software/fastp/fastp -i /hpcdata/Mimir/adrian/research/073_bologna/data/D0BioP2/D0BioP2_1.fq.gz -I /hpcdata/Mimir/adrian/research/073_bologna/data/D0BioP2/D0BioP2_2.fq.gz -o D0BioP2_part1.clean.R1.fq.gz -O D0BioP2_part1.clean.R2.fq.gz --thread 16 --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 -h D0BioP2_part1.fastp.html -j D0BioP2_part1.fastp.json
date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_D0BioP2_reverse -t 16 -b 100 --rf-stranded --verbose D0BioP2_part1.clean.R1.fq.gz D0BioP2_part1.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_D0BioP2_forward -t 16 -b 100 --fr-stranded --verbose D0BioP2_part1.clean.R1.fq.gz D0BioP2_part1.clean.R2.fq.gz
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_D0BioP2_unstranded -t 16 -b 100 --verbose D0BioP2_part1.clean.R1.fq.gz D0BioP2_part1.clean.R2.fq.gz

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir /hpcdata/Mimir/adrian/research/073_bologna/results/D0BioP2_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/073_bologna/results/D0BioP2_processed/.
cp *fastp.html /hpcdata/Mimir/adrian/research/073_bologna/results/D0BioP2_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    