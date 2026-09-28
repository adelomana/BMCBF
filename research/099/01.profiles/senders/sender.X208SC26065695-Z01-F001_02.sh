#!/bin/bash
    
#SBATCH --job-name=X208SC26065695-Z01-F001_02
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=8
#SBATCH --hint=nomultithread  
#SBATCH --output=messages/messages.X208SC26065695-Z01-F001_02.out.txt   
#SBATCH --error=messages/messages.X208SC26065695-Z01-F001_02.err.txt 

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
cp -rf /hpcdata/Mimir/adrian/research/099/data/X208SC26065695-Z01-F001_02 $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/X208SC26065695-Z01-F001_02
ls
date

date

# call kallisto
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_X208SC26065695-Z01-F001_02_reverse -t 16 -b 100 --rf-stranded --verbose 
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_X208SC26065695-Z01-F001_02_forward -t 16 -b 100 --fr-stranded --verbose 
time /users/home/adrian/software/kallisto/kallisto quant -i /users/home/adrian/software/kallisto/108/index.idx -o kallisto_output_X208SC26065695-Z01-F001_02_unstranded -t 16 -b 100 --verbose 

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir /hpcdata/Mimir/adrian/research/099/results/X208SC26065695-Z01-F001_02_processed
cp -rf kallisto_output_* /hpcdata/Mimir/adrian/research/099/results/X208SC26065695-Z01-F001_02_processed/.
cp *fastp.html /hpcdata/Mimir/adrian/research/099/results/X208SC26065695-Z01-F001_02_processed/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    