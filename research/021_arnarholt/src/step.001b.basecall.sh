#!/bin/bash

#SBATCH --job-name=001b
#SBATCH --partition=gpu-1xA100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32      
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.001b.out.txt   
#SBATCH --error=messages/messages.001b.err.txt 

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
# 2. copy files into scratch and cd there
#
echo ""
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/research/021_arnarholt/data/nanopore/Sheepseq_sample2/no_sample/20240808_1504_1D_PAS69274_99e60aae/pod5_pass $tdir/pod5_pass02
date

cd $tdir

#
# 3. basecalling 
#
echo ""
echo "about to basecall sup"
date
time /users/home/adrian/software/dorado/dorado-0.9.1-linux-x64/bin/dorado basecaller sup pod5_pass02 --recursive -v --device cuda:all > sample02.bam
date

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
date
cp -rf sample02.bam /hpcdata/Mimir/adrian/research/021_arnarholt/results/bam/.
date

#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    