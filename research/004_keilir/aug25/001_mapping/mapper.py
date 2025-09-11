import os, sys

def launcher(sample):

    print(sample)

   

    #
    # write sender
    #
    submitter_file = 'senders/sender.{}.sh'.format(sample)

    text = """#!/bin/bash

#SBATCH --job-name={}
#SBATCH --partition=mimir-interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={}       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.{}.out.txt   
#SBATCH --error=messages/messages.{}.err.txt 

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

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
cp -rf {}{} $tdir/.
date


#
# 3. call bowtie2
#
echo ""
echo "about to call bowtie2"

cd $tdir
cd {}

date
time bowtie2 -x /users/home/adrian/software/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as -1 clean.R1.fq.gz -2 clean.R2.fq.gz --end-to-end --very-sensitive --no-mixed --no-discordant --dovetail -p 32 -I 10 -X 700 | samtools view -bS -@ 32 - | samtools sort -@ 32 -o human.bam
time bowtie2 -x /users/home/adrian/software/bowtie2/e_coli/e_coli -1 clean.R1.fq.gz -2 clean.R2.fq.gz --end-to-end --very-sensitive --no-mixed --no-discordant --no-dovetail --no-overlap -p 32 -I 10 -X 700 | samtools view -bS -@ 32 - | samtools sort -@ 32 -o ecoli.bam
date

#
# 5. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
mkdir /hpcdata/Mimir/adrian/research/keilir/results/{}
date
cp -rf *.bam /hpcdata/Mimir/adrian/research/keilir/results/{}/.
date

#
# 6. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    """.format(sample, number_threads, sample, sample, fastq_dir, sample, sample, sample, sample)

    
    with open(submitter_file, 'w') as f:
        f.write(text)
   
   
    #
    # launch sender
    #
    os.system('sbatch {}'.format(submitter_file))

    return None

#
# 0. user-defined variables
#
fastq_dir = '/hpcdata/Mimir/adrian/research/keilir/data/000_trimmed/' 
number_threads = 32

#
# 1. create a directory for senders
#
if os.path.exists('senders') == False:
    os.mkdir('senders')

#
# 2. determine folders
#
all_folders = os.listdir(fastq_dir)
all_folders.sort()

for sample in all_folders: 
    launcher(sample)
    #sys.exit()