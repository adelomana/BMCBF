import os, sys

def launcher(sample):

    print(sample)


    #
    # write sender
    #
    submitter_file = 'senders/sender.{}.sh'.format(sample)

    working_files_names = os.listdir(fastq_dir + sample)
    working_files_names.sort()

    if 'IgG' in sample:
        dovetail_flag = ''
    else:
        dovetail_flag = '--dovetail'

    text = f"""#!/bin/bash

#SBATCH --job-name={sample}
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.{sample}.out.txt   
#SBATCH --error=messages/messages.{sample}.err.txt 

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
cp -rf {fastq_dir}{sample} $tdir/.
date
cd $tdir
cd {sample}

#
# 3. call fastp
#
echo ""
echo "about to call fastp"

date
time fastp -i {working_files_names[0]} -I {working_files_names[1]} -o clean.R1.fq.gz -O clean.R2.fq.gz --thread {threads} --trim_poly_g --detect_adapter_for_pe --allow_gap_overlap_trimming --length_required 25 -h report.html -j report.json
date

#
# 4. call bowtie2
#
echo ""
echo "about to call bowtie2"

date
time bowtie2 -x /users/home/adrian/software/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as -1 clean.R1.fq.gz -2 clean.R2.fq.gz --end-to-end --very-sensitive --no-mixed --no-discordant {dovetail_flag} -p {threads} -I 10 -X 700 | samtools view -bS -@ {threads} - | samtools sort -@ {threads} -o human.bam
date

#
# 5. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
mkdir /hpcdata/Mimir/adrian/research/keilir/results/{sample}

date
cp -rf report* /hpcdata/Mimir/adrian/research/keilir/results/{sample}/.
cp -rf human.bam /hpcdata/Mimir/adrian/research/keilir/results/{sample}/.
date

#
# 6. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
    """

    
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
fastq_dir = '/hpcdata/Mimir/adrian/door/' 
threads = 32

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