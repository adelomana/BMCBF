import os, sys

def launcher(sample):

    print(sample)

    #
    # define java commands
    #
    all_files = os.listdir(raw_fastq_dir + sample)
    working_tags = [element.split('_R')[0] for element in all_files]
    unique_working_tags = list(set(working_tags))
    unique_working_tags.sort()

    trimmo_commands = []; trimmo_outputs = []
    for tag in unique_working_tags:

        working_files = [element for element in all_files if tag in element]

        working_label = working_files[0].split('R1')[0]

        executable='time java -jar {}trimmomatic-0.39.jar PE -threads {} -phred33 '.format(trimmomatic_path,number_threads)
        options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)

        input1 = sample + '/' + working_files[0]
        input2 = sample + '/' + working_files[1]

        output1 = 'clean_fastq/' + working_label + '_R1_clean.fastq.gz'
        output2 = 'clean_fastq/' + working_label + '_R2_clean.fastq.gz'

        trimmo_outputs.append(output1)
        trimmo_outputs.append(output2)

        garbage1 = 'clean_fastq/' + working_label + '_R1_garbage.fastq.gz'
        garbage2 = 'clean_fastq/' + working_label + '_R2_garbage.fastq.gz'

        input_files = input1 + ' ' + input2
        output_files = output1 + ' ' + garbage1 + ' ' + output2 + ' ' + garbage2

        command = executable + input_files + ' ' + output_files + options

        trimmo_commands.append(command)
    trimmo_command = "\n".join(trimmo_commands)

    #
    # define kallisto command
    #    
    executable = 'time {} quant'.format(kallisto_executable)
    fastq_files_string = ' '.join(trimmo_outputs)

    options = ' -i {} -o kallisto_output_a -t {} -b 100 --rf-stranded --verbose '.format(transcriptome_index, number_threads)
    kallisto_cmd_a = executable + options + fastq_files_string

    options = ' -i {} -o kallisto_output_b -t {} -b 100 --fr-stranded --verbose '.format(transcriptome_index, number_threads)
    kallisto_cmd_b = executable + options + fastq_files_string

    options = ' -i {} -o kallisto_output_c -t {} -b 100 --verbose '.format(transcriptome_index, number_threads)
    kallisto_cmd_c = executable + options + fastq_files_string

    #
    # write sender
    #
    submitter_file = 'senders/sender.{}.sh'.format(sample)

    text = """#!/bin/bash

#SBATCH --job-name={}
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={}       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.{}.out.txt   
#SBATCH --error=messages/messages.{}.err.txt 

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
# 3. call Trimmomatic
#
echo ""
echo "about to call trimmomatic"
date
cd $tdir
mkdir $tdir/clean_fastq
{}
date

# 4. call kallisto
{}
{}
{}

#
# 5. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
date
mkdir {}{}_processed
cp -rf kallisto_output_* {}{}_processed/.
date

#
# 6. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date

    """.format(sample, number_threads, sample, sample, raw_fastq_dir, sample, trimmo_command, kallisto_cmd_a, kallisto_cmd_b, kallisto_cmd_c, output_dir, sample, output_dir, sample)

    
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
raw_fastq_dir = '/hpcdata/Mimir/adrian/research/051/data/' 
trimmomatic_path = '/users/home/adrian/software/Trimmomatic-0.39/'
adapter_file = trimmomatic_path + 'adapters/TruSeq3-PE-2.fa'
number_threads = 16
transcriptome_index = '/users/home/adrian/software/kallisto/108/index.idx'
kallisto_executable = '/users/home/adrian/software/kallisto/kallisto'
output_dir = '/hpcdata/Mimir/adrian/research/051/results/'

#
# 1. create a directory for senders
#
if os.path.exists('senders') == False:
    os.mkdir('senders')

#
# 2. determine folders
#
all_folders = os.listdir(raw_fastq_dir)
all_folders.sort()

for sample in all_folders: 
    launcher(sample)
    #sys.exit()
