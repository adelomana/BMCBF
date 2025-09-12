import os

def launcher(sample):

    print(sample)

    #
    # define fastp command
    #
    output_dir = 'clean_fastq/'

    all_files = os.listdir(raw_fastq_dir + sample)
    print(all_files)
    working_files = [element for element in all_files if 'fq.gz' in element]
    working_files.sort()

    # define command
    executable = f'time {fastp_executable}'
    inputs = '-i {} -I {}'.format(working_files[0], working_files[1])
    outputs = '-o clean.R1.fq.gz -O clean.R2.fq.gz'
    options = f'--thread {number_threads} --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_tail --cut_window_size 4 --cut_mean_quality 20 --length_required 36 --thread {number_threads} -h report.html -j report.json'

    fastp_command = executable + ' ' + inputs + ' ' + outputs + ' ' + options

    #
    # define kallisto command
    #    
    executable = 'time /users/home/adrian/software/kallisto/kallisto quant'
    fastq_files_string = 'clean.R1.fq.gz clean.R2.fq.gz'

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

    text = f"""#!/bin/bash
    
#SBATCH --job-name={sample}
#SBATCH --partition=mimir-interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={number_threads}       
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.{sample}.out.txt   
#SBATCH --error=messages/messages.{sample}.err.txt 

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
cp -rf {raw_fastq_dir}{sample} $tdir/.
date

#
# 3. call fastp
#
echo ""
echo "about to call fastp"
cd $tdir/{sample}
ls
date
{fastp_command}
date

# call kallisto
{kallisto_cmd_a}
{kallisto_cmd_b}
{kallisto_cmd_c}

#
# 4. copy results back to my folders
#
echo ""
echo "about to copy results out of scratch to my dirs"
ls
date
mkdir {results_dir}{sample}_processed
cp -rf kallisto_output_* {results_dir}{sample}_processed/.
cp -rf report* {results_dir}{sample}_processed/.
date

#
# 5. clean scratch
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
raw_fastq_dir = '/hpcdata/Mimir/shared/adrian/2025.08.26_adrian2erna/rnaseq/F25A910000195_ANIytzuT/'
results_dir = '/hpcdata/Mimir/shared/adrian/2025.09.12_adrian2julia/'
fastp_executable = '/users/home/adrian/software/fastp/fastp'
number_threads = 16
transcriptome_index = '/users/home/adrian/software/kallisto/108/index.idx'

#
# 1. determine folders
#
dirs = [d for d in os.listdir(raw_fastq_dir) if os.path.isdir(os.path.join(raw_fastq_dir, d))]
dirs.sort()
print(dirs)

#
# 2. make some folders if they were not there yet
#
os.makedirs('senders', exist_ok=True)
os.makedirs('messages', exist_ok=True)

for sample in dirs: 
    launcher(sample)