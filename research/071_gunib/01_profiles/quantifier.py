import os, sys

def launcher(sample):

    print('building code for sample', sample)

    #
    # define fastp command
    #
    all_files = os.listdir(os.path.join(raw_fastq_dir, sample))
    working_files = sorted(f for f in all_files if f.endswith(("fq.gz", "fastq.gz")))

    r1_files = sorted(f for f in working_files if "1.f" in f)
    r2_files = sorted(f for f in working_files if "2.f" in f)

    print('these are the r1 files', r1_files)
    print('these are the r2 files', r2_files)

    fastp_commands_list = []
    clean_r1_files = []
    clean_r2_files = []

    for index, (r1, r2) in enumerate(zip(r1_files, r2_files), start=1):

        r1_path = os.path.join(raw_fastq_dir, sample, r1)
        r2_path = os.path.join(raw_fastq_dir, sample, r2)

        out_r1 = os.path.join(f"{sample}_part{index}.clean.R1.fq.gz")
        out_r2 = os.path.join(f"{sample}_part{index}.clean.R2.fq.gz")

        html_report = os.path.join(f"{sample}_part{index}.fastp.html")
        json_report = os.path.join(f"{sample}_part{index}.fastp.json")

        one_fastp_command = (
            f"time {fastp_executable} "
            f"-i {r1_path} -I {r2_path} "
            f"-o {out_r1} -O {out_r2} "
            f"--thread {number_threads} "
            f"--detect_adapter_for_pe --trim_poly_g --trim_poly_x "
            f"--cut_tail --cut_window_size 4 --cut_mean_quality 20 "
            f"--length_required 36 "
            f"-h {html_report} -j {json_report}"
        )
        fastp_commands_list.append(one_fastp_command)

        clean_r1_files.append(out_r1)
        clean_r2_files.append(out_r2)
    
    fastp_command = "\n".join(fastp_commands_list)
    print()

    #
    # define kallisto command
    #    
    kallisto_executable = "time /users/home/adrian/software/kallisto/kallisto quant"
    fastq_files = []
    for r1, r2 in zip(clean_r1_files, clean_r2_files):
        fastq_files.append(r1)
        fastq_files.append(r2)

    fastq_files_string = " ".join(fastq_files) 

    kallisto_options = f"-i {transcriptome_index} -o kallisto_output_{sample}_reverse -t {number_threads} -b 100 --rf-stranded --verbose"
    kallisto_cmd_a = f"{kallisto_executable} {kallisto_options} {fastq_files_string}"

    kallisto_options = f"-i {transcriptome_index} -o kallisto_output_{sample}_forward -t {number_threads} -b 100 --fr-stranded --verbose"
    kallisto_cmd_b = f"{kallisto_executable} {kallisto_options} {fastq_files_string}"

    kallisto_options = f"-i {transcriptome_index} -o kallisto_output_{sample}_unstranded -t {number_threads} -b 100 --verbose"
    kallisto_cmd_c = f"{kallisto_executable} {kallisto_options} {fastq_files_string}"
    
    #
    # write sender
    #
    submitter_file = 'senders/sender.{}.sh'.format(sample)

    text = f"""#!/bin/bash
    
#SBATCH --job-name={sample}
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node={int(number_threads/2)}
#SBATCH --hint=nomultithread  
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
    #sys.exit()

    return None

#
# 0. user-defined variables
#
raw_fastq_dir = '/hpcdata/Mimir/adrian/door2next/rami_data/'
results_dir = '/hpcdata/Mimir/adrian/research/071_gunib/results/'
fastp_executable = '/users/home/adrian/software/fastp/fastp'
number_threads = 16
transcriptome_index = '/users/home/adrian/software/kallisto/108/index.idx'

#
# 1. determine folders
#
dirs = [d for d in os.listdir(raw_fastq_dir) if os.path.isdir(os.path.join(raw_fastq_dir, d))]
dirs.sort()
print('detected folders...')
print(dirs)
print()

#
# 2. make some folders if they were not there yet
#
os.makedirs('senders', exist_ok=True)
os.makedirs('messages', exist_ok=True)

for sample in dirs: 
    launcher(sample)