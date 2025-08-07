import os
import time
import pandas
import pysradb
import datetime
import sys


# Parameters
OUTDIR = "/hpcdata/Mimir/adrian/research/038_husavik/data/"
FASTERQ_DUMP = "/users/home/adrian/software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump"
PRE = "/users/home/adrian/software/sratoolkit.3.2.1-ubuntu64/bin/prefetch"
NUM_THREADS = 1

#
# read info
#
srr_df = pandas.read_csv('srr_list.txt', sep='\t')

#
# build senders
#
job_dir = "senders"
os.makedirs(job_dir, exist_ok=True)


for index, row in srr_df.iterrows():
    srr = row['SRRs']

    print(srr)

    # create sender
    sra_file = os.path.join('senders', f"{os.getcwd()}/{srr}/{srr}.sra")
 
    job_script = f"""#!/bin/bash

#SBATCH --job-name={srr}
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={NUM_THREADS}
#SBATCH --hint=multithread
#SBATCH --output=messages/{srr}.out.txt
#SBATCH --error=messages/{srr}.err.txt

pwd
date
time {PRE} {srr} --progress
date

date
time {FASTERQ_DUMP} "{sra_file}" --split-3 -e {NUM_THREADS} -O "{OUTDIR}" -vxp
date

echo "all done."

"""

    script_path = os.path.join(job_dir, f"job_{srr}.sh")
    with open(script_path, "w") as f:
        f.write(job_script)

    # launch sender
    os.system(f'sbatch {script_path}')
    


    

