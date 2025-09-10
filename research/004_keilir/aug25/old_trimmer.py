#
# chat gpt says
#

#ILLUMINACLIP:2:30:10 is a solid default for Illumina adapters.
#SLIDINGWINDOW:4:20 gently trims poor tails without shredding short inserts.
#MINLEN:25 keeps ultra-short junk from mapping spuriously (tune to your read length).


# transcriptomics

#options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)


# CnRAP
#https://github.com/mbassalbioinformatics/CnRAP/blob/master/01_cut_n_run_pairedReads_filter_align.py
#output_command = "java -jar " +trimmomatic_path_jar+ " PE -threads " +num_cores+ " -phred33 " +read1_fq_gz+ " " +read2_fq_gz+ " " +trim_folder+ "" +sample_id+ "_1.paired_trimmomatic.fastq.gz " +unapired_trim_folder+ "" +sample_id+ "_1.unpaired.fastq.gz " +trim_folder+ "" +sample_id+ "_2.paired_trimmomatic.fastq.gz " +unapired_trim_folder+ "" +sample_id+ "_2.unpaired.fastq.gz ILLUMINACLIP:" +adapter_path+ "Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25 2> " +logs_folder+ "" +sample_id+ "_trimmomatic.err"

# https://github.com/mbassalbioinformatics/CaRAS


import os



def trim(sample):
    
    print(sample)

    executable = 'time java -jar {}trimmomatic-0.39.jar PE -threads {} -phred33 '.format(trimmomatic_path, number_threads)

    detected = os.listdir(raw_fastq_dir + sample)
    detected.sort()
    working_label = working_files[0].split('.R1')[0]
    print(detected)
    working_files = [element for element in all_files if tag in element]
    print(working_files)
    sys.exit()
    working_label = working_files[0].split('.R1')[0]

    input1 = sample + '/' + working_files[0]
    input2 = sample + '/' + working_files[1]
    
    garbage1 = 'clean_fastq/' + working_label + '_R1_garbage.fastq.gz'
    garbage2 = 'clean_fastq/' + working_label + '_R2_garbage.fastq.gz'

    input_files = input1 + ' ' + input2
    output_files = output1 + ' ' + garbage1 + ' ' + output2 + ' ' + garbage2
    
    command = executable + input_files + ' ' + output_files + options
    print(cmd)
    return None


#
# MAIN
#

# user defined variables

raw_fastq_dir = '/Users/adrian/research/bmcbf/004_keilir/data/' 
trimmomatic_path = '/Users/adrian/software/Trimmomatic-0.39/'
adapter_file = trimmomatic_path + 'adapters/TruSeq3-PE-2.fa'
number_threads = 4
options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)


# run
samples = os.listdir(raw_fastq_dir)
samples.sort()

for sample in samples:
    trim(sample)
