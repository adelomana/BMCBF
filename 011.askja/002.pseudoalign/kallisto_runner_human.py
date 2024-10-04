###
### usage: time python kallisto_runner_human.py &> messages.kallisto.human.txt
###

import sys, datetime, os

def kallisto_caller(label):

    printt('about to quantify {}'.format(label))

    # detect the list of fastq files, either if we have one or multiple samples
    found_labels = os.listdir(clean_fastq_dir)
    working_labels = [found_label for found_label in found_labels if label in found_label]
    fastq_files = []
    for working_label in working_labels:
        working_files = os.listdir(clean_fastq_dir + working_label)
        wf = [element for element in working_files if 'clean' in element]
        wf.sort()
        for element in wf:
            fastq_files.append(clean_fastq_dir + working_label + '/' + element)
    fastq_files_string = ' '.join(fastq_files)

    # create the full command
    sample_output_dir = results_dir + label
    executable = 'time kallisto quant'
    options = ' -i {} -o {} -t {} -b {} {} --verbose '.format(transcriptome_index, sample_output_dir, threads, boots, strand_flag)
    command = executable + options + fastq_files_string

    print('')
    print(command)
    os.system(command)
    print('')

    return None

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

#
# 0. user defined variables
#
clean_fastq_dir = '/Users/adrian/research/011.askja/data/clean_fastq/'
boots = 100
threads = 8
results_dir = '/Users/adrian/research/011.askja/results/kallisto/kallisto.{}/'.format(boots)
transcriptome_index = '/Users/adrian/software/kallisto/human_index_standard/index.idx'

# version 5
strand_flag = '--rf-stranded'    # processed 32,576,754 reads, 28,145,049 reads pseudoaligned
#strand_flag = '--fr-stranded'   # processed 32,576,754 reads, 889,414 reads pseudoaligned
#strand_flag = ''                # processed 32,576,754 reads, 28,666,500 reads pseudoaligned

strand_flag = '--fr-stranded'
results_dir = '/Users/adrian/scratch/kallisto/'
boots = 4
threads = 4


#
# 1. recover labels
#
printt('recover labels...')

found_labels = os.listdir(clean_fastq_dir)
labels = [label for label in found_labels if 'HRH88' in label]

labels.sort()
printt(labels)

#
# 2. call kallisto quant
#
if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

for label in labels:
    kallisto_caller(label)
    #sys.exit()
