###
### usage: time python trimmer.py &> messages.trimmer.txt
###

import os, datetime, sys, re

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

def trimmomatic_caller(sample):

    printt('working with sample {}'.format(sample))

    output_dir = clean_fastq_dir + sample + '/'
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    all_files = os.listdir(raw_fastq_dir + sample)
    working_files = [element for element in all_files if 'fq.gz' in element]
    working_labels = list(set(['_'.join(element.split('_')[:-1]) for element in working_files]))

    for working_label in working_labels:

        executable='time java -jar {}trimmomatic-0.39.jar PE -threads {} -phred33 '.format(trimmomatic_path,number_threads)
        options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)

        input1 = raw_fastq_dir + sample + '/' + working_label + '_1.fq.gz'
        input2 = raw_fastq_dir + sample + '/' + working_label + '_2.fq.gz'

        output1 = output_dir + sample + '_' + working_label + '_R1_clean.fastq.gz'
        output2 = output_dir + sample + '_' + working_label + '_R2_clean.fastq.gz'

        garbage1 = output_dir + sample + '_' + working_label + '_R1_garbage.fastq.gz'
        garbage2 = output_dir + sample + '_' + working_label + '_R2_garbage.fastq.gz'

        input_files = input1 + ' ' + input2
        output_files = output1 + ' ' + garbage1 + ' ' + output2 + ' ' + garbage2

        command = executable + input_files + ' ' + output_files + options

        printt('about to clean {}'.format(working_label))
        print('')
        print(command)
        print('')
        os.system(command)
        print('')

        #sys.exit()

    return None

#
# 0. user defined variablessample
#
raw_fastq_dir = '/Users/adrian/research/egilsstadir/data/X204SC24063187-Z01-F001/01.RawData/'
clean_fastq_dir = '/Users/adrian/research/egilsstadir/data/clean_fastq/'
trimmomatic_path = '/Users/adrian/software/Trimmomatic-0.39/'
adapter_file = trimmomatic_path + 'adapters/TruSeq3-PE-2.fa'
number_threads = 4

#
# 1. recover samples
#
all_folders = os.listdir(raw_fastq_dir)
working_folders = [element for element in all_folders if element != '.DS_Store']
unique_labels = list(set(working_folders))
unique_labels.sort()
print(unique_labels)
print()

#
# 2. iterate Trimmomatic
#
for label in unique_labels:
    trimmomatic_caller(label)
