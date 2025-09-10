import os, sys

data_dir = '/Users/adrian/research/bmcbf/004_keilir/data/'
results_dir = '/Users/adrian/research/bmcbf/004_keilir/results/000_trimmed/'

# read the files
samples = os.listdir(data_dir)
samples.sort()
print(samples)
print()

# iterate samples
for sample in samples:
    print('about to run {}'.format(sample))

    working_files_names = os.listdir(data_dir+sample)
    working_files_names.sort()

    out_dir = results_dir + sample
    if os.path.exists(out_dir) == False:
        os.mkdir(out_dir)  

    # define command
    executable = 'time fastp'


    inputs = '-i {} -I {}'.format(data_dir+sample+'/'+working_files_names[0], data_dir+sample+'/'+working_files_names[1])
    outputs = '-o {}/clean.R1.fq.gz -O {}/clean.R2.fq.gz'.format(out_dir, out_dir)
    options = '--thread 8 --trim_poly_g --detect_adapter_for_pe  --allow_gap_overlap_trimming -h {}/report.html -j {}/report.json'.format(out_dir, out_dir) 
    command = executable + ' ' + inputs + ' ' + outputs + ' ' + options

    print(command)

    print()
    os.system(command)
    print()