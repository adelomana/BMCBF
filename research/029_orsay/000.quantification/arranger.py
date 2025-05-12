import os, shutil

dir = '/hpcdata/Mimir/adrian/research/029_orsay/data/orsay/'
files = os.listdir(dir)
print(files)

labels = []
for file in files:
    v = file.split('.')[0]

    labels.append(v)
uniquelabels = list(set(labels))
uniquelabels.sort()
print(len(uniquelabels), uniquelabels)

for uniquelabel in uniquelabels:
    new_path = dir + uniquelabel
    if os.path.exists(new_path) == False:
        os.mkdir(new_path)

    shutil.move('{}/{}.R1.fastq.gz'.format(dir, uniquelabel), "{}/".format(new_path))
    shutil.move('{}/{}.R2.fastq.gz'.format(dir, uniquelabel), "{}/".format(new_path))