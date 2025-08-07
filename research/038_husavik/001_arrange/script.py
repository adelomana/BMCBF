import os
import shutil
from glob import glob

# Set your working directory
data_dir = "/hpcdata/Mimir/adrian/research/038_husavik/data"

# Get all *_1.fastq files to extract unique SRR IDs
fastq_files = glob(os.path.join(data_dir, "SRR*_1.fastq"))

for file1 in fastq_files:
    base = os.path.basename(file1)
    srr_id = base.split("_")[0]
    file2 = file1.replace("_1.fastq", "_2.fastq")

    # Create a directory for the SRR ID
    dest_dir = os.path.join(data_dir, srr_id)
    os.makedirs(dest_dir, exist_ok=True)

    # Move both files into the directory
    for f in [file1, file2]:
        if os.path.exists(f):
            print(f"Moving {os.path.basename(f)} to {dest_dir}")
            shutil.move(f, os.path.join(dest_dir, os.path.basename(f)))
        else:
            print(f"Warning: File not found: {f}")
