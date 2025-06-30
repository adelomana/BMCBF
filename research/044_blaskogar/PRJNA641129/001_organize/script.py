import pandas as pd
import os
import shutil

# Configuration
info_file = "/hpcdata/Mimir/adrian/research/044_vala_GSE/src/000.retrieve/info.tsv"
base_dir = "/hpcdata/Mimir/adrian/research/044_vala_GSE/data/PRJNA641129"

# Load mapping file
df = pd.read_csv(info_file, sep="\t")

# Group SRRs (run_accession) by GSM (experiment_alias)
grouped = df.groupby("experiment_alias")["run_accession"].apply(list)

# Iterate and move files
for gsm, srr_list in grouped.items():
    gsm_path = os.path.join(base_dir, gsm)
    os.makedirs(gsm_path, exist_ok=True)

    for srr in srr_list:
        srr_path = os.path.join(base_dir, srr)
        if os.path.isdir(srr_path):
            for f in os.listdir(srr_path):
                src = os.path.join(srr_path, f)
                dst = os.path.join(gsm_path, f)
                print(f"Moving {src} → {dst}")
                shutil.move(src, dst)
        else:
            print(f"SRR directory not found: {srr_path}")